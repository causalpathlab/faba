//! The site thresholds two at a time: each panel puts two criteria's values
//! against each other, the sites kept by every threshold apart from those
//! removed by any, with the box where both drawn thresholds pass. A site in
//! the box but removed is removed by a threshold not drawn there, which no
//! histogram of one criterion shows.
//!
//! The pairs shown first are chosen by [`auto_pairs`]: those whose two
//! thresholds split the most sites (pass one, fail the other) on axes that
//! agree least, spread over the criteria.

use super::args::SiteFilterArgs;
use super::density::{Plot, Squares};
pub use super::density::{BINS, KEPT, REMOVED};
use super::sites::{Criterion, SiteTable};
use crate::figure::{Anchor, Canvas, INK, MUTED};

/// Panels a view shows.
pub const PANELS: usize = 4;
/// What an axis already shown is worth, a use each, when pairs are picked.
const REUSE: f64 = 0.15;
/// Sites at most the rank correlation is read on.
const SAMPLE: usize = 5000;

/// The sites of one modality as the view holds them.
#[derive(Clone, Copy)]
pub struct Sites<'a> {
    pub table: &'a SiteTable,
    /// Kept cells per site, when a `_site` matrix gave them.
    pub n_cells: Option<&'a [usize]>,
    /// Per site, a bit (`1 << c as u8`) for every threshold it fails.
    pub fails: &'a [u8],
    /// The criteria that apply.
    pub criteria: &'a [Criterion],
    /// Each axis' values, read once (they never change with a threshold).
    pub cache: &'a Cache,
}

impl Sites<'_> {
    fn len(&self) -> usize {
        self.fails.len()
    }

    fn cells(&self, i: usize) -> Option<usize> {
        self.n_cells.map(|c| c[i])
    }

    /// Axis `a`'s values per site, on its drawn scale, and its range.
    pub fn column(&self, a: Axis) -> &Column {
        self.cache.0[a as usize].get_or_init(|| {
            use rayon::prelude::*;
            let values: Vec<f32> = (0..self.len())
                .into_par_iter()
                .map(|i| a.value(self, i) as f32)
                .collect();
            let range = range_of(&values, a);
            Column { values, range }
        })
    }
}

/// One axis' values per site, on its drawn scale, and the range it is
/// drawn over.
pub struct Column {
    pub values: Vec<f32>,
    pub range: (f64, f64),
}

/// Each axis' [`Column`], filled as it is first asked for.
#[derive(Default)]
pub struct Cache([std::sync::OnceLock<Column>; 7]);

/// What a panel's axis shows: a criterion's raw value, drawn on a scale
/// that spreads it (logs for counts, `-log10` for p-values). The edit ratio
/// carries both its bounds.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash, serde::Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Axis {
    Coverage,
    Converted,
    EditRatio,
    Fold,
    LogOdds,
    Pvalue,
    Cells,
}

impl Axis {
    pub const ALL: [Axis; 7] = [
        Axis::Coverage,
        Axis::Converted,
        Axis::EditRatio,
        Axis::Fold,
        Axis::LogOdds,
        Axis::Pvalue,
        Axis::Cells,
    ];

    /// The thresholds on this axis.
    pub fn criteria(self) -> &'static [Criterion] {
        match self {
            Axis::Coverage => &[Criterion::MinCoverage],
            Axis::Converted => &[Criterion::MinConverted],
            Axis::EditRatio => &[Criterion::MinEditRatio, Criterion::MaxEditRatio],
            Axis::Fold => &[Criterion::MinFold],
            Axis::LogOdds => &[Criterion::MinLogOdds],
            Axis::Pvalue => &[Criterion::MaxPv],
            Axis::Cells => &[Criterion::MinCells],
        }
    }

    /// Its thresholds' bits in a site's failed checks.
    pub fn bits(self) -> u8 {
        self.criteria().iter().fold(0, |b, &c| b | 1 << c as u8)
    }

    /// What it is, as an axis title.
    pub fn name(self) -> &'static str {
        match self {
            Axis::Coverage => "coverage (reads)",
            Axis::Converted => "converted reads",
            Axis::EditRatio => "edit ratio",
            Axis::Fold => "fold over control (log2)",
            Axis::LogOdds => "log odds",
            Axis::Pvalue => "-log10 p",
            Axis::Cells => "cells with a converted read",
        }
    }

    /// In a word or two, for a narrow title.
    pub fn short(self) -> &'static str {
        match self {
            Axis::Coverage => "coverage",
            Axis::Converted => "converted",
            Axis::EditRatio => "edit ratio",
            Axis::Fold => "fold",
            Axis::LogOdds => "log odds",
            Axis::Pvalue => "-log10 p",
            Axis::Cells => "cells",
        }
    }

    /// `raw` (its criterion's raw value) on the drawn scale.
    pub fn scaled(self, raw: f64) -> f64 {
        match self {
            Axis::Coverage => raw.max(1.0).log10(),
            Axis::Converted | Axis::Cells => (raw + 1.0).log10(),
            Axis::EditRatio | Axis::LogOdds => raw,
            Axis::Fold => raw.log2(),
            // A p-value of 0 is past any scale: infinite, drawn at its end.
            Axis::Pvalue if raw <= 0.0 => f64::INFINITY,
            Axis::Pvalue => -raw.log10(),
        }
    }

    /// Its ticks within the drawn range `r`, each a drawn position and its
    /// label: on a count's log scale at 0 (where the scale has it), 1, 10,
    /// 100 and on, as counts read; otherwise about five round values, with
    /// the decimals their step needs.
    pub fn ticks(self, r: (f64, f64)) -> Vec<(f64, String)> {
        let inside = |p: f64| p >= r.0 - 1e-9 && p <= r.1 + 1e-9;
        let decades = |from: i32| (from..10).map(|k| 10f64.powi(k));
        // The values ticked, and where each is drawn: log10(v) for
        // coverage; log10(v + 1) for the rest, from 0 then 10 (1 would
        // crowd 0).
        let (values, at): (Vec<f64>, fn(f64) -> f64) = match self {
            Axis::Coverage => (decades(0).collect(), f64::log10),
            Axis::Converted | Axis::Cells => {
                (std::iter::once(0.0).chain(decades(1)).collect(), |v| {
                    (v + 1.0).log10()
                })
            }
            _ => {
                let (values, step) = round_ticks(r);
                let places = (-step.log10().floor()).max(0.0) as usize;
                return (values.into_iter())
                    .map(|v| (v, format!("{v:.places$}")))
                    .collect();
            }
        };
        let all: Vec<(f64, String)> = (values.into_iter())
            .map(|v| (at(v), crate::tui::ui::compact(v)))
            .filter(|(p, _)| inside(*p))
            .collect();
        // Every other decade when there are many.
        let every = all.len().div_ceil(6).max(1);
        all.into_iter().step_by(every).collect()
    }

    /// Site `i`'s value on the drawn scale: infinite at a clean control,
    /// NaN where there is none.
    pub fn value(self, s: &Sites, i: usize) -> f64 {
        self.scaled(self.criteria()[0].raw(s.table, i, s.cells(i)))
    }

    /// Whether the sites have it: fold and log odds need a control arm,
    /// cells a `_site` matrix.
    pub fn applies(self, s: &Sites) -> bool {
        self.criteria().iter().any(|c| s.criteria.contains(c))
    }

    /// Where its thresholds pass, on the drawn scale: a lower and an upper
    /// end, none for an end no threshold holds (an off one, or none).
    pub fn passes(self, f: &SiteFilterArgs, s: &Sites) -> (Option<f64>, Option<f64>) {
        let (mut lo, mut hi) = (None, None);
        for &c in self.criteria() {
            if !s.criteria.contains(&c) || c.is_off(f) {
                continue;
            }
            let v = self.scaled(c.get(f));
            // A p-value's upper bound is a lower one on -log10.
            match (c.is_max(), self == Axis::Pvalue) {
                (true, false) => hi = Some(v),
                _ => lo = Some(v),
            }
        }
        (lo, hi)
    }
}

/// The axes `s` has, in [`Axis::ALL`]'s order.
pub fn axes(s: &Sites) -> Vec<Axis> {
    Axis::ALL.into_iter().filter(|a| a.applies(s)).collect()
}

/// A pair of axes and how much it shows: the sites failing exactly one of
/// its thresholds, and the axes' rank correlation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct Pair {
    pub x: Axis,
    pub y: Axis,
    pub split: usize,
    pub rho: f64,
}

impl Pair {
    /// The sites it splits, as many as the axes are apart.
    pub fn score(&self) -> f64 {
        self.split as f64 * (1.0 - self.rho.abs())
    }

    /// Why it was picked, in a line.
    pub fn why(&self) -> String {
        format!(
            "auto: {} sites pass one of these thresholds and fail the other; rank correlation {:.2}",
            crate::site_analysis::pileup::fmt_thousands(self.split as i64),
            self.rho
        )
    }
}

/// Every pair of `s`'s axes, scored: one pass over the sites counts each
/// set of failed thresholds, and each pair's split is read off those.
pub fn score_pairs(s: &Sites) -> Vec<Pair> {
    use rayon::prelude::*;
    let axes = axes(s);
    let sampled: Vec<Vec<f32>> = axes.par_iter().map(|&a| sampled(s, a)).collect();
    let masks = (s.fails.par_chunks(1 << 16))
        .map(|chunk| {
            let mut n = [0usize; 256];
            chunk.iter().for_each(|&m| n[m as usize] += 1);
            n
        })
        .reduce(
            || [0; 256],
            |mut a, b| {
                a.iter_mut().zip(b).for_each(|(x, y)| *x += y);
                a
            },
        );
    let mut out = Vec::new();
    for p in 0..axes.len() {
        for q in p + 1..axes.len() {
            let (bx, by) = (axes[p].bits(), axes[q].bits());
            let split = (0..256)
                .filter(|&m| (m as u8 & bx != 0) != (m as u8 & by != 0))
                .map(|m| masks[m])
                .sum();
            out.push(Pair {
                x: axes[p],
                y: axes[q],
                split,
                rho: rank_correlation(&sampled[p], &sampled[q]),
            });
        }
    }
    out
}

/// `k` pairs (fewer when there are fewer), the best first, an axis already
/// shown counting for [`REUSE`] a use, so the panels spread over the axes.
pub fn auto_pairs(s: &Sites, k: usize) -> Vec<Pair> {
    let mut left = score_pairs(s);
    let mut used = rustc_hash::FxHashMap::<Axis, i32>::default();
    let mut out = Vec::new();
    while out.len() < k && !left.is_empty() {
        let worth = |p: &Pair| {
            let n = used.get(&p.x).unwrap_or(&0) + used.get(&p.y).unwrap_or(&0);
            p.score() * REUSE.powi(n)
        };
        let best = (0..left.len())
            .max_by(|&a, &b| worth(&left[a]).total_cmp(&worth(&left[b])))
            .expect("some pair is left");
        let p = left.swap_remove(best);
        *used.entry(p.x).or_default() += 1;
        *used.entry(p.y).or_default() += 1;
        out.push(p);
    }
    out
}

/// Every `step`-th site's value on axis `a`, at most [`SAMPLE`] of them,
/// NaN ranked last (as +inf).
fn sampled(s: &Sites, a: Axis) -> Vec<f32> {
    let step = s.len().div_ceil(SAMPLE).max(1);
    let col = &s.column(a).values;
    (col.iter().step_by(step))
        .map(|&v| if v.is_nan() { f32::INFINITY } else { v })
        .collect()
}

/// The rank correlation of `a` and `b`; 0 when it has none (an axis that
/// does not vary, too few sites).
fn rank_correlation(a: &[f32], b: &[f32]) -> f64 {
    let rho = legume_numeric::matrix::agreement::spearman(a, b) as f64;
    if rho.is_nan() {
        0.0
    } else {
        rho
    }
}

/// The range axis `a` is drawn over: its finite values' 0.1% to 99.9%
/// quantiles (on a sample), widened a little; a value past it is drawn at
/// its end.
fn range_of(values: &[f32], a: Axis) -> (f64, f64) {
    let step = values.len().div_ceil(SAMPLE * 4).max(1);
    let v: Vec<f32> = (values.iter().step_by(step))
        .copied()
        .filter(|x| x.is_finite())
        .collect();
    if v.is_empty() {
        return (0.0, 1.0);
    }
    let q = legume_numeric::matrix::utils::quantiles(&v, &[0.001, 0.999]);
    let (lo, hi) = (q[0] as f64, q[1] as f64);
    let pad = ((hi - lo) * 0.02).max(1e-6);
    match a {
        Axis::EditRatio => (0.0, 1.0),
        Axis::Coverage | Axis::Converted | Axis::Cells | Axis::Pvalue => (lo.min(0.0), hi + pad),
        _ => (lo - pad, hi + pad),
    }
}

/// One panel's sites, binned: kept and removed per square, and the sites
/// in the box where both drawn thresholds pass.
#[derive(Clone, Debug)]
pub struct Grid {
    pub x: Axis,
    pub y: Axis,
    pub xr: (f64, f64),
    pub yr: (f64, f64),
    pub squares: Squares,
    /// The box, on the drawn scales (an open end at the range's).
    pub xbox: (f64, f64),
    pub ybox: (f64, f64),
    /// Sites in the box kept, and removed by a threshold not drawn here.
    pub in_kept: usize,
    pub in_removed: usize,
}

/// The panel of axes `x` and `y` over `s` at the thresholds `f`, binned in
/// parallel from the axes' cached values; the box's counts read off the
/// fail bits (in the box: passing both axes' thresholds).
pub fn grid(s: &Sites, (x, y): (Axis, Axis), f: &SiteFilterArgs, bins: usize) -> Grid {
    use rayon::prelude::*;
    let (cx, cy) = (s.column(x), s.column(y));
    let (xr, yr) = (cx.range, cy.range);
    let open = |(lo, hi): (Option<f64>, Option<f64>)| {
        (lo.unwrap_or(f64::NEG_INFINITY), hi.unwrap_or(f64::INFINITY))
    };
    let ((xlo, xhi), (ylo, yhi)) = (open(x.passes(f, s)), open(y.passes(f, s)));
    // Each value's square along an axis, from a scale and offset.
    let along = |(lo, hi): (f64, f64)| {
        let k = (bins as f64 / (hi - lo)) as f32;
        let last = bins as f32 - 1.0;
        move |v: f32| ((v - lo as f32) * k).floor().clamp(0.0, last) as usize
    };
    let (sx, sy) = (along(xr), along(yr));
    let drawn = x.bits() | y.bits();
    let n = bins * bins;
    // Per chunk: kept and removed per square, then the box's two counts.
    let (kept, removed, in_kept, in_removed) = (0..s.len())
        .into_par_iter()
        .with_min_len(1 << 15)
        .fold(
            || (vec![0u32; n], vec![0u32; n], 0usize, 0usize),
            |mut acc, i| {
                let (vx, vy) = (cx.values[i], cy.values[i]);
                if vx.is_nan() || vy.is_nan() {
                    return acc;
                }
                let m = s.fails[i];
                let at = sy(vy) * bins + sx(vx);
                if m == 0 {
                    acc.0[at] += 1;
                    acc.2 += 1;
                } else {
                    acc.1[at] += 1;
                    acc.3 += usize::from(m & drawn == 0);
                }
                acc
            },
        )
        .reduce(
            || (vec![0u32; n], vec![0u32; n], 0, 0),
            |mut a, b| {
                a.0.iter_mut().zip(&b.0).for_each(|(p, q)| *p += q);
                a.1.iter_mut().zip(&b.1).for_each(|(p, q)| *p += q);
                (a.0, a.1, a.2 + b.2, a.3 + b.3)
            },
        );
    Grid {
        x,
        y,
        xr,
        yr,
        squares: Squares::new(bins, kept, removed),
        xbox: (xlo.max(xr.0), xhi.min(xr.1)),
        ybox: (ylo.max(yr.0), yhi.min(yr.1)),
        in_kept,
        in_removed,
    }
}

impl Grid {
    /// The box's counts, in a line: in full, or for a `narrow` panel.
    pub fn in_box(&self, narrow: bool) -> String {
        let n = |v: usize| crate::site_analysis::pileup::fmt_thousands(v as i64);
        let removed = if narrow {
            "removed by others"
        } else {
            "removed by a threshold not drawn here"
        };
        format!(
            "in the box: {} kept · {} {removed}",
            n(self.in_kept),
            n(self.in_removed)
        )
    }

    /// Drawn in the box `(x, y, w, h)` of `canvas`: the squares, the box
    /// where both thresholds pass dashed, the axes and their titles, and
    /// `title` over it (none when empty).
    pub fn draw(&self, canvas: &mut Canvas, (x, y, w, h): (f64, f64, f64, f64), title: &str) {
        let head = if title.is_empty() { 0.0 } else { 14.0 };
        let plot = Plot::new((x, y + head, w, h - head), self.xr, self.yr);
        if !title.is_empty() {
            canvas.text(plot.px, y + 10.0, title, 10.0, Anchor::Start, INK);
        }
        self.squares.draw(canvas, &plot);
        // The box where both drawn thresholds pass.
        let (x0, y1) = plot.at(self.xbox.0, self.ybox.0);
        let (x1, y0) = plot.at(self.xbox.1, self.ybox.1);
        if x1 > x0 && y1 > y0 {
            for (a, b, c, d) in [
                (x0, y0, x1, y0),
                (x1, y0, x1, y1),
                (x1, y1, x0, y1),
                (x0, y1, x0, y0),
            ] {
                canvas.dashed(a, b, c, d, INK, 1.0);
            }
        }
        let ticks = [self.x.ticks(self.xr), self.y.ticks(self.yr)];
        plot.axes(
            canvas,
            [&ticks[0], &ticks[1]],
            [self.x.name(), self.y.name()],
        );
    }
}

/// About five round ticks within `r` (legume-tui's), and their step.
fn round_ticks((lo, hi): (f64, f64)) -> (Vec<f64>, f64) {
    let values = legume_tui::hist::Axis::linear(lo, hi).ticks();
    let step = match values.as_slice() {
        [a, b, ..] => b - a,
        _ => 1.0,
    };
    // -0 written as 0.
    let values = (values.into_iter())
        .map(|v| if v.abs() < step * 1e-9 { 0.0 } else { v })
        .collect();
    (values, step)
}

/// Auto's `k` pairs over `s` as panels: the axes and why each was picked.
pub fn auto_panels(s: &Sites, k: usize) -> Vec<(Axis, Axis, String)> {
    (auto_pairs(s, k).into_iter())
        .map(|p| (p.x, p.y, p.why()))
        .collect()
}

/// The panels of `pairs` (with why each was picked) over `s` at `f`, as one
/// figure of `w` by `h`: two by two, the legend over them, the panels
/// binned in parallel.
pub fn figure(
    s: &Sites,
    pairs: &[(Axis, Axis, String)],
    f: &SiteFilterArgs,
    (w, h): (f64, f64),
    title: &str,
) -> String {
    use rayon::prelude::*;
    let grids: Vec<Grid> = (pairs.par_iter())
        .map(|(x, y, _)| grid(s, (*x, *y), f, BINS))
        .collect();
    crate::figure::svg(w, h, |c| {
        c.text(12.0, 16.0, title, 11.0, Anchor::Start, INK);
        let key = [
            (12.0, KEPT, "kept by every threshold"),
            (150.0, REMOVED, "removed by at least one"),
        ];
        for (x, colour, said) in key {
            c.rect(x, 24.0, 9.0, 9.0, colour);
            c.text(x + 13.0, 32.0, said, 9.0, Anchor::Start, MUTED);
        }
        let note = "darker: more sites · dashed: where both drawn thresholds pass";
        c.text(300.0, 32.0, note, 9.0, Anchor::Start, MUTED);
        let (top, cols) = (40.0, 2usize);
        let rows = pairs.len().div_ceil(cols).max(1);
        let (pw, ph) = ((w - 12.0) / cols as f64, (h - top) / rows as f64);
        for (k, ((x, y, why), g)) in pairs.iter().zip(&grids).enumerate() {
            let (cx, cy) = (6.0 + (k % cols) as f64 * pw, top + (k / cols) as f64 * ph);
            let title = format!("{} by {}", y.short(), x.short());
            g.draw(c, (cx, cy, pw - 6.0, ph - 22.0), &title);
            c.text(
                cx + 6.0,
                cy + ph - 12.0,
                &g.in_box(false),
                8.0,
                Anchor::Start,
                MUTED,
            );
            c.text(cx + 6.0, cy + ph - 3.0, why, 7.5, Anchor::Start, MUTED);
        }
    })
}
