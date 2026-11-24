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
use super::sites::{Criterion, SiteTable};
use crate::figure::{Anchor, Canvas, FAINT, INK, MUTED};

/// A site that every threshold keeps, drawn in this.
pub const KEPT: &str = "#2a78d6";
/// A site some threshold removes.
pub const REMOVED: &str = "#eb6834";
/// Panels a view shows.
pub const PANELS: usize = 4;
/// Bins per side of a panel's grid.
pub const BINS: usize = 44;
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
        let count = |x: f64| match x {
            x if x >= 1e6 => format!("{}M", x / 1e6),
            x if x >= 1e3 => format!("{}k", x / 1e3),
            x => format!("{x}"),
        };
        let inside = |p: f64| p >= r.0 - 1e-9 && p <= r.1 + 1e-9;
        match self {
            Axis::Coverage | Axis::Converted | Axis::Cells => {
                // Drawn as log10(v) for coverage, log10(v + 1) for the rest.
                let at = |v: f64| match self {
                    Axis::Coverage => v.log10(),
                    _ => (v + 1.0).log10(),
                };
                // 0, then 10 on for a count drawn as log10(v + 1): 1 would
                // crowd 0.
                let first = if self == Axis::Coverage { 0 } else { 1 };
                let decades = (first..10).map(|k| 10f64.powi(k));
                let zero = (self != Axis::Coverage).then_some(0.0);
                let all: Vec<(f64, String)> = (zero.into_iter().chain(decades))
                    .map(|v| (at(v), count(v)))
                    .filter(|(p, _)| inside(*p))
                    .collect();
                // Every other decade when there are many.
                let every = all.len().div_ceil(6).max(1);
                all.into_iter().step_by(every).collect()
            }
            _ => {
                let (values, step) = round_ticks(r);
                let places = (-step.log10().floor()).max(0.0) as usize;
                (values.into_iter())
                    .map(|v| (v, format!("{v:.places$}")))
                    .collect()
            }
        }
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
#[derive(Clone, Copy, Debug, PartialEq, serde::Serialize)]
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
    let ranks: Vec<Vec<f64>> = axes.par_iter().map(|&a| sampled_ranks(s, a)).collect();
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
                rho: pearson(&ranks[p], &ranks[q]),
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

/// Every `step`-th site's rank on axis `a` (ties at their mean rank, NaN
/// last), at most [`SAMPLE`] of them.
fn sampled_ranks(s: &Sites, a: Axis) -> Vec<f64> {
    let step = s.len().div_ceil(SAMPLE).max(1);
    let col = &s.column(a).values;
    let v: Vec<f64> = (0..s.len()).step_by(step).map(|i| col[i] as f64).collect();
    let mut order: Vec<usize> = (0..v.len()).collect();
    let key = |x: f64| if x.is_nan() { f64::INFINITY } else { x };
    order.sort_by(|&p, &q| key(v[p]).total_cmp(&key(v[q])));
    let mut rank = vec![0.0; v.len()];
    let mut i = 0;
    while i < order.len() {
        let mut j = i;
        while j + 1 < order.len() && key(v[order[j + 1]]) == key(v[order[i]]) {
            j += 1;
        }
        let mean = (i + j) as f64 / 2.0;
        order[i..=j].iter().for_each(|&o| rank[o] = mean);
        i = j + 1;
    }
    rank
}

/// The correlation of `a` and `b`; 0 when either does not vary.
fn pearson(a: &[f64], b: &[f64]) -> f64 {
    let n = a.len() as f64;
    if n < 2.0 {
        return 0.0;
    }
    let (ma, mb) = (a.iter().sum::<f64>() / n, b.iter().sum::<f64>() / n);
    let (mut ab, mut aa, mut bb) = (0.0, 0.0, 0.0);
    for (x, y) in a.iter().zip(b) {
        ab += (x - ma) * (y - mb);
        aa += (x - ma).powi(2);
        bb += (y - mb).powi(2);
    }
    if aa == 0.0 || bb == 0.0 {
        0.0
    } else {
        ab / (aa * bb).sqrt()
    }
}

/// The range axis `a` is drawn over: its finite values' 0.1% to 99.9%
/// quantiles (on a sample), widened a little; a value past it is drawn at
/// its end.
pub fn range(s: &Sites, a: Axis) -> (f64, f64) {
    s.column(a).range
}

/// [`range`] of `values` on axis `a`.
fn range_of(values: &[f32], a: Axis) -> (f64, f64) {
    let step = values.len().div_ceil(SAMPLE * 4).max(1);
    let mut v: Vec<f64> = (values.iter().step_by(step))
        .map(|&x| x as f64)
        .filter(|x| x.is_finite())
        .collect();
    if v.is_empty() {
        return (0.0, 1.0);
    }
    v.sort_by(f64::total_cmp);
    let at = |q: f64| v[((v.len() - 1) as f64 * q).round() as usize];
    let (lo, hi) = (at(0.001), at(0.999));
    let pad = ((hi - lo) * 0.02).max(1e-6);
    match a {
        Axis::EditRatio => (0.0, 1.0),
        Axis::Coverage | Axis::Converted | Axis::Cells | Axis::Pvalue => (lo.min(0.0), hi + pad),
        _ => (lo - pad, hi + pad),
    }
}

/// One panel's sites, binned: kept and removed per square (row-major from
/// the bottom left), and the sites in the box where both drawn thresholds
/// pass.
#[derive(Clone, Debug)]
pub struct Grid {
    pub x: Axis,
    pub y: Axis,
    pub xr: (f64, f64),
    pub yr: (f64, f64),
    pub bins: usize,
    pub kept: Vec<u32>,
    pub removed: Vec<u32>,
    /// The box, on the drawn scales (an open end at the range's).
    pub xbox: (f64, f64),
    pub ybox: (f64, f64),
    /// Sites in the box kept, and removed by a threshold not drawn here.
    pub in_kept: usize,
    pub in_removed: usize,
    /// The busiest square's sites, which the deepest colour stands for.
    pub top: u32,
}

/// The panel of axes `x` and `y` over `s` at the thresholds `f`, binned in
/// parallel from the axes' cached values; the box's counts read off the
/// fail bits (in the box: passing both axes' thresholds).
pub fn grid(s: &Sites, (x, y): (Axis, Axis), f: &SiteFilterArgs, bins: usize) -> Grid {
    use rayon::prelude::*;
    let (cx, cy) = (s.column(x), s.column(y));
    let (xr, yr) = (cx.range, cy.range);
    let open = |(lo, hi): (Option<f64>, Option<f64>), r: (f64, f64)| {
        (
            lo.unwrap_or(f64::NEG_INFINITY),
            hi.unwrap_or(f64::INFINITY),
            r,
        )
    };
    let (xlo, xhi, _) = open(x.passes(f, s), xr);
    let (ylo, yhi, _) = open(y.passes(f, s), yr);
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
    let top = (kept.iter().zip(&removed))
        .map(|(k, r)| k + r)
        .max()
        .unwrap_or(0)
        .max(1);
    Grid {
        x,
        y,
        xr,
        yr,
        bins,
        top,
        kept,
        removed,
        xbox: (xlo.max(xr.0), xhi.min(xr.1)),
        ybox: (ylo.max(yr.0), yhi.min(yr.1)),
        in_kept,
        in_removed,
    }
}

impl Grid {
    /// The box's counts, in a line.
    pub fn in_box(&self) -> String {
        format!(
            "in the box: {} kept · {} removed by a threshold not drawn here",
            crate::site_analysis::pileup::fmt_thousands(self.in_kept as i64),
            crate::site_analysis::pileup::fmt_thousands(self.in_removed as i64)
        )
    }

    /// [`Self::in_box`] for a narrow panel.
    pub fn in_box_short(&self) -> String {
        format!(
            "in the box: {} kept · {} removed by others",
            crate::site_analysis::pileup::fmt_thousands(self.in_kept as i64),
            crate::site_analysis::pileup::fmt_thousands(self.in_removed as i64)
        )
    }

    /// Square `at`'s colour, as `(r, g, b)` over a white ground: kept or
    /// removed by which most of its sites are, deeper as it holds more of
    /// the busiest square's; none when empty.
    pub fn colour(&self, at: usize) -> Option<(u8, u8, u8)> {
        let (k, r) = (self.kept[at], self.removed[at]);
        let total = k + r;
        if total == 0 {
            return None;
        }
        let depth = 0.18 + 0.82 * (total as f64).ln_1p() / (self.top as f64).ln_1p();
        let base = rgb(if k >= r { KEPT } else { REMOVED });
        let mix = |c: u8| (255.0 + (c as f64 - 255.0) * depth).round() as u8;
        Some((mix(base.0), mix(base.1), mix(base.2)))
    }

    /// Drawn in the box `(x, y, w, h)` of `canvas`: the squares, the box
    /// where both thresholds pass dashed, the axes and their titles, and
    /// `title` over it.
    pub fn draw(&self, canvas: &mut Canvas, (x, y, w, h): (f64, f64, f64, f64), title: &str) {
        let size = (w.min(h) / 26.0).clamp(7.0, 10.0);
        let (l, r, t, b) = (size * 4.6, size, size * 2.0, size * 3.6);
        let (pw, ph) = ((w - l - r).max(1.0), (h - t - b).max(1.0));
        let (px, py) = (x + l, y + t);
        let sx = |v: f64| px + (v - self.xr.0) / (self.xr.1 - self.xr.0) * pw;
        let sy = |v: f64| py + ph - (v - self.yr.0) / (self.yr.1 - self.yr.0) * ph;
        if !title.is_empty() {
            canvas.text(x + l, y + size * 1.2, title, size, Anchor::Start, INK);
        }
        let (cw, ch) = (pw / self.bins as f64, ph / self.bins as f64);
        for by in 0..self.bins {
            for bx in 0..self.bins {
                if let Some((r, g, b)) = self.colour(by * self.bins + bx) {
                    let fill = format!("#{r:02x}{g:02x}{b:02x}");
                    let (cx, cy) = (px + bx as f64 * cw, py + ph - (by + 1) as f64 * ch);
                    canvas.rect(cx, cy, cw + 0.05, ch + 0.05, &fill);
                }
            }
        }
        // The box where both drawn thresholds pass.
        let (x0, x1) = (sx(self.xbox.0), sx(self.xbox.1));
        let (y0, y1) = (sy(self.ybox.1), sy(self.ybox.0));
        if x1 > x0 && y1 > y0 {
            for (a, b2, c, d) in [
                (x0, y0, x1, y0),
                (x1, y0, x1, y1),
                (x1, y1, x0, y1),
                (x0, y1, x0, y0),
            ] {
                canvas.dashed(a, b2, c, d, INK, 1.0);
            }
        }
        for (a, b2, c, d) in [(px, py + ph, px + pw, py + ph), (px, py, px, py + ph)] {
            canvas.line(a, b2, c, d, FAINT, 0.8);
        }
        for (v, label) in self.x.ticks(self.xr) {
            canvas.line(sx(v), py + ph, sx(v), py + ph + 3.0, FAINT, 0.8);
            let at = (sx(v), py + ph + size * 1.4);
            canvas.text(at.0, at.1, &label, size * 0.9, Anchor::Middle, MUTED);
        }
        for (v, label) in self.y.ticks(self.yr) {
            canvas.line(px - 3.0, sy(v), px, sy(v), FAINT, 0.8);
            let at = (px - size * 0.5, sy(v) + size * 0.35);
            canvas.text(at.0, at.1, &label, size * 0.9, Anchor::End, MUTED);
        }
        canvas.text(
            px + pw / 2.0,
            y + h - size * 0.6,
            self.x.name(),
            size,
            Anchor::Middle,
            INK,
        );
        canvas.vtext(x + size * 0.9, py + ph / 2.0, self.y.name(), size, INK);
    }
}

/// About five round ticks within `(lo, hi)`, and their step.
fn round_ticks((lo, hi): (f64, f64)) -> (Vec<f64>, f64) {
    let span = hi - lo;
    if span.is_nan() || span <= 0.0 {
        return (vec![lo], 1.0);
    }
    let raw = span / 4.0;
    let mag = 10f64.powf(raw.log10().floor());
    let step = [1.0, 2.0, 5.0, 10.0]
        .iter()
        .map(|m| m * mag)
        .find(|&s| s >= raw)
        .unwrap_or(10.0 * mag);
    let first = (lo / step).ceil() * step;
    let values = (0..)
        .map(|k| first + k as f64 * step)
        .take_while(|&v| v <= hi + 1e-9)
        // -0 written as 0.
        .map(|v| if v.abs() < step * 1e-9 { 0.0 } else { v })
        .collect();
    (values, step)
}

/// `#rrggbb` as its parts.
fn rgb(hex: &str) -> (u8, u8, u8) {
    let p = |i: usize| u8::from_str_radix(&hex[i..i + 2], 16).unwrap_or(0);
    (p(1), p(3), p(5))
}

/// The panels of `pairs` (with why each was picked) over `s` at `f`, as one
/// figure of `w` by `h`: two by two, the legend over them.
pub fn figure(
    s: &Sites,
    pairs: &[(Axis, Axis, String)],
    f: &SiteFilterArgs,
    (w, h): (f64, f64),
    title: &str,
) -> String {
    crate::figure::svg(w, h, |c| {
        c.text(12.0, 16.0, title, 11.0, Anchor::Start, INK);
        c.rect(12.0, 24.0, 9.0, 9.0, KEPT);
        c.text(
            25.0,
            32.0,
            "kept by every threshold",
            9.0,
            Anchor::Start,
            MUTED,
        );
        c.rect(150.0, 24.0, 9.0, 9.0, REMOVED);
        c.text(
            163.0,
            32.0,
            "removed by at least one",
            9.0,
            Anchor::Start,
            MUTED,
        );
        c.text(
            300.0,
            32.0,
            "darker: more sites · dashed: where both drawn thresholds pass",
            9.0,
            Anchor::Start,
            MUTED,
        );
        let (top, cols) = (40.0, 2usize);
        let rows = pairs.len().div_ceil(cols).max(1);
        let (pw, ph) = ((w - 12.0) / cols as f64, (h - top) / rows as f64);
        for (k, (x, y, why)) in pairs.iter().enumerate() {
            let g = grid(s, (*x, *y), f, BINS);
            let (cx, cy) = (6.0 + (k % cols) as f64 * pw, top + (k / cols) as f64 * ph);
            g.draw(
                c,
                (cx, cy, pw - 6.0, ph - 22.0),
                &format!("{} by {}", y.short(), x.short()),
            );
            c.text(
                cx + 6.0,
                cy + ph - 12.0,
                &g.in_box(),
                8.0,
                Anchor::Start,
                MUTED,
            );
            c.text(cx + 6.0, cy + ph - 3.0, why, 7.5, Anchor::Start, MUTED);
        }
    })
}
