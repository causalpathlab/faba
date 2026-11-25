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
pub use super::density::BINS;
use super::density::{Plot, Squares};
use super::scale::{along, edge, spread, Ruler, Scale, Scales};
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
    /// Per site, a bit ([`Criterion::bit`]) for every threshold it fails.
    pub fails: &'a [u8],
    /// The criteria that apply.
    pub criteria: &'a [Criterion],
    /// Each axis' values, read once (they never change with a threshold).
    pub cache: &'a Cache,
}

impl<'a> Sites<'a> {
    fn len(&self) -> usize {
        self.fails.len()
    }

    fn cells(&self, i: usize) -> Option<usize> {
        self.n_cells.map(|c| c[i])
    }

    /// Axis `a`'s values per site, as shown and along its line (on the
    /// scale the cache holds for it).
    pub fn column(&self, a: Axis) -> &'a Column {
        self.cache.columns[a as usize].get_or_init(|| {
            use rayon::prelude::*;
            let shown: Vec<f32> = (0..self.len())
                .into_par_iter()
                .map(|i| a.value(self, i) as f32)
                .collect();
            let scale = self.cache.scales.of(a);
            let values = (shown.par_iter())
                .map(|&v| along(scale, v as f64) as f32)
                .collect();
            Column {
                spread: spread(&shown),
                shown,
                values,
            }
        })
    }

    /// Axis `a`'s ruler at the thresholds `f`, reaching those that apply
    /// and are on: what the panels `faba site-plot` draws read it on.
    pub fn ruler(&self, a: Axis, f: &SiteFilterArgs) -> Ruler {
        let on = |c: Criterion| self.criteria.contains(&c) && !c.is_off(f);
        let scale = self.cache.scales.of(a);
        Ruler::cut((a, scale), self.column(a).spread, f, on, self.cache.bins).0
    }
}

/// One axis' values per site: as shown (see [`Axis::shown`]) and along its
/// line (see [`super::scale::along`]); and where the shown ones spread,
/// read once (see [`spread`]).
pub struct Column {
    pub shown: Vec<f32>,
    pub values: Vec<f32>,
    pub spread: Option<(f64, f64)>,
}

/// Each axis' [`Column`], filled as it is first asked for, on `scales`;
/// its rulers in about `bins` bins.
pub struct Cache {
    columns: [std::sync::OnceLock<Column>; 7],
    pub scales: Scales,
    pub bins: usize,
}

impl Default for Cache {
    /// On the default scales, in a histogram's own bins.
    fn default() -> Self {
        let width = super::site_plot::HIST_OWN.0;
        Self::new(Scales::default(), data_beans::figure::bins_for(width))
    }
}

impl Cache {
    /// Empty, each axis to be read on its scale in `scales`, its ruler in
    /// about `bins` bins.
    pub fn new(scales: Scales, bins: usize) -> Self {
        Self {
            columns: Default::default(),
            scales,
            bins,
        }
    }

    /// The columns read so far, by `Axis as usize`.
    pub fn into_columns(self) -> [Option<Column>; 7] {
        self.columns.map(std::sync::OnceLock::into_inner)
    }
}

/// What a panel's axis shows: a criterion's raw value as shown (a p-value
/// as `-log10 p`), on its quantity's [`Scale`]. The edit ratio carries both
/// its bounds. Every site figure takes these as its quantities.
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

    /// The quantity criterion `c` thresholds.
    pub fn of(c: Criterion) -> Self {
        (Axis::ALL.into_iter())
            .find(|a| a.criteria().contains(&c))
            .expect("every criterion thresholds a quantity")
    }

    /// Its thresholds' bits in a site's failed checks.
    pub fn bits(self) -> u8 {
        self.criteria().iter().fold(0, |b, &c| b | c.bit())
    }

    /// What it is, as an axis title.
    pub fn name(self) -> &'static str {
        match self {
            Axis::Coverage => "coverage (reads)",
            Axis::Converted => "converted reads",
            Axis::EditRatio => "edit ratio",
            Axis::Fold => "fold over control",
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

    /// Its name on the command line and in file names.
    pub fn key(self) -> &'static str {
        match self {
            Axis::Coverage => "coverage",
            Axis::Converted => "converted",
            Axis::EditRatio => "edit_ratio",
            Axis::Fold => "fold",
            Axis::LogOdds => "log_odds",
            Axis::Pvalue => "pv",
            Axis::Cells => "cells",
        }
    }

    /// The quantity named `s`: its [`Axis::key`], with `-` for `_`, or
    /// `pvalue`.
    pub fn parse(s: &str) -> anyhow::Result<Self> {
        let s = s.trim().to_ascii_lowercase().replace('-', "_");
        let s = if s == "pvalue" { "pv" } else { s.as_str() };
        Axis::ALL.into_iter().find(|a| a.key() == s).ok_or_else(|| {
            let known: Vec<&str> = Axis::ALL.iter().map(|a| a.key()).collect();
            anyhow::anyhow!("no quantity {s}: one of {}", known.join(", "))
        })
    }

    /// Its scale unless `--scale` says otherwise: `-log10 p` linear (it is
    /// already a log, and the strong sites spread out on it); fold on a log,
    /// as a ratio reads (1 near the start, the long tail folded in); counts
    /// on a log, their tails long; the edit ratio (0 to 1) and the log odds
    /// (signed, so it has no log) linear.
    pub fn default_scale(self) -> Scale {
        match self {
            Axis::Coverage | Axis::Converted | Axis::Cells | Axis::Fold => Scale::Log,
            Axis::EditRatio | Axis::LogOdds | Axis::Pvalue => Scale::Linear,
        }
    }

    /// Whether a log scale can show it: not the signed log odds.
    pub fn can_log(self) -> bool {
        self != Axis::LogOdds
    }

    /// Whether it counts reads or cells.
    pub fn integer(self) -> bool {
        matches!(self, Axis::Coverage | Axis::Converted | Axis::Cells)
    }

    /// Its title on a paper's panel: short.
    pub fn title(self) -> &'static str {
        match self {
            Axis::Coverage => "coverage",
            Axis::Converted => "converted reads",
            Axis::EditRatio => "edit ratio",
            Axis::Fold => "fold over control",
            Axis::LogOdds => "log odds",
            Axis::Pvalue => "-log10 p",
            Axis::Cells => "cells",
        }
    }

    /// `raw` (its criterion's raw value) as shown: a p-value as `-log10 p`
    /// (infinite at 0, NaN where it is, so the site is left out), the rest
    /// as they are.
    pub fn shown(self, raw: f64) -> f64 {
        match self {
            Axis::Pvalue if raw <= 0.0 => f64::INFINITY,
            Axis::Pvalue if raw.is_nan() => raw,
            Axis::Pvalue => (-raw.log10()).max(0.0),
            _ => raw,
        }
    }

    /// Site `i`'s value as shown: infinite at a clean control, NaN where
    /// there is none.
    pub fn value(self, s: &Sites, i: usize) -> f64 {
        self.shown(self.criteria()[0].raw(s.table, i, s.cells(i)))
    }

    /// Whether the sites have it: fold and log odds need a control arm,
    /// cells a `_site` matrix.
    pub fn applies(self, s: &Sites) -> bool {
        self.criteria().iter().any(|c| s.criteria.contains(c))
    }

    /// Where its thresholds pass, along its line, each at its edge (see
    /// [`super::scale::edge`]): a lower and an upper end, none for an end
    /// no threshold holds (an off one, or none).
    pub fn passes(self, f: &SiteFilterArgs, s: &Sites) -> (Option<f64>, Option<f64>) {
        let (mut lo, mut hi) = (None, None);
        let scale = s.cache.scales.of(self);
        for &c in self.criteria() {
            if !s.criteria.contains(&c) || c.is_off(f) {
                continue;
            }
            let v = along(scale, edge(self, scale, c, c.get(f)));
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

/// One panel's sites, binned: kept and removed per square, and the sites
/// in the box where both drawn thresholds pass.
#[derive(Clone, Debug)]
pub struct Grid {
    pub x: Axis,
    pub y: Axis,
    pub xr: (f64, f64),
    pub yr: (f64, f64),
    pub squares: Squares,
    /// Each axis' ticks along its line, from its ruler.
    pub ticks: [Vec<(f64, String)>; 2],
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
    let (rx, ry) = (s.ruler(x, f), s.ruler(y, f));
    let (xr, yr) = (rx.span(), ry.span());
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
        ticks: [rx.ticks(), ry.ticks()],
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
        let plot = Plot::on((x, y + head, w, h - head), self.xr, self.yr);
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
        plot.axes(
            canvas,
            [&self.ticks[0], &self.ticks[1]],
            [self.x.name(), self.y.name()],
        );
    }
}

/// Auto's `k` pairs over `s` as panels: the axes and why each was picked.
pub fn auto_panels(s: &Sites, k: usize) -> Vec<(Axis, Axis, String)> {
    (auto_pairs(s, k).into_iter())
        .map(|p| (p.x, p.y, p.why()))
        .collect()
}

/// The panels of `pairs` (with why each was picked) over `s` at `f`, as one
/// figure of `w` by `h` in `style`: two by two, the legend over them, the
/// panels binned in parallel. Alone (`--panel-only`), the panels' squares,
/// boxes, axes and axis titles only.
pub fn figure(
    s: &Sites,
    pairs: &[(Axis, Axis, String)],
    f: &SiteFilterArgs,
    (w, h): (f64, f64),
    title: &str,
    style: &crate::figure::style::Style,
) -> String {
    use rayon::prelude::*;
    let grids: Vec<Grid> = (pairs.par_iter())
        .map(|(x, y, _)| grid(s, (*x, *y), f, BINS))
        .collect();
    let mut canvas = Canvas::styled(w, h, style);
    // Laid out at 8 points on the page the font scales.
    let (w, h) = canvas.page();
    let c = &mut canvas;
    let alone = style.panel_only;
    if !alone {
        c.text(12.0, 16.0, title, 11.0, Anchor::Start, INK);
        let roles = style.palette.roles();
        let key = [
            (12.0, roles.kept, "kept by every threshold"),
            (150.0, roles.removed, "removed by at least one"),
        ];
        for (x, colour, said) in key {
            c.rect(x, 24.0, 9.0, 9.0, colour);
            c.text(x + 13.0, 32.0, said, 9.0, Anchor::Start, MUTED);
        }
        let note = "darker: more sites · dashed: where both drawn thresholds pass";
        c.text(300.0, 32.0, note, 9.0, Anchor::Start, MUTED);
    }
    // Room over the panels for the key, and under each for its notes.
    let (top, notes) = if alone { (0.0, 0.0) } else { (40.0, 22.0) };
    let cols = 2usize;
    let rows = pairs.len().div_ceil(cols).max(1);
    let (pw, ph) = ((w - 12.0) / cols as f64, (h - top) / rows as f64);
    for (k, ((x, y, why), g)) in pairs.iter().zip(&grids).enumerate() {
        let (cx, cy) = (6.0 + (k % cols) as f64 * pw, top + (k / cols) as f64 * ph);
        if alone {
            g.draw(c, (cx, cy, pw - 6.0, ph), "");
            continue;
        }
        let title = format!("{} by {}", y.short(), x.short());
        g.draw(c, (cx, cy, pw - 6.0, ph - notes), &title);
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
    canvas.finish()
}
