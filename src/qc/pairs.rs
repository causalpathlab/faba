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
use super::arms::{self, ArmColumn, ArmLook, Front};
pub use super::density::BINS;
use super::density::{Plot, Squares};
use super::scale::{along, edge, spread, Ruler, Scale, Scales};
use super::sites::{Criterion, SiteTable};
use crate::figure::{mixed, rgb, Anchor, Canvas, INK, MUTED};
use data_beans::figure::panel::Marginals;

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

    pub(crate) fn cells(&self, i: usize) -> Option<usize> {
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
        self.ruler_over(a, f, self.column(a).spread)
    }

    /// [`Sites::ruler`] over values spread as `spread`.
    fn ruler_over(&self, a: Axis, f: &SiteFilterArgs, spread: Option<(f64, f64)>) -> Ruler {
        let on = |c: Criterion| self.criteria.contains(&c) && !c.is_off(f);
        let scale = self.cache.scales.of(a);
        Ruler::cut((a, scale), spread, f, on, self.cache.bins).0
    }

    /// Axis `a` on each arm (see [`arms::per_arm`]), read once; none where
    /// it has no value on the control arm.
    pub fn arm_column(&self, a: Axis) -> Option<&'a ArmColumn> {
        (arms::per_arm(a, self.table))
            .then(|| self.cache.arms[a as usize].get_or_init(|| ArmColumn::of(self, a)))
    }

    /// Axis `a`'s values on the signal arm along its line, as the view
    /// style's 2D panels place them: its own reads where it has a value on
    /// each arm (see [`arms::wt`]), else its column's.
    pub fn wt_values(&self, a: Axis) -> &'a [f32] {
        match self.arm_column(a) {
            Some(arm) => &arm.wt,
            None => &self.column(a).values,
        }
    }

    /// Axis `a`'s ruler as the view style's 2D panels read it: over its
    /// column and, where it has a value on each arm, both arms' values too
    /// (the coverage's column, signal and control added, spans either
    /// arm's own already).
    pub fn arm_ruler(&self, a: Axis, f: &SiteFilterArgs) -> Ruler {
        let own = self.column(a).spread;
        let spread = match (own, self.arm_column(a).and_then(|c| c.spread)) {
            (Some((l1, h1)), Some((l2, h2))) => Some((l1.min(l2), h1.max(h2))),
            (own, arms) => own.or(arms),
        };
        self.ruler_over(a, f, spread)
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
    /// Each axis on each arm, where it has a value on both.
    arms: [std::sync::OnceLock<ArmColumn>; 7],
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
            arms: Default::default(),
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

/// One panel's sites, binned: the signal arm's kept and removed per
/// square, the control arm's where both axes have a value on it, and the
/// sites in the box where both drawn thresholds pass.
#[derive(Clone, Debug)]
pub struct Grid {
    pub x: Axis,
    pub y: Axis,
    pub xr: (f64, f64),
    pub yr: (f64, f64),
    /// The signal arm (WT): kept by every threshold, and the rest.
    pub squares: Squares,
    /// The control arm's (MUT) sites per square, as `squares` is laid out;
    /// none unless both axes have a value on it (see [`arms::per_arm`]).
    pub mutant: Option<Vec<u32>>,
    /// Whether each axis reads each arm's own reads (see [`arms::wt`]).
    pub per_arm: [bool; 2],
    /// Each axis' ticks along its line, from its ruler.
    pub ticks: [Vec<(f64, String)>; 2],
    /// Each tick's value as shown, for labelling it another way (a
    /// p-value as p).
    pub tick_values: [Vec<f64>; 2],
    /// The box, on the drawn scales (an open end at the range's).
    pub xbox: (f64, f64),
    pub ybox: (f64, f64),
    /// Sites in the box kept, and removed by a threshold not drawn here.
    pub in_kept: usize,
    pub in_removed: usize,
}

/// The panel of axes `x` and `y` over `s` at the thresholds `f`, binned in
/// parallel: the signal arm's values (its own reads where an axis has a
/// value on each arm, see [`Sites::wt_values`]) kept apart from removed,
/// and the control arm's where both axes have one there; the box's counts
/// read off the fail bits (in the box: passing both axes' thresholds).
pub fn grid(s: &Sites, (x, y): (Axis, Axis), f: &SiteFilterArgs, bins: usize) -> Grid {
    use rayon::prelude::*;
    let (cx, cy) = (s.wt_values(x), s.wt_values(y));
    let (rx, ry) = (s.arm_ruler(x, f), s.arm_ruler(y, f));
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
                let (vx, vy) = (cx[i], cy[i]);
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
    // The control arm, where both axes have a value on it.
    let mutant = match (s.arm_column(x), s.arm_column(y)) {
        (Some(ax), Some(ay)) => Some(
            (0..s.len())
                .into_par_iter()
                .with_min_len(1 << 15)
                .fold(
                    || vec![0u32; n],
                    |mut acc, i| {
                        let (vx, vy) = (ax.mutant[i], ay.mutant[i]);
                        if !(vx.is_nan() || vy.is_nan()) {
                            acc[sy(vy) * bins + sx(vx)] += 1;
                        }
                        acc
                    },
                )
                .reduce(
                    || vec![0u32; n],
                    |mut a, b| {
                        a.iter_mut().zip(&b).for_each(|(p, q)| *p += q);
                        a
                    },
                ),
        ),
        _ => None,
    };
    let values = |r: &Ruler| r.axis.ticks();
    Grid {
        x,
        y,
        xr,
        yr,
        squares: Squares::new(bins, kept, removed),
        mutant,
        per_arm: [x, y].map(|a| s.arm_column(a).is_some()),
        ticks: [rx.ticks(), ry.ticks()],
        tick_values: [values(&rx), values(&ry)],
        xbox: (xlo.max(xr.0), xhi.min(xr.1)),
        ybox: (ylo.max(yr.0), yhi.min(yr.1)),
        in_kept,
        in_removed,
    }
}

/// How a square of one arm is drawn: filled, or white and outlined in
/// this colour.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Tone {
    Fill((u8, u8, u8)),
    Outlined((u8, u8, u8)),
}

/// `#rrggbb` as its parts.
fn rgb_of(hex: &str) -> (u8, u8, u8) {
    let (r, g, b) = rgb(hex);
    (r as u8, g as u8, b as u8)
}

impl Tone {
    /// As `#rrggbb`.
    fn hex((r, g, b): (u8, u8, u8)) -> String {
        format!("#{r:02x}{g:02x}{b:02x}")
    }
}

/// How the view style draws a 2D panel beyond its squares: the arms (which
/// in front, how the p-value reads) and the axis titles asked for.
#[derive(Clone, Copy, Default)]
pub struct PairLook<'a> {
    pub arms: ArmLook,
    pub titles: Option<&'a super::site_plot::Titles>,
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

    /// The busiest square of either arm, at least 1.
    fn top(&self) -> u32 {
        let mutant = self.mutant.iter().flatten().copied().max().unwrap_or(0);
        self.squares.top.max(mutant).max(1)
    }

    /// Square `at` of each arm as drawn, none where it is empty: on the
    /// signal arm, in `palette`'s kept colour where most of it is kept,
    /// else white outlined, deeper as it holds more of the busiest square
    /// of either arm; on the control arm, light gray.
    pub fn tones(&self, at: usize, palette: crate::figure::style::Palette) -> [Option<Tone>; 2] {
        let top = (self.top() as f64).ln_1p();
        // The signal arm's from well past the control's light gray, so the
        // sparsest kept square still reads dark.
        let depth = |n: u32| 0.6 + 0.4 * (n as f64).ln_1p() / top;
        let white = (255.0, 255.0, 255.0);
        let (k, o) = (self.squares.kept[at], self.squares.other[at]);
        let wt = (k + o > 0).then(|| {
            let d = depth(k + o);
            match k >= o {
                true => Tone::Fill(mixed(white, rgb(palette.roles().kept), d)),
                false => Tone::Outlined(mixed(white, rgb(arms::OUTLINE), d)),
            }
        });
        let m = self.mutant.as_ref().map_or(0, |v| v[at]);
        let mutant = (m > 0).then(|| Tone::Fill(rgb_of(arms::MUT_FILL)));
        [wt, mutant]
    }

    /// Square `at` as one tone, the arm in front first: what a terminal
    /// cell shows of it.
    pub fn front_tone(
        &self,
        at: usize,
        palette: crate::figure::style::Palette,
        front: Front,
    ) -> Option<Tone> {
        let [wt, mutant] = self.tones(at, palette);
        match front {
            Front::Wt => wt.or(mutant),
            Front::Mut => mutant.or(wt),
        }
    }

    /// The squares over `plot`, in the canvas' palette: the arm behind
    /// first, then the one in front over it.
    fn draw_squares(&self, c: &mut Canvas, plot: &Plot, front: Front) {
        let palette = c.style().palette;
        let n = self.squares.bins;
        let (cw, ch) = (plot.pw / n as f64, plot.ph / n as f64);
        let order = match front {
            Front::Wt => [1, 0],
            Front::Mut => [0, 1],
        };
        for arm in order {
            for by in 0..n {
                for bx in 0..n {
                    let Some(tone) = self.tones(by * n + bx, palette)[arm] else {
                        continue;
                    };
                    let (x, y) = (
                        plot.px + bx as f64 * cw,
                        plot.py + plot.ph - (by + 1) as f64 * ch,
                    );
                    match tone {
                        Tone::Fill(rgb) => c.rect(x, y, cw + 0.05, ch + 0.05, &Tone::hex(rgb)),
                        Tone::Outlined(rgb) => {
                            let w = (0.12 * cw.min(ch)).clamp(0.15, 0.5);
                            let (x, y, sw, sh) = (x + w / 2.0, y + w / 2.0, cw - w, ch - w);
                            c.outlined((x, y, sw, sh), arms::FILTERED_FILL, (&Tone::hex(rgb), w));
                        }
                    }
                }
            }
        }
    }

    /// Drawn in the box `(x, y, w, h)` of `canvas` as `how` says: the
    /// squares (the arm behind, then the one in front), the box where both
    /// thresholds pass dashed, the axes and their titles, `title` over it
    /// (none when empty); and `marginals`, each `size` of the plot's
    /// shorter side deep (see [`framed`]): every square's column above,
    /// its row to the right, the signal arm's kept under its removed and
    /// the control arm's behind (or in front), each scaled to the busiest.
    pub fn draw_with(
        &self,
        canvas: &mut Canvas,
        (x, y, w, h): (f64, f64, f64, f64),
        title: &str,
        (marginals, size): (Marginals, f64),
        how: &PairLook,
    ) {
        let head = if title.is_empty() { 0.0 } else { 14.0 };
        let lay = framed(
            (x, y + head, w, h - head),
            self.xr,
            self.yr,
            marginals,
            size,
        );
        let plot = &lay.plot;
        if !title.is_empty() {
            canvas.text(plot.px, y + 10.0, title, 10.0, Anchor::Start, INK);
        }
        self.draw_squares(canvas, plot, how.arms.front);
        self.marginals(canvas, &lay, how.arms.front);
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
        let p_axis = how.arms.p_axis;
        let ticks = |k: usize, a: Axis| -> Vec<(f64, String)> {
            match a {
                Axis::Pvalue => (self.ticks[k].iter().zip(&self.tick_values[k]))
                    .map(|((at, _), &v)| (*at, p_axis.label(v)))
                    .collect(),
                _ => self.ticks[k].clone(),
            }
        };
        // Each axis' own title: the p-value as it reads, each arm's own
        // reads said so, as the histograms say them.
        let own = |a: Axis, per_arm: bool| match a {
            Axis::Pvalue => p_axis.name(),
            Axis::Coverage if per_arm => "reads per arm",
            Axis::Converted if per_arm => "converted reads per arm",
            _ => a.name(),
        };
        let own = [own(self.x, self.per_arm[0]), own(self.y, self.per_arm[1])];
        let titles = match how.titles {
            Some(t) => [t.x_or(self.x, own[0]), t.y_or(self.y, own[1])],
            None => own.map(String::from),
        };
        plot.axes(
            canvas,
            [&ticks(0, self.x), &ticks(1, self.y)],
            [&titles[0], &titles[1]],
        );
    }
}

impl Grid {
    /// The marginals `lay` has room for: each column's sites above, each
    /// row's to the right, the signal arm's kept under its removed (white,
    /// outlined), the control arm's light gray behind (or in front).
    fn marginals(&self, c: &mut Canvas, lay: &Framed, front: Front) {
        let roles = c.style().palette.roles();
        let n = self.squares.bins;
        let (cols, rows) = self.margins();
        let (mcols, mrows) = self.mutant_margins();
        // Each marginal scaled to its busiest bar of either arm.
        let busiest = |wt: &[Margin], m: &[f64]| {
            (wt.iter().map(|(k, o)| k + o))
                .chain(m.iter().copied())
                .fold(1.0, f64::max)
        };
        let (most_cols, most_rows) = (busiest(&cols, &mcols), busiest(&rows, &mrows));
        let plot = &lay.plot;
        let outline = 0.3;
        // A bar's foot and depth, from 0 to 1 of the marginal's, placed
        // by `place`; each arm's, the one behind first.
        type Place<'p> = &'p dyn Fn(f64, f64) -> (f64, f64, f64, f64);
        let bars = |c: &mut Canvas, (wt, m): (Margin, f64), most: f64, place: Place| {
            let (kept, other) = wt;
            let draw_wt = |c: &mut Canvas| {
                let (x, y, w, h) = place(0.0, kept / most);
                c.rect(x, y, w, h, roles.kept);
                // Too thin to show white inside its outline: none.
                let (_, _, ow, oh) = place(kept / most, other / most);
                if ow.min(oh) > 2.0 * outline {
                    let at = place(kept / most, other / most);
                    c.outlined(at, arms::FILTERED_FILL, (arms::OUTLINE, outline));
                }
            };
            let draw_mut = |c: &mut Canvas| {
                let (x, y, w, h) = place(0.0, m / most);
                c.rect(x, y, w, h, arms::MUT_FILL);
            };
            match front {
                Front::Wt => {
                    draw_mut(c);
                    draw_wt(c);
                }
                Front::Mut => {
                    draw_wt(c);
                    draw_mut(c);
                }
            }
        };
        if let Some((_, ty, _, th)) = lay.top {
            let cw = plot.pw / n as f64;
            for (k, (wt, m)) in cols.into_iter().zip(mcols).enumerate() {
                let bx = plot.px + k as f64 * cw;
                // From the marginal's foot up: `from` to `from + depth`.
                let place =
                    |from: f64, depth: f64| (bx, ty + th - (from + depth) * th, cw, depth * th);
                bars(c, (wt, m), most_cols, &place);
            }
        }
        if let Some((rx, _, rw, _)) = lay.right {
            let ch = plot.ph / n as f64;
            for (j, (wt, m)) in rows.into_iter().zip(mrows).enumerate() {
                let by = plot.py + plot.ph - (j + 1) as f64 * ch;
                let place = |from: f64, depth: f64| (rx + from * rw, by, depth * rw, ch);
                bars(c, (wt, m), most_rows, &place);
            }
        }
    }

    /// The signal arm's kept and removed per column (along x, what the top
    /// marginal draws) and per row (along y, the right one's); square `i`
    /// is in row `i / n`, column `i % n`.
    pub fn margins(&self) -> (Vec<Margin>, Vec<Margin>) {
        let n = self.squares.bins;
        let (mut cols, mut rows) = (vec![(0.0, 0.0); n], vec![(0.0, 0.0); n]);
        let squares = self.squares.kept.iter().zip(&self.squares.other);
        for (i, (&k, &o)) in squares.enumerate() {
            for at in [&mut cols[i % n], &mut rows[i / n]] {
                *at = (at.0 + k as f64, at.1 + o as f64);
            }
        }
        (cols, rows)
    }

    /// The control arm's sites per column and per row, as
    /// [`Grid::margins`] lays them out; none at all without it.
    pub fn mutant_margins(&self) -> (Vec<f64>, Vec<f64>) {
        let n = self.squares.bins;
        let (mut cols, mut rows) = (vec![0.0; n], vec![0.0; n]);
        for (i, &m) in self.mutant.iter().flatten().enumerate() {
            cols[i % n] += m as f64;
            rows[i / n] += m as f64;
        }
        (cols, rows)
    }
}

/// A marginal bar's sites, `(kept, removed)`.
pub type Margin = (f64, f64);

/// A 2D panel's plot and its marginals' boxes `(x, y, w, h)`.
pub struct Framed {
    pub plot: Plot,
    pub top: Option<(f64, f64, f64, f64)>,
    pub right: Option<(f64, f64, f64, f64)>,
}

/// Room between the plot and a marginal, in points of the 8-point layout.
const MARGINAL_GAP: f64 = 2.0;

/// The plot over `xr` by `yr` in `area` as [`Plot::on`] lays it out, with
/// room taken from it for `marginals`, each `size` of its shorter side
/// deep: the top one over its whole width, the right one beside its whole
/// height.
pub fn framed(
    area: (f64, f64, f64, f64),
    xr: (f64, f64),
    yr: (f64, f64),
    marginals: Marginals,
    size: f64,
) -> Framed {
    let mut plot = Plot::on(area, xr, yr);
    let side = size * plot.pw.min(plot.ph);
    let room = |on: bool| if on { side + MARGINAL_GAP } else { 0.0 };
    plot.pw = (plot.pw - room(marginals.right)).max(1.0);
    plot.ph = (plot.ph - room(marginals.top)).max(1.0);
    plot.py += room(marginals.top);
    let top = (marginals.top).then_some((plot.px, plot.py - MARGINAL_GAP - side, plot.pw, side));
    let right =
        (marginals.right).then_some((plot.px + plot.pw + MARGINAL_GAP, plot.py, side, plot.ph));
    Framed { plot, top, right }
}

/// `y` against `x` over `s` at `f` as the view draws its 2D panel
/// ([`Grid::draw_with`]), `bins` a side, the whole of a figure `w` by `h`
/// in `style`, with `marginals` each `size` deep, the arms and the axis
/// titles as `how` says; titled as the view's 2D figure titles it, but
/// alone (`--panel-only`).
pub fn view_pair(
    s: &Sites,
    (x, y): (Axis, Axis),
    f: &SiteFilterArgs,
    bins: usize,
    ((w, h), style): ((f64, f64), &crate::figure::style::Style),
    (marginals, how): ((Marginals, f64), &PairLook),
) -> String {
    let g = grid(s, (x, y), f, bins);
    let mut canvas = Canvas::styled(w, h, style);
    let (pw, ph) = canvas.page();
    let title = match style.panel_only {
        true => String::new(),
        false => format!("{} by {}", y.short(), x.short()),
    };
    g.draw_with(&mut canvas, (0.0, 0.0, pw, ph), &title, marginals, how);
    canvas.finish()
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
/// boxes, axes and axis titles only. The arms drawn as `arms` says.
pub fn figure(
    s: &Sites,
    pairs: &[(Axis, Axis, String)],
    f: &SiteFilterArgs,
    (w, h): (f64, f64),
    (title, arms): (&str, ArmLook),
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
        c.rect(12.0, 24.0, 9.0, 9.0, roles.kept);
        c.text(
            25.0,
            32.0,
            "WT kept by every threshold",
            9.0,
            Anchor::Start,
            MUTED,
        );
        c.outlined(
            (150.5, 24.5, 8.0, 8.0),
            arms::FILTERED_FILL,
            (arms::OUTLINE, 0.6),
        );
        c.text(
            163.0,
            32.0,
            "WT removed by at least one",
            9.0,
            Anchor::Start,
            MUTED,
        );
        let mut x = 300.0;
        if s.table.has_control() {
            c.rect(x, 24.0, 9.0, 9.0, arms::MUT_FILL);
            c.text(
                x + 13.0,
                32.0,
                "MUT (reads and edit ratio)",
                9.0,
                Anchor::Start,
                MUTED,
            );
            x += 140.0;
        }
        let note = "darker: more sites · dashed: both drawn thresholds pass";
        c.text(x, 32.0, note, 9.0, Anchor::Start, MUTED);
    }
    // Room over the panels for the key, and under each for its notes.
    let (top, notes) = if alone { (0.0, 0.0) } else { (40.0, 22.0) };
    let cols = 2usize;
    let rows = pairs.len().div_ceil(cols).max(1);
    let (pw, ph) = ((w - 12.0) / cols as f64, (h - top) / rows as f64);
    for (k, ((x, y, why), g)) in pairs.iter().zip(&grids).enumerate() {
        let (cx, cy) = (6.0 + (k % cols) as f64 * pw, top + (k / cols) as f64 * ph);
        let how = PairLook { arms, titles: None };
        let none = (Marginals::default(), 0.0);
        if alone {
            g.draw_with(c, (cx, cy, pw - 6.0, ph), "", none, &how);
            continue;
        }
        let title = format!("{} by {}", y.short(), x.short());
        g.draw_with(c, (cx, cy, pw - 6.0, ph - notes), &title, none, &how);
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
