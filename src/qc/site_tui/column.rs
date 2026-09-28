//! Each criterion's display rules, and the focused column binned.

use super::*;

/// How each criterion is drawn.
impl Criterion {
    pub(super) fn label(self) -> &'static str {
        match self {
            Criterion::MaxPv => "p-value",
            Criterion::MinLogOdds => "log odds",
            Criterion::MinFold => "fold",
            Criterion::MinCoverage => "coverage",
            Criterion::MinConverted => "converted",
            Criterion::MinEditRatio => "min edit ratio",
            Criterion::MaxEditRatio => "max edit ratio",
            Criterion::MinCells => "cells",
        }
    }

    /// What the histogram's x axis shows.
    pub(super) fn axis(self) -> &'static str {
        match self {
            Criterion::MaxPv => "-log10 p",
            Criterion::MinLogOdds => "raw log odds ratio",
            Criterion::MinFold => "signal / control edit rate",
            Criterion::MinCoverage => "signal + control reads",
            Criterion::MinConverted => "converted signal reads",
            Criterion::MinEditRatio => "converted / coverage: low is weak editing",
            Criterion::MaxEditRatio => "converted / coverage: high is variant-like",
            Criterion::MinCells => "kept cells with a converted read",
        }
    }

    /// Whether the site keeps the high end of the displayed axis (all but
    /// the upper edit-ratio bound: a p-value is shown as -log10 p).
    pub(super) fn keeps_high(self) -> bool {
        self != Criterion::MaxEditRatio
    }

    /// Whether bin `k` lies on the dropped side of the threshold's bin.
    pub(super) fn drops(self, k: i32, threshold: Option<i32>) -> bool {
        match threshold {
            Some(p) if self.keeps_high() => k < p,
            Some(p) => k > p,
            None => false,
        }
    }

    /// The threshold as shown: `≥ v` or `≤ v`, `None` when off.
    pub(super) fn shown(self, f: &SiteFilterArgs) -> Option<String> {
        let op = if self.is_max() { "≤" } else { "≥" };
        (!self.is_off(f)).then(|| format!("{op} {}", self.fmt_short(self.get(f))))
    }

    pub(super) fn is_integer(self) -> bool {
        matches!(
            self,
            Criterion::MinCoverage | Criterion::MinConverted | Criterion::MinCells
        )
    }

    pub(super) fn scales(self) -> &'static [Scale] {
        match self {
            Criterion::MinLogOdds => &[Scale::Linear],
            // Log bins resolve the p ~ 0.05 region the decision is made in.
            Criterion::MaxPv => &[Scale::Log, Scale::Sqrt, Scale::Linear],
            Criterion::MinEditRatio | Criterion::MaxEditRatio => &[Scale::Linear, Scale::Sqrt],
            _ => &[Scale::Log, Scale::Sqrt, Scale::Linear],
        }
    }

    /// Position of a raw value on the displayed axis (may be infinite).
    pub(super) fn display(self, raw: f64) -> f64 {
        match self {
            Criterion::MaxPv if raw <= 0.0 => f64::INFINITY,
            Criterion::MaxPv => (-raw.log10()).max(0.0),
            _ => raw,
        }
    }

    /// Exact, for the flags.
    pub(super) fn fmt(self, v: f64) -> String {
        if self.is_integer() {
            format!("{}", v as u64)
        } else {
            format!("{}", v as f32)
        }
    }

    /// Three significant digits, for the screen.
    pub(super) fn fmt_short(self, v: f64) -> String {
        let a = v.abs();
        if self.is_integer() || !v.is_finite() || a == 0.0 {
            self.fmt(v)
        } else if !(1e-3..1e5).contains(&a) {
            format!("{v:.2e}")
        } else {
            let decimals = (2 - a.log10().floor() as i32).max(0) as usize;
            let s = format!("{v:.decimals$}");
            if s.contains('.') {
                s.trim_end_matches('0').trim_end_matches('.').to_string()
            } else {
                s
            }
        }
    }

    pub(super) fn bit(self) -> u8 {
        1 << self as u8
    }
}

/// The focused criterion's column, binned.
pub(super) struct Column {
    /// Display value per site, infinities clamped to the finite range.
    pub(super) display: Vec<f32>,
    pub(super) n_inf: usize,
    pub(super) lo: f64,
    pub(super) hi: f64,
    pub(super) hist: Binned,
    /// Histogram slot per site.
    pub(super) slot: Vec<u16>,
    /// Thresholds the bin steps visit, ascending on the displayed axis:
    /// `(bin, raw)` with `raw` the loosest threshold that keeps the whole bin
    /// and drops the bins on the far side. Stops that keep everything are
    /// left out: they read as off.
    pub(super) stops: Vec<(i32, f64)>,
    /// The smallest display value in each bin, for tick labels.
    pub(super) lowest: Vec<f64>,
}

impl Column {
    pub(super) fn new(view: &SiteView, c: Criterion, scale: Scale) -> Self {
        let mut display: Vec<f32> = (0..view.table.len())
            .map(|i| c.display(c.raw(view.table, i, view.cells(i))) as f32)
            .collect();
        let (mut lo, mut hi, mut n_inf) = (f64::INFINITY, f64::NEG_INFINITY, 0);
        for &d in &display {
            if d.is_finite() {
                lo = lo.min(d as f64);
                hi = hi.max(d as f64);
            } else {
                n_inf += 1;
            }
        }
        if lo > hi {
            (lo, hi) = (0.0, 0.0);
        }
        display
            .iter_mut()
            .for_each(|d| *d = d.clamp(lo as f32, hi as f32));
        let (hist, slot, stops, lowest) = bin(view, c, scale, &display, lo, hi);
        Self {
            display,
            n_inf,
            lo,
            hi,
            hist,
            slot,
            stops,
            lowest,
        }
    }

    pub(super) fn rebin(&mut self, view: &SiteView, c: Criterion, scale: Scale) {
        (self.hist, self.slot, self.stops, self.lowest) =
            bin(view, c, scale, &self.display, self.lo, self.hi);
    }

    /// Bin key of a raw threshold, clamped to the histogram.
    pub(super) fn key_of(&self, c: Criterion, raw: f64) -> i32 {
        let d = c.display(raw).clamp(self.lo, self.hi);
        let kmax = self.hist.kmax();
        self.hist.bins.key(d).clamp(self.hist.kmin, kmax)
    }
}

/// A column binned: the histogram, each site's slot, the stops, and each
/// bin's smallest value.
pub(super) type ColumnBins = (Binned, Vec<u16>, Vec<(i32, f64)>, Vec<f64>);

/// Bin `display` on `scale`: the histogram, each site's slot, and the stops.
pub(super) fn bin(
    view: &SiteView,
    c: Criterion,
    scale: Scale,
    display: &[f32],
    lo: f64,
    hi: f64,
) -> ColumnBins {
    // A signed axis spans -max..max; widen so it still gets ~50 bins.
    let max = if lo < 0.0 {
        2.0 * lo.abs().max(hi.abs())
    } else {
        hi
    };
    let bins = Binning::new(scale, max, c.is_integer());
    let kmin = bins.key(lo);
    let n = (bins.key(hi) - kmin + 1).clamp(1, MAX_BINS) as usize;
    let mut counts = vec![0; n];
    let mut slot = Vec::with_capacity(display.len());
    let mut edge: Vec<Option<f64>> = vec![None; n];
    let mut lowest = vec![f64::INFINITY; n];
    for (i, &d) in display.iter().enumerate() {
        let b = (bins.key(d as f64) - kmin).clamp(0, n as i32 - 1) as usize;
        counts[b] += 1;
        lowest[b] = lowest[b].min(d as f64);
        slot.push(b as u16);
        let raw = c.raw(view.table, i, view.cells(i));
        edge[b] = Some(match edge[b] {
            None => raw,
            Some(e) if c.is_max() => e.max(raw),
            Some(e) => e.min(raw),
        });
    }
    let stops = edge
        .iter()
        .enumerate()
        .filter_map(|(b, e)| e.map(|raw| (kmin + b as i32, raw)))
        .filter(|&(_, raw)| {
            let mut f = SiteFilterArgs::permissive();
            c.set(&mut f, raw);
            !c.is_off(&f)
        })
        .collect();
    (Binned { bins, kmin, counts }, slot, stops, lowest)
}
