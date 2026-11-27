//! Shared genomic binning for the Miami figure and the ASCII pileup.
//!
//! [`BinEdges`] fixes a `[min_pos, max_pos]` extent over `num_bins` equal
//! columns. Every track (matrix epi signal, BAM read depth, site markers)
//! maps positions through the SAME edges so the stacked tracks line up.
//! The mapping rule is the historical pileup rule (`pos_to_col` /
//! `bin_positions_with_extent` in `pileup.rs`), lifted here so both the
//! ASCII path and the figure share one definition.

use crate::figure::{scaled, unscaled};
use crate::tui::ui::Scale;

/// A fixed binning grid over an inclusive genomic extent: equal columns
/// of the bases from `min_pos`, or of their square root or log10(1 + bp)
/// (`x_scale`), so the left end is spread and the right end pressed.
#[derive(Clone, Copy, Debug)]
pub struct BinEdges {
    pub min_pos: i64,
    pub max_pos: i64,
    pub num_bins: usize,
    pub x_scale: Scale,
}

impl BinEdges {
    pub fn new(min_pos: i64, max_pos: i64, num_bins: usize) -> Self {
        Self {
            min_pos,
            max_pos,
            num_bins: num_bins.max(1),
            x_scale: Scale::Linear,
        }
    }

    /// The same edges on `x_scale`: the bases from the left end on it.
    pub fn on(self, x_scale: Scale) -> Self {
        Self { x_scale, ..self }
    }

    /// Whether its columns are equal in bases.
    pub fn is_linear(&self) -> bool {
        self.x_scale == Scale::Linear
    }

    /// How far across `pos` falls, 0 at `min_pos` and 1 at `max_pos`, on
    /// its scale; clamped to the extent.
    pub fn frac(&self, pos: i64) -> f64 {
        let rel = (pos - self.min_pos).clamp(0, self.span() as i64) as f64;
        scaled(self.x_scale, rel) / scaled(self.x_scale, self.span() as f64)
    }

    /// The first position at least `frac` of the way across (see
    /// [`BinEdges::frac`]): a column's first base, as [`BinEdges::col_of`]
    /// bins it.
    pub fn pos_at(&self, frac: f64) -> i64 {
        let top = scaled(self.x_scale, self.span() as f64);
        let rel = unscaled(self.x_scale, frac.clamp(0.0, 1.0) * top);
        // Within a rounding error of a whole base, that base.
        self.min_pos + (rel - 1e-9).ceil().max(0.0) as i64
    }

    /// Genomic span in bp (at least 1, so a degenerate single-position
    /// extent doesn't divide by zero).
    #[inline]
    pub fn span(&self) -> u64 {
        (self.max_pos - self.min_pos).max(1) as u64
    }

    /// Map a coordinate to its bin column. Positions at or below `min_pos`
    /// land in column 0; `max_pos` lands in the last column.
    #[inline]
    pub fn col_of(&self, pos: i64) -> usize {
        if !self.is_linear() {
            let col = (self.frac(pos) * self.num_bins as f64) as usize;
            return col.min(self.num_bins - 1);
        }
        let rel = (pos - self.min_pos).max(0) as u64;
        (rel * self.num_bins as u64 / self.span()).min(self.num_bins as u64 - 1) as usize
    }

    /// Genomic range `[start, stop)` of column `col`: on a sqrt or log x, a
    /// column at the spread left end narrower than a base holds none
    /// (`start == stop`).
    pub fn col_range(&self, col: usize) -> (i64, i64) {
        let (span, n) = (self.span() as i64, self.num_bins as i64);
        let at = |i: i64| match self.is_linear() {
            true => self.min_pos + i * span / n,
            false => self.pos_at(i as f64 / n as f64),
        };
        (at(col as i64), at(col as i64 + 1))
    }

    /// Bin a position/value list into `num_bins` summed columns. Positions
    /// outside `[min_pos, max_pos]` are dropped. When `log_transform`,
    /// each bin becomes `log10(1 + sum)`.
    pub fn bin(&self, positions: &[(i64, f64)], log_transform: bool) -> Vec<f64> {
        let mut bins = vec![0.0f64; self.num_bins];
        for &(pos, val) in positions {
            if pos < self.min_pos || pos > self.max_pos {
                continue;
            }
            bins[self.col_of(pos)] += val;
        }
        if log_transform {
            for v in bins.iter_mut() {
                *v = (1.0 + *v).log10();
            }
        }
        bins
    }

    /// Pixel x of a genomic position within `[x_left, x_left + plot_w]`,
    /// clamped to the plot band. Shared by every track + the gene model so
    /// they align horizontally.
    #[inline]
    pub fn x_px(&self, pos: i64, x_left: f32, plot_w: f32) -> f32 {
        x_left + self.frac(pos) as f32 * plot_w
    }
}

/// 99th-percentile of positive values — a spike-robust scale so one
/// outlier bin doesn't flatten every panel. Mirrors
/// `senna .../strand/place.rs::robust_max`: quickselect, O(n).
pub fn robust_max(values: impl Iterator<Item = f64>) -> f64 {
    let mut v: Vec<f64> = values.filter(|x| *x > 0.0).collect();
    if v.is_empty() {
        return 0.0;
    }
    let idx = (((v.len() - 1) as f64) * 0.99).round() as usize;
    *v.select_nth_unstable_by(idx, |a, b| a.total_cmp(b)).1
}

#[cfg(test)]
mod tests;
