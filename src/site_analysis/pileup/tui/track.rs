//! A track's values and how they bin onto the view's columns.

use super::*;

/// What a track holds.
pub enum Values<'a> {
    /// Sorted `(position, value)` pairs, optionally in front of a per-position total.
    Reads {
        front: &'a [(i64, f64)],
        behind: Option<&'a [(i64, f64)]>,
    },
    /// Genomic bins `(start, end, value)`; each column shows the bin covering it.
    Ranges(&'a [(i64, i64, f64)]),
}

/// One track of the browser.
pub struct Track<'a> {
    pub label: &'a str,
    /// What the values are, for the panel title.
    pub name: &'a str,
    pub values: Values<'a>,
    /// Bins become `log10(1 + sum)`, as the printed pileup does.
    pub log: bool,
}

/// A track binned over the window: front, and the total when there is one.
pub(super) type Binned = (Vec<f64>, Option<Vec<f64>>);

impl<'a> Track<'a> {
    /// A read track with no total.
    pub fn single(label: &'a str, name: &'a str, front: &'a [(i64, f64)], log: bool) -> Self {
        Track {
            label,
            name,
            values: Values::Reads {
                front,
                behind: None,
            },
            log,
        }
    }

    /// This read track drawn in front of `total`.
    pub fn with_total(mut self, total: &'a [(i64, f64)]) -> Self {
        if let Values::Reads { behind, .. } = &mut self.values {
            *behind = Some(total);
        }
        self
    }

    /// A read-depth track over genomic bins.
    pub fn depth(label: &'a str, ranges: &'a [(i64, i64, f64)]) -> Self {
        Track {
            label,
            name: "reads per depth bin",
            values: Values::Ranges(ranges),
            log: false,
        }
    }

    /// The read positions; none for a depth track.
    pub(super) fn front(&self) -> &'a [(i64, f64)] {
        match self.values {
            Values::Reads { front, .. } => front,
            Values::Ranges(_) => &[],
        }
    }

    /// The total behind the reads, when there is one.
    pub(super) fn behind(&self) -> Option<&'a [(i64, f64)]> {
        match self.values {
            Values::Reads { behind, .. } => behind,
            Values::Ranges(_) => None,
        }
    }

    pub(super) fn is_reads(&self) -> bool {
        matches!(self.values, Values::Reads { .. })
    }

    /// Binned over `edges`.
    pub(super) fn bin(&self, edges: &BinEdges) -> Binned {
        match self.values {
            Values::Reads { front, behind } => (
                edges.bin(front, self.log),
                behind.map(|b| edges.bin(b, self.log)),
            ),
            Values::Ranges(ranges) => (ranges_per_column(ranges, edges), None),
        }
    }

    /// Distinct positions inside `lo..=hi`.
    pub(super) fn sites_in(&self, lo: i64, hi: i64) -> Vec<i64> {
        let front = self.front();
        let a = front.partition_point(|p| p.0 < lo);
        let b = front.partition_point(|p| p.0 <= hi);
        distinct_positions(&front[a..b])
    }
}

/// Every track binned over one window, once per frame.
pub(super) struct Bins {
    pub(super) edges: BinEdges,
    pub(super) tracks: Vec<Binned>,
    /// Shared y-axis top over read tracks (not depth): the tallest bar, total included.
    pub(super) shared: Option<f64>,
}

/// Per column, the value of the genomic bin covering the column's middle.
pub(super) fn ranges_per_column(ranges: &[(i64, i64, f64)], edges: &BinEdges) -> Vec<f64> {
    (0..edges.num_bins)
        .map(|k| {
            let (start, stop) = edges.col_range(k);
            let mid = (start + stop) / 2;
            let i = ranges.partition_point(|r| r.1 <= mid);
            ranges.get(i).filter(|r| r.0 <= mid).map_or(0.0, |r| r.2)
        })
        .collect()
}

/// How the contrast row compares two tracks' converted fractions per bar.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Measure {
    /// Fraction A minus fraction B, in percentage points.
    Difference,
    /// log2 of fraction A over fraction B, half a read of pseudocount.
    Log2Fold,
}

impl Measure {
    pub(super) fn name(self, on: &str) -> String {
        match self {
            Measure::Difference => format!("{on} fraction difference (pp)"),
            Measure::Log2Fold => format!("log2 fold of {on} fraction"),
        }
    }

    pub(super) fn next(self) -> Self {
        match self {
            Measure::Difference => Measure::Log2Fold,
            Measure::Log2Fold => Measure::Difference,
        }
    }

    /// The measure for converted `(ma, na)` of A and `(mb, nb)` of B; `None` without reads.
    pub fn of(self, (ma, na): (f64, f64), (mb, nb): (f64, f64)) -> Option<f64> {
        if na <= 0.0 || nb <= 0.0 {
            return None;
        }
        Some(match self {
            Measure::Difference => 100.0 * (ma / na - mb / nb),
            Measure::Log2Fold => {
                ((ma + 0.5) / (na + 1.0)).log2() - ((mb + 0.5) / (nb + 1.0)).log2()
            }
        })
    }

    pub(super) fn label(self, v: f64) -> String {
        match self {
            Measure::Difference => format!("{v:+.1}"),
            Measure::Log2Fold => format!("{v:+.2}"),
        }
    }
}

/// Two tracks compared bar by bar, `a` against `b`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) struct Contrast {
    pub(super) a: usize,
    pub(super) b: usize,
    pub(super) measure: Measure,
}
