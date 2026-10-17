//! The browser's rows: which there are, what each is titled and draws,
//! and the genes' lanes.

use super::*;

/// A row of the browser.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Row {
    Track(usize),
    Contrast(Contrast),
    /// One track above a zero line and another below it.
    Mirror(usize, usize),
    /// The annotation's genes in view.
    Genes,
}

impl PileupView<'_> {
    /// The rows top to bottom: contrast, tracks, then genes; a track
    /// flagged to mirror takes the next into its row, unless split.
    pub(super) fn rows(&self) -> Vec<Row> {
        let mut rows: Vec<Row> = self.contrast.map(Row::Contrast).into_iter().collect();
        let mut k = 0;
        while k < self.tracks.len() {
            if self.mirror && self.tracks[k].mirror_with_next {
                rows.push(Row::Mirror(k, k + 1));
                k += 2;
            } else {
                rows.push(Row::Track(k));
                k += 1;
            }
        }
        rows.extend((!self.genes.is_empty()).then_some(Row::Genes));
        rows
    }

    /// The top of `row`'s scale: the shared one, but for a track on its own
    /// scale, locked or not.
    pub(super) fn row_max(&self, row: Row, bins: &Bins) -> Option<f64> {
        match row {
            Row::Mirror(k, _) | Row::Track(k) if self.tracks[k].own_scale => None,
            _ => bins.shared,
        }
    }

    /// Whether some track shares a row with the next (`v` toggles it).
    pub(super) fn mirrorable(&self) -> bool {
        self.tracks.iter().any(|t| t.mirror_with_next)
    }

    pub(super) fn row_title(&self, row: Row) -> String {
        match row {
            Row::Track(i) => {
                let t = &self.tracks[i];
                format!("{} · {}", t.label, self.y_title(i))
            }
            Row::Contrast(c) => {
                let (a, b) = (self.tracks[c.a].label, self.tracks[c.b].label);
                format!("{a} vs {b} · {}", c.measure.name(self.on))
            }
            Row::Mirror(a, b) => {
                let y = self.y_title(a);
                let (a, b) = (&self.tracks[a], &self.tracks[b]);
                format!("{} above, {} below · {y}", a.label, b.label)
            }
            Row::Genes => "genes".into(),
        }
    }

    /// Genes overlapping the window, each with a lane so labels do not collide.
    pub(super) fn gene_lanes(&self) -> Vec<(usize, &GeneModel)> {
        let (lo, hi) = self.window;
        let gap = ((hi - lo) as f64 * LANE_GAP) as i64;
        let mut ends: Vec<i64> = Vec::new();
        let mut out = Vec::new();
        for g in self.genes.iter().filter(|g| g.hi > lo && g.lo <= hi) {
            let lane = match ends.iter().position(|&end| end + gap < g.lo) {
                Some(l) => l,
                None if ends.len() < MAX_LANES => {
                    ends.push(i64::MIN);
                    ends.len() - 1
                }
                None => continue,
            };
            ends[lane] = g.hi;
            out.push((lane, g));
        }
        out
    }

    /// Gene lanes in view (at least one).
    pub(super) fn lanes(&self) -> usize {
        let lanes = self.gene_lanes().iter().map(|(l, _)| l + 1).max();
        lanes.unwrap_or(1)
    }

    /// Every track binned over the window, and as drawn.
    pub(super) fn bins(&self) -> Bins {
        let edges = self.edges();
        let tracks: Vec<Binned> = self.tracks.iter().map(|t| t.bin(&edges)).collect();
        let shown: Vec<Shown> = (0..self.tracks.len())
            .map(|k| self.shown(k, &edges, &tracks[k]))
            .collect();
        let tallest = self
            .tracks
            .iter()
            .zip(&shown)
            .filter(|(t, _)| t.is_reads() && !t.own_scale)
            .map(|(_, s)| s.values.iter().fold(0.0, |m: f64, &v| m.max(v)))
            .fold(0.0, f64::max);
        Bins {
            edges,
            tracks,
            shown,
            shared: self.y_lock.or((tallest > 0.0).then_some(tallest)),
        }
    }

    /// Track `k` as [`Show`] asks, from its binned reads. A track with no
    /// total (or a depth track) draws as it is, but for its sites; what a
    /// mode needs beyond the binned reads is binned only in that mode.
    pub(super) fn shown(&self, k: usize, edges: &BinEdges, (front, behind): &Binned) -> Shown {
        let t = &self.tracks[k];
        let shown = |values: Vec<f64>, front: Option<Vec<f64>>, accent: bool| Shown {
            values,
            front,
            accent: accent || t.accent,
        };
        match (self.show, behind) {
            (Show::Sites, _) if t.is_reads() => {
                let ones: Vec<(i64, f64)> = distinct_positions(t.front())
                    .into_iter()
                    .map(|p| (p, 1.0))
                    .collect();
                shown(edges.bin(&ones, false), None, false)
            }
            (Show::Converted, Some(_)) => shown(front.clone(), None, true),
            (Show::Unconverted, Some(_)) => {
                let total = t.behind().unwrap_or_default();
                shown(
                    edges.bin(&unconverted(t.front(), total), t.log),
                    None,
                    false,
                )
            }
            (Show::Both, Some(total)) => shown(total.clone(), Some(front.clone()), true),
            _ => shown(front.clone(), None, false),
        }
    }

    /// What track `k`'s bars count, for its title and y axis.
    pub(super) fn y_title(&self, k: usize) -> String {
        let t = &self.tracks[k];
        let sites = self.show == Show::Sites && t.is_reads();
        if sites || t.behind().is_some() {
            self.show.label(self.on, self.off)
        } else {
            t.name.to_string()
        }
    }

    /// The contrast per bar, from raw sums (rebinned when the tracks show logs).
    pub(super) fn contrast_values(&self, c: Contrast, bins: &Bins) -> Vec<Option<f64>> {
        let counts = |k: usize| -> (Vec<f64>, Vec<f64>) {
            let t = &self.tracks[k];
            match &bins.tracks[k] {
                (front, Some(total)) if !t.log => (front.clone(), total.clone()),
                _ => (
                    bins.edges.bin(t.front(), false),
                    bins.edges.bin(t.behind().unwrap_or_default(), false),
                ),
            }
        };
        let ((ma, na), (mb, nb)) = (counts(c.a), counts(c.b));
        (0..ma.len())
            .map(|k| c.measure.of((ma[k], na[k]), (mb[k], nb[k])))
            .collect()
    }

    pub(super) fn plot_count(rows: &[Row]) -> usize {
        rows.iter().filter(|&&r| r != Row::Genes).count().max(1)
    }
}
