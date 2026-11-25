//! Each batch's cell QC as a figure: its barcodes by rank against their
//! counts (the knee), and genes detected against counts, each barcode
//! coloured by the verdict data-beans' cell QC gave it (kept, a MAD outlier
//! or near empty).

use super::density::{decades, Plot, Squares, BINS};
use crate::figure::style::{Palette, Style};
use crate::figure::{Anchor, Canvas, INK, MUTED};

/// A cells figure's own size, in points of its 8-point layout, when the
/// page sets neither side.
pub const SIZE: (f64, f64) = (720.0, 320.0);
use data_beans::qc_lib::QcReport;

/// The verdict on a barcode.
#[derive(Clone, Copy, Debug, PartialEq, Eq, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "snake_case")]
pub enum Verdict {
    Kept,
    /// Dropped as an outlier on detected genes or counts.
    Outlier,
    /// Under the near-empty floor of detected genes.
    NearEmpty,
}

impl Verdict {
    pub const ALL: [Verdict; 3] = [Verdict::Kept, Verdict::Outlier, Verdict::NearEmpty];

    /// Its colour in `palette`: kept and removed as the sites', near empty
    /// a mid gray.
    pub fn colour(self, palette: Palette) -> &'static str {
        let roles = palette.roles();
        match self {
            Verdict::Kept => roles.kept,
            Verdict::Outlier => roles.removed,
            Verdict::NearEmpty => roles.near_empty,
        }
    }

    pub fn said(self) -> &'static str {
        match self {
            Verdict::Kept => "kept",
            Verdict::Outlier => "MAD outlier",
            Verdict::NearEmpty => "near empty",
        }
    }
}

/// Points along the knee's log rank axis, at most (and where it turns).
pub const KNEE_STEPS: usize = 400;

/// The cells' density: its ranges (log10), and kept and the rest per
/// square.
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct Grid {
    pub xr: (f64, f64),
    pub yr: (f64, f64),
    pub squares: Squares,
}

/// A batch's barcodes: each one's counts, detected genes and verdict.
#[derive(Clone, Debug, serde::Serialize, serde::Deserialize)]
pub struct Cells {
    pub batch: Box<str>,
    pub counts: Vec<f32>,
    pub genes: Vec<f32>,
    pub verdict: Vec<Verdict>,
}

impl Cells {
    /// From the cell QC's report on `batch`.
    pub fn of(batch: &str, r: &QcReport) -> Self {
        let verdict = (0..r.total_counts.len())
            .map(|c| match (r.near_empty.get(c), r.train_keep.get(c)) {
                (Some(true), _) => Verdict::NearEmpty,
                (_, Some(false)) => Verdict::Outlier,
                _ => Verdict::Kept,
            })
            .collect();
        Self {
            batch: batch.into(),
            counts: r.total_counts.clone(),
            genes: r.n_genes.clone(),
            verdict,
        }
    }

    /// Its figure's file stem, under `qc_plots/`.
    pub fn stem(&self) -> String {
        format!("{}_cells", self.batch)
    }

    /// Barcodes with each verdict, in [`Verdict::ALL`]'s order.
    pub fn tally(&self) -> [usize; 3] {
        let mut n = [0; 3];
        for v in &self.verdict {
            n[*v as usize] += 1;
        }
        n
    }

    /// The figure, `w` by `h` in `style`: the knee on the left, genes by
    /// counts on the right, the verdicts' key and counts over them; alone
    /// (`--panel-only`), the two plots only. Laid out at 8 points on the
    /// page the style's font scales.
    pub fn figure(&self, (w, h): (f64, f64), style: &Style) -> String {
        let mut c = Canvas::styled(w, h, style);
        let (w, h) = c.page();
        let half = (w - 12.0) / 2.0;
        let top = if style.panel_only { 0.0 } else { 40.0 };
        if !style.panel_only {
            self.key(&mut c, half, style.palette);
        }
        let palette = style.palette;
        self.knee(&mut c, (6.0, top, half - 6.0, h - top - 6.0), palette);
        self.density(&mut c, (6.0 + half, top, half - 6.0, h - top - 6.0));
        c.finish()
    }

    /// The title, the verdicts' key and what the right plot's colours say,
    /// over plots `half` wide.
    fn key(&self, c: &mut Canvas, half: f64, palette: Palette) {
        let [kept, outlier, empty] = self.tally();
        let title = format!(
            "{} cells: {kept} kept, {outlier} MAD outliers, {empty} near empty, of {}",
            self.batch,
            self.counts.len()
        );
        c.text(12.0, 16.0, &title, 11.0, Anchor::Start, INK);
        let mut x = 12.0;
        for v in Verdict::ALL {
            c.rect(x, 24.0, 9.0, 9.0, v.colour(palette));
            c.text(x + 13.0, 32.0, v.said(), 9.0, Anchor::Start, MUTED);
            x += 100.0;
        }
        let r = palette.roles();
        let note = format!(
            "right: a square {} where most of its barcodes are kept, {} where not",
            r.kept_name, r.removed_name
        );
        c.text(12.0 + half, 32.0, &note, 9.0, Anchor::Start, MUTED);
    }

    /// The knee as drawn: each barcode's log10 rank and log10 counts,
    /// thinned to a point each 1/`steps` of the log rank axis and wherever
    /// the verdict turns, the last always; and that axis' end.
    pub fn knee_points(&self, steps: usize) -> (Vec<(f64, f64, Verdict)>, f64) {
        let mut order: Vec<usize> = (0..self.counts.len()).collect();
        order.sort_by(|&a, &b| self.counts[b].total_cmp(&self.counts[a]));
        let end = (order.len().max(2) as f64).log10();
        let gap = end / steps.max(1) as f64;
        let mut out: Vec<(f64, f64, Verdict)> = Vec::new();
        for (r, &i) in order.iter().enumerate() {
            let x = (r as f64 + 1.0).log10();
            let p = (x, (self.counts[i] as f64).max(1.0).log10(), self.verdict[i]);
            let keep = match out.last() {
                None => true,
                Some(q) => x - q.0 >= gap || q.2 != p.2 || r + 1 == order.len(),
            };
            if keep {
                out.push(p);
            }
        }
        (out, end)
    }

    /// Counts by rank, both on log scales, the curve coloured by verdict.
    fn knee(&self, c: &mut Canvas, area: (f64, f64, f64, f64), palette: Palette) {
        let (points, end) = self.knee_points(KNEE_STEPS);
        let top = points.first().map_or(1.0, |p| p.1);
        let plot = Plot::on(area, (0.0, end), (0.0, top.max(1.0)));
        for pair in points.windows(2) {
            let (a, b) = (plot.at(pair[0].0, pair[0].1), plot.at(pair[1].0, pair[1].1));
            c.line(a.0, a.1, b.0, b.1, pair[1].2.colour(palette), 1.6);
        }
        let ticks = [decades(plot.xr), decades(plot.yr)];
        plot.axes(c, [&ticks[0], &ticks[1]], ["barcode rank", "counts"]);
    }

    /// Genes detected by counts (log10, from the decade below the fewest
    /// to the most), binned [`BINS`] a side: the ranges, then kept and the
    /// rest per square (row-major from the bottom left).
    pub fn density_grid(&self) -> Grid {
        let lg =
            |v: &[f32]| -> Vec<f64> { v.iter().map(|&x| (x as f64).max(1.0).log10()).collect() };
        let (lx, ly) = (lg(&self.counts), lg(&self.genes));
        // From the decade below the fewest to the most: no empty decades.
        let span = |v: &[f64]| {
            let lo = v.iter().copied().fold(f64::INFINITY, f64::min);
            let hi = v.iter().copied().fold(1.0, f64::max);
            (lo.floor().clamp(0.0, (hi - 1.0).max(0.0)), hi)
        };
        let (xr, yr) = (span(&lx), span(&ly));
        let (mut kept, mut other) = (vec![0u32; BINS * BINS], vec![0u32; BINS * BINS]);
        let slot = |v: f64, (lo, hi): (f64, f64)| {
            (((v - lo) / (hi - lo) * BINS as f64) as usize).min(BINS - 1)
        };
        for i in 0..lx.len() {
            let at = slot(ly[i], yr) * BINS + slot(lx[i], xr);
            if self.verdict[i] == Verdict::Kept {
                kept[at] += 1;
            } else {
                other[at] += 1;
            }
        }
        Grid {
            xr,
            yr,
            squares: Squares::new(BINS, kept, other),
        }
    }

    /// Genes detected by counts (log scales), a square kept or not by
    /// which most of its barcodes are, deeper as it holds more.
    fn density(&self, c: &mut Canvas, area: (f64, f64, f64, f64)) {
        let g = self.density_grid();
        let plot = Plot::on(area, g.xr, g.yr);
        g.squares.draw(c, &plot);
        let ticks = [decades(g.xr), decades(g.yr)];
        plot.axes(c, [&ticks[0], &ticks[1]], ["counts", "genes detected"]);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn report() -> QcReport {
        QcReport {
            train_keep: vec![true, true, false, true],
            near_empty: vec![false, false, false, true],
            feature_keep: vec![],
            n_genes: vec![900.0, 700.0, 9000.0, 1.0],
            total_counts: vec![3000.0, 2000.0, 60000.0, 2.0],
            mito_frac: None,
            ribo_frac: None,
            n_cells_dropped: 1,
            n_features_dropped: 0,
        }
    }

    /// The verdicts read off the report: near empty first, then an
    /// outlier, else kept; the figure names them and their counts.
    #[test]
    fn the_verdicts_and_figure_follow_the_report() {
        let c = Cells::of("wt_rep1", &report());
        assert_eq!(
            c.verdict,
            [
                Verdict::Kept,
                Verdict::Kept,
                Verdict::Outlier,
                Verdict::NearEmpty
            ]
        );
        assert_eq!(c.tally(), [2, 1, 1]);
        let svg = c.figure((720.0, 320.0), &Default::default());
        assert!(svg.contains("wt_rep1 cells: 2 kept, 1 MAD outliers, 1 near empty, of 4"));
        assert!(svg.contains("barcode rank") && svg.contains("genes detected"));
        assert!(crate::figure::png(&svg, 1.0).is_ok());
    }
}
