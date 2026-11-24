//! Each batch's cell QC as a figure: its barcodes by rank against their
//! counts (the knee), and genes detected against counts, each barcode
//! coloured by the verdict data-beans' cell QC gave it (kept, a MAD outlier
//! or near empty).

use crate::figure::{Anchor, Canvas, FAINT, INK, MUTED};
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

    pub fn colour(self) -> &'static str {
        match self {
            Verdict::Kept => crate::qc::pairs::KEPT,
            Verdict::Outlier => crate::qc::pairs::REMOVED,
            Verdict::NearEmpty => "#8a8f8b",
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

    /// Barcodes with each verdict, in [`Verdict::ALL`]'s order.
    pub fn tally(&self) -> [usize; 3] {
        let mut n = [0; 3];
        for v in &self.verdict {
            n[*v as usize] += 1;
        }
        n
    }

    /// The figure, `w` by `h`: the knee on the left, genes by counts on the
    /// right, the verdicts' key and counts over them.
    pub fn figure(&self, (w, h): (f64, f64)) -> String {
        crate::figure::svg(w, h, |c| {
            let [kept, outlier, empty] = self.tally();
            let title = format!(
                "{} cells: {kept} kept, {outlier} MAD outliers, {empty} near empty, of {}",
                self.batch,
                self.counts.len()
            );
            c.text(12.0, 16.0, &title, 11.0, Anchor::Start, INK);
            let mut x = 12.0;
            for v in Verdict::ALL {
                c.rect(x, 24.0, 9.0, 9.0, v.colour());
                c.text(x + 13.0, 32.0, v.said(), 9.0, Anchor::Start, MUTED);
                x += 100.0;
            }
            let half = (w - 12.0) / 2.0;
            let note = "right: a square blue where most of its barcodes are kept, orange where not";
            c.text(12.0 + half, 32.0, note, 9.0, Anchor::Start, MUTED);
            self.knee(c, (6.0, 40.0, half - 6.0, h - 46.0));
            self.density(c, (6.0 + half, 40.0, half - 6.0, h - 46.0));
        })
    }

    /// Counts by rank, both on log scales, the curve coloured by verdict.
    fn knee(&self, c: &mut Canvas, (x, y, w, h): (f64, f64, f64, f64)) {
        let mut order: Vec<usize> = (0..self.counts.len()).collect();
        order.sort_by(|&a, &b| self.counts[b].total_cmp(&self.counts[a]));
        let n = order.len().max(2) as f64;
        let top = self.counts.iter().copied().fold(1.0f32, f32::max) as f64;
        let frame = Frame::new((x, y, w, h), (0.0, n.log10()), (0.0, top.max(10.0).log10()));
        frame.axes(c, ("barcode rank", "counts"), true);
        let mut last: Option<(f64, f64)> = None;
        for (r, &i) in order.iter().enumerate() {
            let p = frame.at(
                (r as f64 + 1.0).log10(),
                (self.counts[i] as f64).max(1.0).log10(),
            );
            // A point a pixel on from the last, or where the verdict turns.
            if let Some(q) = last {
                let turn = r > 0 && self.verdict[order[r - 1]] != self.verdict[i];
                if (p.0 - q.0).abs() < 0.6 && !turn && r + 1 < order.len() {
                    continue;
                }
                c.line(q.0, q.1, p.0, p.1, self.verdict[i].colour(), 1.6);
            }
            last = Some(p);
        }
    }

    /// Genes detected by counts (log scales), a square kept or not by
    /// which most of its barcodes are, deeper as it holds more.
    fn density(&self, c: &mut Canvas, (x, y, w, h): (f64, f64, f64, f64)) {
        const BINS: usize = 44;
        let lg = |v: f32| (v as f64).max(1.0).log10();
        // From the decade below the fewest to the most: no empty decades.
        let span = |v: &[f32]| {
            let lo = v.iter().map(|&x| lg(x)).fold(f64::INFINITY, f64::min);
            let hi = v.iter().map(|&x| lg(x)).fold(1.0, f64::max);
            (lo.floor().clamp(0.0, (hi - 1.0).max(0.0)), hi)
        };
        let (xr, yr) = (span(&self.counts), span(&self.genes));
        let frame = Frame::new((x, y, w, h), xr, yr);
        let (mut kept, mut other) = (vec![0u32; BINS * BINS], vec![0u32; BINS * BINS]);
        let slot = |v: f64, (lo, hi): (f64, f64)| {
            (((v - lo) / (hi - lo) * BINS as f64) as usize).min(BINS - 1)
        };
        for i in 0..self.counts.len() {
            let at = slot(lg(self.genes[i]), yr) * BINS + slot(lg(self.counts[i]), xr);
            if self.verdict[i] == Verdict::Kept {
                kept[at] += 1;
            } else {
                other[at] += 1;
            }
        }
        let top = (0..kept.len())
            .map(|j| kept[j] + other[j])
            .max()
            .unwrap_or(1)
            .max(1);
        let (cw, ch) = (frame.pw / BINS as f64, frame.ph / BINS as f64);
        for by in 0..BINS {
            for bx in 0..BINS {
                let (k, o) = (kept[by * BINS + bx], other[by * BINS + bx]);
                if k + o == 0 {
                    continue;
                }
                let depth = 0.18 + 0.82 * ((k + o) as f64).ln_1p() / (top as f64).ln_1p();
                let base = if k >= o {
                    Verdict::Kept
                } else {
                    Verdict::Outlier
                };
                let fill = shade(base.colour(), depth);
                let px = frame.px + bx as f64 * cw;
                let py = frame.py + frame.ph - (by + 1) as f64 * ch;
                c.rect(px, py, cw + 0.05, ch + 0.05, &fill);
            }
        }
        frame.axes(c, ("counts", "genes detected"), true);
    }
}

/// `hex` over white at `depth` (0 white, 1 itself).
fn shade(hex: &str, depth: f64) -> String {
    let p = |i: usize| u8::from_str_radix(&hex[i..i + 2], 16).unwrap_or(0) as f64;
    let m = |v: f64| (255.0 + (v - 255.0) * depth).round() as u8;
    format!("#{:02x}{:02x}{:02x}", m(p(1)), m(p(3)), m(p(5)))
}

/// A plot's box and scales.
struct Frame {
    px: f64,
    py: f64,
    pw: f64,
    ph: f64,
    xr: (f64, f64),
    yr: (f64, f64),
    x: f64,
    y: f64,
    h: f64,
}

impl Frame {
    fn new((x, y, w, h): (f64, f64, f64, f64), xr: (f64, f64), yr: (f64, f64)) -> Self {
        let (l, b) = (46.0, 30.0);
        Self {
            px: x + l,
            py: y + 6.0,
            pw: (w - l - 8.0).max(1.0),
            ph: (h - b - 6.0).max(1.0),
            xr,
            yr,
            x,
            y,
            h,
        }
    }

    fn at(&self, vx: f64, vy: f64) -> (f64, f64) {
        (
            self.px + (vx - self.xr.0) / (self.xr.1 - self.xr.0) * self.pw,
            self.py + self.ph - (vy - self.yr.0) / (self.yr.1 - self.yr.0) * self.ph,
        )
    }

    /// The axes, ticks at each power of ten (`log`: the values are
    /// log10), and their titles.
    fn axes(&self, c: &mut Canvas, (xt, yt): (&str, &str), log: bool) {
        c.line(
            self.px,
            self.py + self.ph,
            self.px + self.pw,
            self.py + self.ph,
            FAINT,
            0.8,
        );
        c.line(self.px, self.py, self.px, self.py + self.ph, FAINT, 0.8);
        let label = |v: f64| {
            let n = if log { 10f64.powf(v) } else { v };
            if n >= 1e6 {
                format!("{:.0}M", n / 1e6)
            } else if n >= 1e3 {
                format!("{:.0}k", n / 1e3)
            } else {
                format!("{n:.0}")
            }
        };
        for k in (self.xr.0.ceil() as i32)..=(self.xr.1.floor() as i32) {
            let (tx, ty) = self.at(k as f64, self.yr.0);
            c.text(tx, ty + 11.0, &label(k as f64), 8.0, Anchor::Middle, MUTED);
        }
        for k in (self.yr.0.ceil() as i32)..=(self.yr.1.floor() as i32) {
            let (tx, ty) = self.at(self.xr.0, k as f64);
            c.text(
                tx - 4.0,
                ty + 3.0,
                &label(k as f64),
                8.0,
                Anchor::End,
                MUTED,
            );
        }
        c.text(
            self.px + self.pw / 2.0,
            self.y + self.h - 2.0,
            xt,
            9.0,
            Anchor::Middle,
            INK,
        );
        c.vtext(self.x + 9.0, self.py + self.ph / 2.0, yt, 9.0, INK);
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
        let svg = c.figure((720.0, 320.0));
        assert!(svg.contains("wt_rep1 cells: 2 kept, 1 MAD outliers, 1 near empty, of 4"));
        assert!(svg.contains("barcode rank") && svg.contains("genes detected"));
        assert!(crate::figure::png(&svg, 1.0).is_ok());
    }
}
