//! The view saved as a figure: the window's rows drawn on a canvas, the
//! gene models as the figure draws them.

use super::{fmt_thousands, Bins, PileupView, Row, GENE_BAND, GENE_LANE};
use crate::figure::{self, Anchor, Bars, Canvas, INK, MUTED};
use crate::site_analysis::miami::genemodel::draw_gene_model;

impl PileupView<'_> {
    /// The visible window as a figure, one panel per track.
    pub(super) fn figure(&self) -> String {
        let (lo, hi) = self.window;
        // Turned position labels take room under the axis: the bars keep theirs.
        let room = figure::tick_room(&self.tick_list(), true);
        let rows = self.rows();
        let genes_h = 24.0 + GENE_LANE * self.lanes() as f64;
        let n_plots = Self::plot_count(&rows);
        // The shape the view has on screen, grown alike on both sides until
        // each plot has room for its axis, the plots sharing what the title
        // and the genes leave; else a page of its own.
        let genes = if rows.contains(&Row::Genes) {
            genes_h
        } else {
            0.0
        };
        // The plots are alike: each at least `least`, sharing alike what a
        // taller page adds.
        let (fixed, least) = (50.0 + genes, 100.0 + room);
        let page = self
            .controls
            .page(self.drawn, (360.0, fixed + least * n_plots as f64));
        let panel_h = page.map_or(170.0 + room, |(_, h)| (h - fixed) / n_plots.max(1) as f64);
        let heights: Vec<f64> = rows
            .iter()
            .map(|&r| if r == Row::Genes { genes_h } else { panel_h })
            .collect();
        let (width, height) = page.unwrap_or((720.0, 50.0 + heights.iter().sum::<f64>()));
        let mut c = Canvas::new(width, height);
        c.bold(16.0, 22.0, &self.title, 12.0, Anchor::Start, INK);
        c.text(
            16.0,
            38.0,
            &format!(
                "{}:{}-{}, {} bp per bar",
                self.chr,
                fmt_thousands(lo),
                fmt_thousands(hi),
                self.bin_width()
            ),
            9.0,
            Anchor::Start,
            MUTED,
        );
        let bins = self.bins();
        let mut y = 44.0;
        for (&row, &h) in rows.iter().zip(&heights) {
            let title = self.row_title(row);
            self.draw_row(row, &bins, &mut c, (8.0, y, width - 16.0, h), None, title);
            y += h;
        }

        c.finish()
    }

    /// Draw one row (track, contrast, mirror or genes) in the box.
    pub(super) fn draw_row(
        &self,
        row: Row,
        bins: &Bins,
        c: &mut Canvas,
        bbox: (f64, f64, f64, f64),
        pointer: Option<usize>,
        title: String,
    ) {
        let (x, y, w, h) = bbox;
        let x_title = format!("{} position", self.chr);
        match row {
            Row::Track(i) => {
                let t = &self.tracks[i];
                let (lo, hi) = self.window;
                let shown = &bins.shown[i];
                let marks = t.sites_in(lo, hi).into_iter();
                Bars {
                    rotated: true,
                    values: &shown.values,
                    front: shown.front.as_deref(),
                    accent: &|_| shown.accent,
                    colour: &|_| figure::BAR,
                    y_scale: self.y_scale,
                    y_max: self.row_max(row, bins),
                    ticks: self.tick_list(),
                    pointer,
                    marks: marks.map(|s| bins.edges.col_of(s)).collect(),
                    dividers: Vec::new(),
                    title,
                    x_title,
                    y_title: self.y_title(i),
                }
                .draw(c, x, y, w, h);
            }
            Row::Contrast(contrast) => {
                let values = self.contrast_values(contrast, bins);
                let measure = contrast.measure;
                figure::Diverging {
                    rotated: true,
                    values: &values,
                    ticks: self.tick_list(),
                    pointer,
                    title,
                    x_title,
                    y_title: measure.name(self.on),
                    label: &|v| measure.label(v),
                }
                .draw(c, x, y, w, h);
            }
            Row::Mirror(a, b) => {
                let half = |k: usize| {
                    let shown = &bins.shown[k];
                    // Converted alone draws all in the converted colour.
                    let front = match (&shown.front, shown.accent) {
                        (Some(f), _) => Some(f.as_slice()),
                        (None, true) => Some(shown.values.as_slice()),
                        (None, false) => None,
                    };
                    figure::Half {
                        values: &shown.values,
                        front,
                        colour: figure::BAR,
                        name: self.tracks[k].label.to_string(),
                    }
                };
                figure::Mirror {
                    rotated: true,
                    up: half(a),
                    down: half(b),
                    y_scale: self.y_scale,
                    y_max: self.row_max(row, bins),
                    y_labels: Default::default(),
                    ticks: self.tick_list(),
                    pointer,
                    title,
                    x_title,
                    y_title: self.y_title(a),
                }
                .draw(c, x, y, w, h);
            }
            Row::Genes => self.draw_genes(c, bbox, title),
        }
    }

    /// About five labelled genomic ticks across the window.
    fn tick_list(&self) -> Vec<(usize, String)> {
        let every = (self.columns / 5).max(1);
        (0..self.columns)
            .step_by(every)
            .map(|k| (k, fmt_thousands(self.bar_range(k).0)))
            .collect()
    }

    /// The genes row in the box, as the figure draws gene models.
    fn draw_genes(&self, c: &mut Canvas, (x, y, w, _h): (f64, f64, f64, f64), title: String) {
        let (left, right) = (52.0, 12.0);
        c.bold(x + left, y + 13.0, &title, 11.0, Anchor::Start, INK);
        let edges = self.edges();
        for (lane, g) in self.gene_lanes() {
            let mid = y + 26.0 + GENE_LANE * lane as f64;
            draw_gene_model(c, g, &edges, x + left, w - left - right, mid, GENE_BAND);
        }
    }
}
