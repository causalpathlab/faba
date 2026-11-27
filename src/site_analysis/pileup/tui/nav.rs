//! Moving about the view: the bars' edges, the cursor, panning, zooming
//! and the `/` search's loci.

use super::{BinEdges, Exit, PileupView, MAX_SPAN};
use crate::figure::term::PlotImage;
use crate::site_analysis::pileup::{parse_query, window_at_sites, Query};
use crate::tui::ui::Scale;
use genomic_data::coordinates::chr_eq;

/// What a figure's window spans.
#[derive(clap::ValueEnum, Clone, Copy, Debug, PartialEq, Eq)]
pub enum Window {
    /// The whole gene, TSS to TES, as the browser opens.
    Gene,
    /// At the gene's sites, of any track: centred on them, 1 kb at least,
    /// or padded as asked (`0` in the browser).
    Sites,
}

impl Window {
    pub fn name(self) -> &'static str {
        match self {
            Window::Gene => "gene",
            Window::Sites => "sites",
        }
    }
}

impl PileupView<'_> {
    /// The first and last site of any track (a read or its total), when
    /// there are any.
    fn sites_span(&self) -> Option<(i64, i64)> {
        let reads = self.tracks.iter().flat_map(|t| {
            let behind = t.behind().unwrap_or_default();
            t.front().iter().chain(behind)
        });
        let at = reads.filter(|p| p.1 > 0.0).map(|p| p.0);
        at.fold(None, |span, p| match span {
            None => Some((p, p)),
            Some((lo, hi)) => Some((lo.min(p), hi.max(p))),
        })
    }

    /// The window at the sites: from the first to the last and the pad
    /// either side, else around them (see [`around`]); within the extent,
    /// a base a bar at least; the whole extent when there are none.
    pub(super) fn sites_window(&self) -> (i64, i64) {
        match self.sites_span() {
            Some(sites) => window_at_sites(sites, self.pad, self.extent, self.columns),
            None => self.extent,
        }
    }

    /// The window `w` asks for: the whole gene or its sites; the cursor
    /// kept in it.
    pub(super) fn set_window(&mut self, w: Window) {
        self.window = match w {
            Window::Gene => self.extent,
            Window::Sites => self.sites_window(),
        };
        self.cursor = self.cursor.clamp(self.window.0, self.window.1);
        self.plots.iter_mut().for_each(PlotImage::invalidate);
    }

    /// Which of the two the window is, when it is either.
    pub(super) fn window_is(&self) -> Option<Window> {
        match self.window {
            w if w == self.extent => Some(Window::Gene),
            w if w == self.sites_window() => Some(Window::Sites),
            _ => None,
        }
    }

    /// `0`: the whole gene, or from the whole gene its sites.
    pub(super) fn whole_or_sites(&mut self) {
        match self.window == self.extent {
            true => self.set_window(Window::Sites),
            false => self.set_window(Window::Gene),
        }
    }

    pub(super) fn edges(&self) -> BinEdges {
        BinEdges::new(self.window.0, self.window.1, self.columns).on(self.x_scale)
    }

    /// What a bar spans, said: its bases ("120 bp per bar"), or on a sqrt
    /// or log x scale that scale ("x: sqrt of bp from the left").
    pub(super) fn per_bar(&self) -> String {
        match self.x_scale {
            Scale::Linear => format!("{} bp per bar", self.bin_width()),
            scale => format!("x: {} of bp from the left", scale.name()),
        }
    }

    /// Bases per bar at the current zoom.
    pub(super) fn bin_width(&self) -> i64 {
        (self.edges().span().div_ceil(self.columns.max(1) as u64) as i64).max(1)
    }

    /// Genomic range `[start, stop)` of bar `col`.
    pub(super) fn bar_range(&self, col: usize) -> (i64, i64) {
        self.edges().col_range(col)
    }

    /// Slide the window so it holds the cursor, keeping its width.
    fn follow_cursor(&mut self) {
        let (lo, hi) = self.window;
        let w = hi - lo;
        if self.cursor < lo {
            self.window = (self.cursor, self.cursor + w);
        } else if self.cursor > hi {
            self.window = (self.cursor - w, self.cursor);
        }
        self.clamp_window();
    }

    pub(super) fn clamp_window(&mut self) {
        let (flo, fhi) = self.extent;
        let w = (self.window.1 - self.window.0).min(fhi - flo);
        let lo = self.window.0.clamp(flo, fhi - w);
        self.window = (lo, lo + w);
    }

    /// Slide the window `bars` bars along, keeping its width and the cursor
    /// in it.
    /// Whether it moved: at an end of the view's span it cannot.
    pub(super) fn pan(&mut self, bars: i64) -> bool {
        // As far as the bars travel: the window's span a column, not the
        // whole bases a bar rounds up to.
        let span = self.window.1 - self.window.0;
        let shift = bars * span / self.columns.max(1) as i64;
        let before = self.window;
        self.window = (self.window.0 + shift, self.window.1 + shift);
        self.clamp_window();
        if self.window == before {
            return false;
        }
        self.cursor = self.cursor.clamp(self.window.0, self.window.1);
        self.plots.iter_mut().for_each(PlotImage::invalidate);
        true
    }

    /// Whether screen cell `(x, y)` is over the plots' bars.
    pub(super) fn over_bars(&self, x: u16, y: u16) -> bool {
        let (top, bottom) = self.chart_rows;
        self.bar_at(x).is_some() && (top..bottom).contains(&y)
    }

    /// The bar at screen column `x`, when it is over the bars.
    pub(super) fn bar_at(&self, x: u16) -> Option<usize> {
        let col = usize::from(x.checked_sub(self.chart_x)?);
        (col < self.columns).then_some(col)
    }

    pub(super) fn move_cursor(&mut self, bars: i64) {
        let (lo, hi) = self.extent;
        // On a sqrt or log x scale the bars differ in bases: to the bar
        // that many over, while the cursor is in the window.
        let edges = self.edges();
        let (wlo, whi) = self.window;
        let step = match edges.is_linear() || !(wlo..=whi).contains(&self.cursor) {
            true => bars * self.bin_width(),
            false => {
                // The next bar that way holding a base (one narrower than a
                // base, at the spread left end, holds none).
                let from = edges.col_of(self.cursor) as i64;
                let held = (bars.abs()..).map(|k| from + bars.signum() * k);
                let next = (held.take_while(|&c| (0..self.columns as i64).contains(&c)))
                    .map(|c| edges.col_range(c as usize))
                    .find(|&(start, stop)| stop > start);
                match next {
                    Some((start, stop)) => start + (stop - start - 1) / 2 - self.cursor,
                    None => bars * self.bin_width(),
                }
            }
        };
        self.cursor = (self.cursor + step).clamp(lo, hi);
        self.follow_cursor();
    }

    pub(super) fn jump_site(&mut self, forward: bool) {
        let next = if forward {
            self.sites.iter().find(|&&s| s > self.cursor)
        } else {
            self.sites.iter().rev().find(|&&s| s < self.cursor)
        };
        if let Some(&s) = next {
            self.cursor = s;
            self.follow_cursor();
        }
    }

    /// Zoom by `factor` (< 1 in, > 1 out) around the cursor; a bar spans at least one base.
    pub(super) fn zoom(&mut self, factor: f64) {
        let (lo, hi) = self.window;
        let min_w = self.columns.max(1) as i64;
        let w = (((hi - lo) as f64 * factor).round() as i64)
            .max(min_w)
            .min(self.extent.1 - self.extent.0);
        let left = ((self.cursor - lo) as f64 / (hi - lo).max(1) as f64 * w as f64) as i64;
        self.window = (self.cursor - left, self.cursor - left + w);
        self.clamp_window();
    }

    /// Zoom out twice as wide; with all of the view's own span in view,
    /// ask the caller for twice that span around it, flanks and their
    /// genes included, up to [`MAX_SPAN`].
    pub(super) fn zoom_out(&mut self) {
        let (lo, hi) = self.extent;
        if self.window != self.extent {
            return self.zoom(2.0);
        }
        if hi - lo >= MAX_SPAN {
            self.status = Some("zoomed out as far as it goes".into());
            return;
        }
        let pad = ((hi - lo) / 2).max(self.columns as i64);
        self.exit = Some(Exit::Locus {
            lo: (lo - pad).max(1),
            hi: hi + pad,
            at: self.cursor,
        });
    }

    /// A `/` search: a locus in this view moves there; anything else exits to the caller.
    pub(super) fn submit(&mut self, query: &str) {
        match parse_query(query) {
            Some(Query::Locus(r, single)) if chr_eq(&r.chr, &self.chr) => {
                let (flo, fhi) = self.extent;
                if r.ub < flo || r.lb > fhi {
                    self.exit = Some(Exit::Search(query.to_string()));
                } else if single {
                    let w = self.window.1 - self.window.0;
                    self.cursor = r.lb.clamp(flo, fhi);
                    self.window = (self.cursor - w / 2, self.cursor - w / 2 + w);
                    self.clamp_window();
                } else {
                    let (lo, hi) = (r.lb.max(flo), r.ub.min(fhi));
                    let min_w = self.columns.max(1) as i64;
                    let pad = (min_w - (hi - lo)).max(0) / 2;
                    self.window = (lo - pad, hi + pad);
                    self.cursor = (lo + hi) / 2;
                    self.clamp_window();
                }
            }
            Some(_) => self.exit = Some(Exit::Search(query.to_string())),
            None => {}
        }
    }

    pub(super) fn cursor_col(&self) -> usize {
        self.edges().col_of(self.cursor)
    }
}
