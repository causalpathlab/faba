//! Rows that do not all fit: each shown keeps room to draw, the contrast
//! and the first track's row stay on top, and the rows under them scroll
//! with `[` and `]`, those hidden counted in the titles. The first shown is
//! kept as a track, so it is the same split or mirrored.

use super::{PileupView, Row};
use crate::figure::term::PlotImage;
use crate::tui::{hidden_marks, moved};

/// Lines a plot row needs to draw: a mirrored row's bar above and below
/// its zero line, the x axis and its labels, and the frame.
const MIN_ROW: u16 = 7;

/// Where the last fit left the rows that scroll.
#[derive(Debug, Clone, Copy, Default, PartialEq, Eq)]
pub(super) struct Fit {
    /// How many of them show.
    shown: usize,
    /// The first shown, among them.
    at: usize,
    /// How many there are.
    of: usize,
}

/// The first and last track a row draws.
fn tracks(row: Row) -> Option<(usize, usize)> {
    match row {
        Row::Track(k) => Some((k, k)),
        Row::Mirror(a, b) => Some((a, b)),
        Row::Contrast(_) | Row::Genes => None,
    }
}

impl PileupView<'_> {
    /// How many of `rows` are plots, and how many of those stay on top.
    fn plots_fixed(&self, rows: &[Row]) -> (usize, usize) {
        let plots = rows.iter().filter(|&&r| r != Row::Genes).count();
        (plots, (usize::from(self.contrast.is_some()) + 1).min(plots))
    }

    /// The rows that scroll.
    fn scrolled_rows(&self) -> Vec<Row> {
        let rows = self.rows();
        let (plots, fixed) = self.plots_fixed(&rows);
        rows[fixed..plots].to_vec()
    }

    /// Where in `scrolled` the first shown is: the first row holding the
    /// track scrolled to or one after it, no further than the last rows
    /// filling the room.
    fn scrolled_at(&self, scrolled: &[Row]) -> usize {
        let holds = |r: &Row| tracks(*r).is_some_and(|(_, last)| last >= self.track_from);
        let at = scrolled.iter().position(holds).unwrap_or(scrolled.len());
        at.min(scrolled.len().saturating_sub(self.fit.shown))
    }

    /// Scrolled to the first track of scrolled row `at`, if there is one.
    fn scroll_to(&mut self, scrolled: &[Row], at: usize) {
        if let Some((first, _)) = scrolled.get(at).copied().and_then(tracks) {
            self.track_from = first;
        }
    }

    /// Room for the plot rows of `rows` in `room` lines: all when each
    /// gets [`MIN_ROW`], else as many as do, at least one; the first shown
    /// kept within reach of the last.
    pub(super) fn fit_rows(&mut self, rows: &[Row], room: u16) {
        let (plots, fixed) = self.plots_fixed(rows);
        let fit = if plots as u32 * MIN_ROW as u32 <= room as u32 {
            plots
        } else {
            (room / MIN_ROW).max(1) as usize
        };
        let of = plots - fixed;
        self.fit.shown = fit.saturating_sub(fixed).min(of);
        self.fit.at = self.scrolled_at(&rows[fixed..plots]);
        self.fit.of = of;
    }

    /// Whether the last fit left rows out.
    pub(super) fn rows_hidden(&self) -> bool {
        self.fit.shown < self.fit.of
    }

    /// Of `rows`, those the last fit shows: those on top, the scrolled ones
    /// from the first shown, then the genes.
    pub(super) fn shown_rows(&self, rows: &[Row]) -> Vec<Row> {
        let (plots, fixed) = self.plots_fixed(rows);
        let from = (fixed + self.fit.at).min(plots);
        let to = (from + self.fit.shown).min(plots);
        [&rows[..fixed], &rows[from..to], &rows[plots..]].concat()
    }

    /// The title of shown row `k` (of `shown`), with the rows hidden above
    /// counted before the first scrolled row's and those below after the
    /// last's; on the last row on top when none of them is shown.
    pub(super) fn marked_title(&self, k: usize, shown: &[Row], title: String) -> String {
        let (plots, fixed) = self.plots_fixed(shown);
        let Fit { shown: fit, at, of } = self.fit;
        let below = of - (at + fit).min(of);
        let first = if plots > fixed {
            fixed
        } else {
            plots.saturating_sub(1)
        };
        let mut parts = Vec::new();
        if k == first && at > 0 {
            parts.push(hidden_marks(at, 0));
        }
        parts.push(title);
        if k + 1 == plots && below > 0 {
            parts.push(hidden_marks(0, below));
        }
        parts.join(" ")
    }

    /// `[` and `]`: the scrolled rows `by` places along, no further than
    /// the last ones filling the room.
    /// The plots are drawn anew only when that moves them.
    pub(super) fn scroll_rows(&mut self, by: isize) {
        let scrolled = self.scrolled_rows();
        let last = scrolled.len().saturating_sub(self.fit.shown);
        let at = moved(self.scrolled_at(&scrolled), by, last);
        let was = self.track_from;
        self.scroll_to(&scrolled, at);
        if self.track_from != was {
            self.plots.iter_mut().for_each(PlotImage::invalidate);
        }
    }

    /// Mirrored rows (`on`) or a row a track, the first row shown still
    /// the one holding the same track.
    pub(super) fn set_mirror(&mut self, on: bool) {
        if on == self.mirror {
            return;
        }
        let scrolled = self.scrolled_rows();
        self.scroll_to(&scrolled, self.scrolled_at(&scrolled));
        self.mirror = on;
    }
}
