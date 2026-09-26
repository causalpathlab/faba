//! Full-screen gene-body browser for `faba pileup --interactive`.
//!
//! The matrix track, and the site track when a site table was given, share
//! one genomic axis. Every frame re-bins the raw positions of the visible
//! window, one bar per terminal column, with the same rule as the printed
//! pileup ([`BinEdges`]), so zooming in resolves single sites. The cursor is a
//! genomic coordinate, so it stays put through zooms and resizes.

use data_beans::interactive::ui::{
    header, help_line, panel, run_screen, Binning, HistPlot, Scale, Screen, DIM, HIGHLIGHT, PLAIN,
};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use super::{fmt_thousands, BinEdges};

/// Width of HistPlot's y gutter.
const GUTTER: u16 = 6;

/// Columns between labelled ticks, at least.
const TICK_SPACING: usize = 14;

/// One track: raw `(position, value)` pairs, sorted by position.
pub struct Track<'a> {
    pub label: &'static str,
    pub signal: &'static str,
    pub positions: &'a [(i64, f64)],
    /// Bins become `log10(1 + sum)`, as the printed pileup does.
    pub log: bool,
}

impl Track<'_> {
    fn bin(&self, edges: &BinEdges) -> Vec<f64> {
        edges.bin(self.positions, self.log)
    }

    /// Distinct positions inside `lo..=hi`.
    fn sites_in(&self, lo: i64, hi: i64) -> Vec<i64> {
        let a = self.positions.partition_point(|p| p.0 < lo);
        let b = self.positions.partition_point(|p| p.0 <= hi);
        let mut out: Vec<i64> = self.positions[a..b].iter().map(|p| p.0).collect();
        out.dedup();
        out
    }
}

/// State of the browser, independent of the terminal so it can be tested.
pub struct PileupView<'a> {
    title: String,
    chr: String,
    tracks: Vec<Track<'a>>,
    /// Every distinct site position over all tracks, for n/p jumps.
    sites: Vec<i64>,
    extent: (i64, i64),
    /// The visible window, inclusive.
    window: (i64, i64),
    cursor: i64,
    /// Bars per track in the last frame (one per chart column).
    columns: usize,
    y_scale: Scale,
    done: bool,
}

impl<'a> PileupView<'a> {
    pub fn new(title: &str, chr: &str, tracks: Vec<Track<'a>>, extent: (i64, i64)) -> Self {
        let mut sites: Vec<i64> = tracks
            .iter()
            .flat_map(|t| t.positions.iter().map(|p| p.0))
            .collect();
        sites.sort_unstable();
        sites.dedup();
        let (lo, hi) = (extent.0.min(extent.1), extent.0.max(extent.1));
        Self {
            title: title.to_string(),
            chr: chr.to_string(),
            tracks,
            cursor: sites.first().copied().unwrap_or(lo).clamp(lo, hi),
            sites,
            extent: (lo, hi),
            window: (lo, hi),
            columns: 80,
            y_scale: Scale::Linear,
            done: false,
        }
    }

    fn edges(&self) -> BinEdges {
        BinEdges::new(self.window.0, self.window.1, self.columns)
    }

    /// Bases per bar at the current zoom.
    fn bin_width(&self) -> i64 {
        (self.edges().span().div_ceil(self.columns.max(1) as u64) as i64).max(1)
    }

    /// Genomic range `[start, stop)` of bar `col`.
    fn bar_range(&self, col: usize) -> (i64, i64) {
        let e = self.edges();
        let span = e.span() as i64;
        let n = e.num_bins as i64;
        let at = |i: i64| e.min_pos + i * span / n;
        (at(col as i64), at(col as i64 + 1))
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

    fn clamp_window(&mut self) {
        let (flo, fhi) = self.extent;
        let w = (self.window.1 - self.window.0).min(fhi - flo);
        let lo = self.window.0.clamp(flo, fhi - w);
        self.window = (lo, lo + w);
    }

    fn move_cursor(&mut self, bars: i64) {
        let (lo, hi) = self.extent;
        self.cursor = (self.cursor + bars * self.bin_width()).clamp(lo, hi);
        self.follow_cursor();
    }

    fn jump_site(&mut self, forward: bool) {
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

    /// Zoom by `factor` (< 1 in, > 1 out) around the cursor; a bar never
    /// gets narrower than one base.
    fn zoom(&mut self, factor: f64) {
        let (lo, hi) = self.window;
        let min_w = self.columns.max(1) as i64;
        let w = (((hi - lo) as f64 * factor).round() as i64)
            .max(min_w)
            .min(self.extent.1 - self.extent.0);
        let left = ((self.cursor - lo) as f64 / (hi - lo).max(1) as f64 * w as f64) as i64;
        self.window = (self.cursor - left, self.cursor - left + w);
        self.clamp_window();
    }

    fn cursor_col(&self) -> usize {
        self.edges().col_of(self.cursor)
    }

    fn readout(&self) -> Line<'static> {
        let dim = |t: String| Span::styled(t, DIM);
        let col = self.cursor_col();
        let (start, stop) = self.bar_range(col);
        let mut spans = vec![
            Span::styled(
                format!(" {}:{}", self.chr, fmt_thousands(self.cursor)),
                HIGHLIGHT,
            ),
            dim(format!(
                "   bar {}-{}",
                fmt_thousands(start),
                fmt_thousands(stop)
            )),
        ];
        let edges = self.edges();
        for t in &self.tracks {
            let v = t.bin(&edges)[col];
            spans.push(dim(format!("   {} ", t.label)));
            spans.push(Span::raw(format!("{v:.2}")));
        }
        let here = self
            .sites
            .iter()
            .filter(|&&s| s >= start && s < stop.max(start + 1))
            .count();
        spans.push(dim(format!("   {here} site(s) in bar")));
        Line::from(spans)
    }
}

impl Screen for PileupView<'_> {
    fn done(&self) -> bool {
        self.done
    }

    fn interrupt(&mut self) {
        self.done = true;
    }

    fn handle_key(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Left | KeyCode::Char('h') => self.move_cursor(-1),
            KeyCode::Right | KeyCode::Char('l') => self.move_cursor(1),
            KeyCode::PageUp => self.move_cursor(-(self.columns as i64) / 2),
            KeyCode::PageDown => self.move_cursor(self.columns as i64 / 2),
            KeyCode::Char('n') => self.jump_site(true),
            KeyCode::Char('p') => self.jump_site(false),
            KeyCode::Char('+' | '=') => self.zoom(0.5),
            KeyCode::Char('-') => self.zoom(2.0),
            KeyCode::Char('0') => self.window = self.extent,
            KeyCode::Char('y') => self.y_scale = self.y_scale.next(),
            KeyCode::Char('q') | KeyCode::Esc | KeyCode::Enter => self.done = true,
            _ => {}
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        let n_tracks = self.tracks.len().max(1) as u32;
        let mut rows = vec![Constraint::Length(1)];
        rows.extend((0..n_tracks).map(|_| Constraint::Ratio(1, n_tracks)));
        rows.extend([Constraint::Length(1), Constraint::Length(1)]);
        let areas = Layout::vertical(rows).split(frame.area());
        let (top, readout, footer) = (areas[0], areas[areas.len() - 2], areas[areas.len() - 1]);

        // One bar per chart column, re-binned from the raw positions.
        let inner_width = areas[1].width.saturating_sub(2 + GUTTER);
        self.columns = (inner_width as usize).max(1);
        self.clamp_window();
        let edges = self.edges();
        let (lo, hi) = self.window;
        let extra = format!(
            "{}:{}-{}  {} bp/bar · y {}",
            self.chr,
            fmt_thousands(lo),
            fmt_thousands(hi),
            self.bin_width(),
            self.y_scale.name()
        );
        frame.render_widget(header("pileup", &self.title, &extra), top);

        let every = TICK_SPACING.max(self.columns / 5);
        let label = |k: i32| {
            let k = k as usize;
            (k.is_multiple_of(every) && k + 10 < self.columns)
                .then(|| fmt_thousands(self.bar_range(k).0))
        };
        let pointer = Some(self.cursor_col() as i32);
        for (t, &area) in self.tracks.iter().zip(&areas[1..areas.len() - 2]) {
            let block = panel(format!(" {} · {} ", t.label, t.signal), true);
            let inner = block.inner(area);
            frame.render_widget(block, area);
            let values = t.bin(&edges);
            let marks = t
                .sites_in(lo, hi)
                .into_iter()
                .map(|s| (edges.col_of(s) as i32, "+", DIM))
                .collect();
            HistPlot {
                bins: Binning::with_width(Scale::Linear, 1.0),
                kmin: 0,
                counts: &values,
                style: &|_| PLAIN,
                subset: None,
                y_scale: self.y_scale,
                pointer,
                marks,
                x_label: Some(&label),
                tick_every: Some(1),
            }
            .render(frame.buffer_mut(), inner);
        }

        frame.render_widget(Paragraph::new(self.readout()), readout);
        frame.render_widget(
            help_line(&[
                ("←/→", "bar"),
                ("n/p", "next/prev site"),
                ("+/-", "zoom"),
                ("0", "whole"),
                ("y", "scale"),
                ("q", "quit"),
            ]),
            footer,
        );
    }
}

/// Browse full screen until the user quits.
pub fn show_pileup(
    title: &str,
    chr: &str,
    tracks: Vec<Track>,
    extent: (i64, i64),
) -> anyhow::Result<()> {
    run_screen(&mut PileupView::new(title, chr, tracks, extent))
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
