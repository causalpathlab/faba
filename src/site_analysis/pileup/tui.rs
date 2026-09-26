//! Full-screen gene-body browser for `faba pileup --interactive`.
//!
//! The matrix track, and the site track when a site table was given, share
//! one genomic axis. Every frame re-bins the raw positions of the visible
//! window, one bar per terminal column, with the same rule as the printed
//! pileup ([`BinEdges`]), so zooming in resolves single sites. The cursor is a
//! genomic coordinate, so it stays put through zooms and resizes.

use data_beans::interactive::ui::{
    header, help_line, input_line, panel, run_screen, Binning, HistPlot, Scale, Screen, DIM,
    HIGHLIGHT, PLAIN,
};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::figure::term::{self, PlotImage};
use crate::figure::{self, Anchor, Bars, Canvas, SavePrompt, INK, MUTED};
use ratatui_image::picker::Picker;

use genomic_data::coordinates::chr_eq;

use super::{fmt_thousands, BinEdges};

/// Width of HistPlot's y gutter.
const GUTTER: u16 = 6;

/// Columns between labelled ticks, at least.
const TICK_SPACING: usize = 14;

/// One track: raw `(position, value)` pairs, sorted by position.
pub struct Track<'a> {
    pub label: &'a str,
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

/// How the user left the browser.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Exit {
    Quit,
    /// `g`: back to the gene list.
    Genes,
    /// `/`: a gene or locus the view cannot show itself.
    Search(String),
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
    save: SavePrompt,
    /// Whether `g` goes back to a gene list.
    back: bool,
    /// The `/` search being typed.
    search: Option<String>,
    /// A message for the footer, until the next key.
    status: Option<String>,
    /// The terminal's image support; `use_images` draws the tracks with it.
    images: Option<Picker>,
    use_images: bool,
    plots: Vec<PlotImage>,
    exit: Option<Exit>,
}

impl<'a> PileupView<'a> {
    pub fn new(
        title: &str,
        chr: &str,
        tracks: Vec<Track<'a>>,
        extent: (i64, i64),
        back: bool,
    ) -> Self {
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
            save: SavePrompt::new(&format!("pileup_{title}")),
            back,
            search: None,
            status: None,
            images: None,
            use_images: false,
            plots: Vec::new(),
            exit: None,
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

    /// Act on a `/` search: a locus inside this view moves the window
    /// there; anything else leaves the view for the caller to open.
    fn submit(&mut self, query: &str) {
        match super::parse_query(query) {
            Some(super::Query::Locus(r, single)) if chr_eq(&r.chr, &self.chr) => {
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

    /// Show `msg` on the footer until the next key.
    pub fn set_status(&mut self, msg: Option<String>) {
        self.status = msg;
    }

    fn cursor_col(&self) -> usize {
        self.edges().col_of(self.cursor)
    }

    /// The visible window as a figure, one panel per track.
    fn figure(&self) -> String {
        let (lo, hi) = self.window;
        let panel_h = 170.0;
        let mut c = Canvas::new(720.0, 50.0 + panel_h * self.tracks.len() as f64);
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
        for (i, t) in self.tracks.iter().enumerate() {
            let title = format!("{} ({})", t.label, t.signal);
            let panel = (8.0, 44.0 + panel_h * i as f64, 704.0, panel_h);
            self.draw_track(i, &mut c, panel, None, title);
        }
        c.finish()
    }

    /// Track `i` over the visible window in the box `(x, y, w, h)`.
    fn draw_track(
        &self,
        i: usize,
        c: &mut Canvas,
        (x, y, w, h): (f64, f64, f64, f64),
        pointer: Option<usize>,
        title: String,
    ) {
        let t = &self.tracks[i];
        let edges = self.edges();
        let (lo, hi) = self.window;
        let every = (self.columns / 5).max(1);
        let ticks = (0..self.columns)
            .step_by(every)
            .map(|k| (k, fmt_thousands(self.bar_range(k).0)))
            .collect();
        let values = t.bin(&edges);
        let marks = t
            .sites_in(lo, hi)
            .into_iter()
            .map(|s| edges.col_of(s))
            .collect();
        Bars {
            values: &values,
            front: None,
            accent: &|_| false,
            y_scale: self.y_scale,
            ticks,
            pointer,
            marks,
            title,
            x_title: format!("{} position", self.chr),
            y_title: t.signal.into(),
        }
        .draw(c, x, y, w, h);
    }

    /// Track `i` with the cursor, sized for an on-screen area.
    fn track_svg(&self, i: usize, w: f64, h: f64) -> String {
        let mut c = Canvas::new(w, h);
        let pointer = Some(self.cursor_col());
        self.draw_track(i, &mut c, (0.0, 0.0, w, h), pointer, String::new());
        c.finish()
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
        self.exit.is_some()
    }

    fn interrupt(&mut self) {
        self.exit = Some(Exit::Quit);
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.plots.iter_mut().for_each(PlotImage::invalidate);
        if self.save.active() {
            if let Some(prefix) = self.save.handle(key) {
                let result = figure::save(&self.figure(), &prefix);
                self.save.report(result);
            }
            return;
        }
        self.save.dismiss();
        self.status = None;
        if let Some(buf) = &mut self.search {
            match key.code {
                KeyCode::Char(ch) if buf.len() < 128 => buf.push(ch),
                KeyCode::Backspace => {
                    buf.pop();
                }
                KeyCode::Esc => self.search = None,
                KeyCode::Enter => {
                    let query = self.search.take().unwrap_or_default();
                    self.submit(query.trim());
                }
                _ => {}
            }
            return;
        }
        match key.code {
            KeyCode::Char('/') => self.search = Some(String::new()),
            KeyCode::Char('s') => self.save.open(),
            KeyCode::Char('i') if self.images.is_some() => self.use_images ^= true,
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
            KeyCode::Char('g') if self.back => self.exit = Some(Exit::Genes),
            KeyCode::Char('q') | KeyCode::Esc | KeyCode::Enter => self.exit = Some(Exit::Quit),
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
        let mut extra = format!(
            "{}:{}-{}  {} bp/bar · y {}",
            self.chr,
            fmt_thousands(lo),
            fmt_thousands(hi),
            self.bin_width(),
            self.y_scale.name()
        );
        if self.images.is_some() {
            extra += if self.use_images {
                " · image"
            } else {
                " · text"
            };
        }
        frame.render_widget(header("pileup", &self.title, &extra), top);

        let every = TICK_SPACING.max(self.columns / 5);
        let ticks: Vec<Option<String>> = (0..self.columns)
            .map(|k| {
                (k.is_multiple_of(every) && k + 10 < self.columns)
                    .then(|| fmt_thousands(self.bar_range(k).0))
            })
            .collect();
        let label = |k: i32| ticks.get(k as usize).cloned().flatten();
        let pointer = Some(self.cursor_col() as i32);
        let picker = self.images.clone().filter(|_| self.use_images);
        let mut plots = std::mem::take(&mut self.plots);
        plots.resize_with(self.tracks.len(), PlotImage::default);
        for (i, &area) in areas[1..areas.len() - 2].iter().enumerate() {
            let t = &self.tracks[i];
            let block = panel(format!(" {} · {} ", t.label, t.signal), true);
            let inner = block.inner(area);
            frame.render_widget(block, area);
            if let Some(picker) = &picker {
                let this = &*self;
                if plots[i].render(frame, inner, picker, |w, h| this.track_svg(i, w, h)) {
                    continue;
                }
            }
            let t = &self.tracks[i];
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
        self.plots = plots;

        frame.render_widget(Paragraph::new(self.readout()), readout);
        let help = if let Some(buf) = &self.search {
            input_line(
                "search gene or chr:start-end: ",
                buf,
                &[("Enter", "go"), ("Esc", "back")],
            )
        } else if let Some(msg) = &self.status {
            Line::from(Span::styled(format!(" {msg}"), HIGHLIGHT))
        } else {
            self.save.footer().unwrap_or_else(|| self.help())
        };
        frame.render_widget(help, footer);
    }
}

impl PileupView<'_> {
    fn help(&self) -> Line<'static> {
        let mut keys = vec![
            ("←/→", "bar"),
            ("n/p", "next/prev site"),
            ("+/-", "zoom"),
            ("0", "whole"),
            ("/", "search"),
            ("y", "scale"),
            ("s", "save"),
        ];
        if self.images.is_some() {
            keys.push(("i", "image/text"));
        }
        if self.back {
            keys.push(("g", "genes"));
        }
        keys.push(("q", "quit"));
        help_line(&keys)
    }
}

/// Browse full screen until the user quits, asks for the gene list (with
/// `back`), or searches for something outside the view. `status` shows on
/// the footer first.
pub fn show_pileup(
    title: &str,
    chr: &str,
    tracks: Vec<Track>,
    extent: (i64, i64),
    back: bool,
    status: Option<String>,
) -> anyhow::Result<Exit> {
    let mut view = PileupView::new(title, chr, tracks, extent, back);
    view.set_status(status);
    view.images = term::picker();
    view.use_images = view.images.is_some();
    run_screen(&mut view)?;
    Ok(view.exit.unwrap_or(Exit::Quit))
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
