//! Full-screen metagene profile for `faba metagene --interactive`.
//!
//! The coding track draws 5'UTR, CDS and 3'UTR bins end to end with the
//! region boundaries labelled; the ncRNA track, when profiled, is its own
//! axis. Adjacent bins merge (never across a region boundary) so the profile
//! fits the terminal, or to a factor the user picks. A cursor reads out a
//! bar's region, bins, MetaPlotR coordinates and count.

use data_beans::interactive::ui::{
    header, help_line, panel, run_screen, Binning, HistPlot, Scale, Screen, DIM, HIGHLIGHT, PLAIN,
};
use data_beans::qc::pct;
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use super::{GeneFeatureHistogram, CDS, NCRNA, UTR3, UTR5};

/// Width of HistPlot's y gutter.
const GUTTER: u16 = 6;

/// On-screen region names (the TSV's `FEATURE_LABELS` avoid apostrophes).
const REGION_NAMES: [&str; 4] = ["5'UTR", "CDS", "3'UTR", "ncRNA"];

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Track {
    Coding,
    NonCoding,
}

impl Track {
    fn regions(self) -> &'static [usize] {
        match self {
            Track::Coding => &[UTR5, CDS, UTR3],
            Track::NonCoding => &[NCRNA],
        }
    }
}

/// One drawn bar: bins `first..=last` of one region, summed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
struct Bar {
    region: usize,
    first: usize,
    last: usize,
    count: usize,
}

/// Bars of `track`, `merge` bins per bar within each region.
fn bars(hist: &GeneFeatureHistogram, track: Track, merge: usize) -> Vec<Bar> {
    let merge = merge.max(1);
    let mut out = Vec::new();
    for &region in track.regions() {
        let counts = &hist.counts[region];
        for first in (0..counts.len()).step_by(merge) {
            let last = (first + merge).min(counts.len()) - 1;
            out.push(Bar {
                region,
                first,
                last,
                count: counts[first..=last].iter().sum(),
            });
        }
    }
    out
}

/// The smallest merge factor whose bars fit `width` columns; when even one
/// bar per region is too wide, one bar per region.
fn fitting_merge(hist: &GeneFeatureHistogram, track: Track, width: usize) -> usize {
    let lens = track.regions().iter().map(|&r| hist.counts[r].len());
    let (n, longest) = lens.fold((0, 1), |(n, m), l| (n + l, m.max(l)));
    let mut merge = n.div_ceil(width.max(1)).clamp(1, longest);
    while merge < longest && bars(hist, track, merge).len() > width {
        merge += 1;
    }
    merge
}

/// State of the profile view, independent of the terminal so it can be tested.
pub struct MetageneView<'a> {
    title: String,
    hist: &'a GeneFeatureHistogram,
    track: Track,
    /// Bins per bar; `None` fits the terminal.
    merge: Option<usize>,
    /// The merge factor the last frame used.
    shown_merge: usize,
    /// The bin under the cursor, as `(region, bin)`, so it survives merging.
    cursor: (usize, usize),
    y_scale: Scale,
    done: bool,
}

impl<'a> MetageneView<'a> {
    pub fn new(title: &str, hist: &'a GeneFeatureHistogram) -> Self {
        let coding = !Track::Coding
            .regions()
            .iter()
            .all(|&r| hist.counts[r].is_empty());
        let track = if coding {
            Track::Coding
        } else {
            Track::NonCoding
        };
        let mut view = Self {
            title: title.to_string(),
            hist,
            track,
            merge: None,
            shown_merge: 1,
            cursor: (0, 0),
            y_scale: Scale::Linear,
            done: false,
        };
        view.cursor_to_start();
        view
    }

    fn has_non_coding(&self) -> bool {
        !self.hist.counts[NCRNA].is_empty()
    }

    fn bars(&self) -> Vec<Bar> {
        bars(self.hist, self.track, self.shown_merge)
    }

    fn cursor_to_start(&mut self) {
        let first = self
            .track
            .regions()
            .iter()
            .copied()
            .find(|&r| !self.hist.counts[r].is_empty());
        self.cursor = (first.unwrap_or(self.track.regions()[0]), 0);
    }

    /// Index of the bar holding the cursor.
    fn cursor_bar(&self, bars: &[Bar]) -> usize {
        let (region, bin) = self.cursor;
        bars.iter()
            .position(|b| b.region == region && (b.first..=b.last).contains(&bin))
            .unwrap_or(0)
    }

    fn move_cursor(&mut self, delta: isize) {
        let bars = self.bars();
        if bars.is_empty() {
            return;
        }
        let i = self.cursor_bar(&bars) as isize + delta;
        let b = bars[i.clamp(0, bars.len() as isize - 1) as usize];
        self.cursor = (b.region, b.first);
    }

    fn switch_track(&mut self) {
        if self.has_non_coding() && self.track == Track::Coding {
            self.track = Track::NonCoding;
        } else if self.track == Track::NonCoding
            && Track::Coding
                .regions()
                .iter()
                .any(|&r| !self.hist.counts[r].is_empty())
        {
            self.track = Track::Coding;
        }
        self.cursor_to_start();
    }

    fn set_merge(&mut self, delta: isize) {
        let m = (self.shown_merge as isize + delta).max(1) as usize;
        self.merge = Some(m);
        self.shown_merge = m;
    }

    fn stats_line(&self, bars: &[Bar]) -> Line<'static> {
        let dim = |t: String| Span::styled(t, DIM);
        let Some(b) = bars.get(self.cursor_bar(bars)) else {
            return Line::from(dim(" no bins".into()));
        };
        let total: usize = bars.iter().map(|b| b.count).sum();
        let (lo, _) = self.hist.bin_edges(b.region, b.first);
        let (_, hi) = self.hist.bin_edges(b.region, b.last);
        let bins = if b.first == b.last {
            format!("bin {}", b.first)
        } else {
            format!("bins {}-{}", b.first, b.last)
        };
        Line::from(vec![
            Span::styled(format!(" {}", REGION_NAMES[b.region]), HIGHLIGHT),
            dim(format!("  {bins}   coordinate {lo:.3}-{hi:.3}   ")),
            Span::raw(format!("{} sites", b.count)),
            dim(format!(" ({:.2}% of the track)", pct(b.count, total))),
        ])
    }

    fn summary_line(&self) -> Line<'static> {
        let dim = |t: String| Span::styled(t, DIM);
        let mut spans = vec![dim(" sites per region:".into())];
        for &r in self.track.regions() {
            let n: usize = self.hist.counts[r].iter().sum();
            spans.push(dim(format!("  {} ", REGION_NAMES[r])));
            spans.push(Span::raw(n.to_string()));
        }
        if self.track == Track::Coding {
            let m = self.hist.scale.median();
            spans.push(dim(format!(
                "   median nt {:.0}/{:.0}/{:.0}",
                m[UTR5], m[CDS], m[UTR3]
            )));
        }
        Line::from(spans)
    }
}

impl Screen for MetageneView<'_> {
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
            KeyCode::Home => self.cursor_to_start(),
            KeyCode::Tab | KeyCode::BackTab => self.switch_track(),
            KeyCode::Char(']' | '+' | '=') => self.set_merge(1),
            KeyCode::Char('[' | '-') => self.set_merge(-1),
            KeyCode::Char('a') => self.merge = None,
            KeyCode::Char('y') => self.y_scale = self.y_scale.next(),
            KeyCode::Char('q') | KeyCode::Esc | KeyCode::Enter => self.done = true,
            _ => {}
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        let [top, body, stats, summary, footer] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Fill(1),
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Length(1),
        ])
        .areas(frame.area());

        let track_name = match self.track {
            Track::Coding => "coding",
            Track::NonCoding => "ncRNA",
        };
        let block = panel(format!(" {track_name} track "), true);
        let inner = block.inner(body);
        let width = inner.width.saturating_sub(GUTTER) as usize;
        self.shown_merge = self
            .merge
            .unwrap_or_else(|| fitting_merge(self.hist, self.track, width));
        let bars = self.bars();

        let extra = format!(
            "{} bins/bar{} · y {}",
            self.shown_merge,
            if self.merge.is_none() { " (fit)" } else { "" },
            self.y_scale.name()
        );
        frame.render_widget(header("metagene", &self.title, &extra), top);
        frame.render_widget(block, body);

        let counts: Vec<usize> = bars.iter().map(|b| b.count).collect();
        let starts: Vec<Option<&str>> = bars
            .iter()
            .enumerate()
            .map(|(i, b)| {
                (i == 0 || bars[i - 1].region != b.region).then(|| REGION_NAMES[b.region])
            })
            .collect();
        let label = |k: i32| {
            starts
                .get(k as usize)
                .copied()
                .flatten()
                .map(str::to_string)
        };
        HistPlot {
            bins: Binning::with_width(Scale::Linear, 1.0),
            kmin: 0,
            counts: &counts,
            style: &|_| PLAIN,
            subset: None,
            y_scale: self.y_scale,
            pointer: (!bars.is_empty()).then(|| self.cursor_bar(&bars) as i32),
            marks: Vec::new(),
            x_label: Some(&label),
            tick_every: Some(1),
        }
        .render(frame.buffer_mut(), inner);

        frame.render_widget(Paragraph::new(self.stats_line(&bars)), stats);
        frame.render_widget(Paragraph::new(self.summary_line()), summary);
        let mut keys = vec![
            ("←/→", "bar"),
            ("[/]", "merge"),
            ("a", "fit"),
            ("y", "scale"),
        ];
        if self.has_non_coding() {
            keys.push(("Tab", "track"));
        }
        keys.push(("q", "quit"));
        frame.render_widget(help_line(&keys), footer);
    }
}

/// Show the profile full screen until the user quits.
pub fn show_metagene(title: &str, hist: &GeneFeatureHistogram) -> anyhow::Result<()> {
    run_screen(&mut MetageneView::new(title, hist))
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
