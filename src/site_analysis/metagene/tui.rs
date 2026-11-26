//! Full-screen metagene profile for `faba metagene`.
//!
//! The coding track draws 5'UTR, CDS and 3'UTR bins end to end with the
//! region boundaries labelled; the ncRNA track, when profiled, is its own
//! axis. Adjacent bins merge (never across a region boundary) so the profile
//! fits the terminal, or to a factor the user picks. A cursor reads out a
//! bar's region, bins, MetaPlotR coordinates and count. The bars count the
//! sites' converted reads when the site table gives them (`c` cycles to the
//! unconverted reads, both, and the sites), each region in its colour, as
//! `faba qc` draws its metagene.

use crate::figure::style::{AxisTitles, Sheet};
use crate::tui::ui::{header, help_line, panel, Scale, Screen, DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::figure::term::PlotImage;
use crate::figure::{self, Anchor, Bars, Canvas, Controls, Key, INK, MUTED};
use crate::site_analysis::show::Show;

use super::{
    region_fill, region_key, GeneFeatureHistogram, RegionGlyphs, CDS, NCRNA, REGION_NAMES, UTR3,
    UTR5,
};

use crate::figure::GUTTER;

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

/// One drawn bar: bins `first..=last` of one region, summed; `front`, the
/// converted reads in front of it when both show.
#[derive(Debug, Clone, Copy, PartialEq)]
struct Bar {
    region: usize,
    first: usize,
    last: usize,
    count: f64,
    front: Option<f64>,
}

/// Bars of `track` as `show` draws them, `merge` bins per bar within each
/// region.
fn bars(hist: &GeneFeatureHistogram, track: Track, merge: usize, show: Show) -> Vec<Bar> {
    let merge = merge.max(1);
    let mut out = Vec::new();
    for &region in track.regions() {
        let (counts, front) = hist.shown(show, region);
        for first in (0..counts.len()).step_by(merge) {
            let last = (first + merge).min(counts.len()) - 1;
            out.push(Bar {
                region,
                first,
                last,
                count: counts[first..=last].iter().sum(),
                front: front.map(|f| f[first..=last].iter().sum()),
            });
        }
    }
    out
}

/// Whether any coding region has bins.
fn has_coding(hist: &GeneFeatureHistogram) -> bool {
    Track::Coding
        .regions()
        .iter()
        .any(|&r| !hist.counts[r].is_empty())
}

/// Each bar that starts a region, with the region's name.
fn region_starts(bars: &[Bar]) -> impl Iterator<Item = (usize, &'static str)> + '_ {
    bars.iter()
        .enumerate()
        .filter(|(i, b)| *i == 0 || bars[i - 1].region != b.region)
        .map(|(i, b)| (i, REGION_NAMES[b.region]))
}

/// Each region's name at its middle bar, between the dividers.
fn region_middles(bars: &[Bar]) -> Vec<(usize, String)> {
    let starts: Vec<(usize, &str)> = region_starts(bars).collect();
    let ends = starts.iter().skip(1).map(|&(i, _)| i).chain([bars.len()]);
    starts
        .iter()
        .zip(ends)
        .map(|(&(start, name), end)| (start + (end - start) / 2, name.to_string()))
        .collect()
}

/// A count as shown: whole when it is, else to one decimal.
fn count_text(n: f64) -> String {
    super::count_text(n, 1)
}

/// The smallest merge factor whose bars fit `width` columns; when even one
/// bar per region is too wide, one bar per region.
fn fitting_merge(hist: &GeneFeatureHistogram, track: Track, width: usize) -> usize {
    let lens = track.regions().iter().map(|&r| hist.counts[r].len());
    let (n, longest) = lens.fold((0, 1), |(n, m), l| (n + l, m.max(l)));
    let mut merge = n.div_ceil(width.max(1)).clamp(1, longest);
    let count = |merge: usize| -> usize {
        track
            .regions()
            .iter()
            .map(|&r| hist.counts[r].len().div_ceil(merge))
            .sum()
    };
    while merge < longest && count(merge) > width {
        merge += 1;
    }
    merge
}

/// State of the profile view, independent of the terminal so it can be tested.
pub struct MetageneView<'a> {
    /// Where the last frame drew the view, for a saved figure's shape.
    drawn: Rect,
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
    /// What the bars count: the sites' reads when the table gives them.
    show: Show,
    controls: Controls,
    plot: PlotImage,
    /// The axis titles asked for, in place of the profile's own.
    titles: AxisTitles,
    done: bool,
}

impl<'a> MetageneView<'a> {
    pub fn new(title: &str, hist: &'a GeneFeatureHistogram) -> Self {
        let track = if has_coding(hist) {
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
            show: if hist.has_reads() {
                Show::Converted
            } else {
                Show::Sites
            },
            controls: Controls::new("metagene"),
            drawn: Rect::default(),
            plot: PlotImage::default(),
            titles: AxisTitles::default(),
            done: false,
        };
        view.cursor_to_start();
        view
    }

    fn has_non_coding(&self) -> bool {
        !self.hist.counts[NCRNA].is_empty()
    }

    /// Whether there are a coding and an ncRNA track to switch between:
    /// only then does Tab act.
    fn two_tracks(&self) -> bool {
        self.has_non_coding() && has_coding(self.hist)
    }

    fn bars(&self) -> Vec<Bar> {
        bars(self.hist, self.track, self.shown_merge, self.show)
    }

    /// What the bars count, for the axis and the key.
    fn unit(&self) -> String {
        self.show.label("converted", "unconverted")
    }

    /// Per region, the shown count summed over its bins.
    fn region_total(&self, region: usize) -> f64 {
        self.hist.shown(self.show, region).0.iter().sum()
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

    /// The cursor to the last bar.
    fn cursor_to_end(&mut self) {
        if let Some(b) = self.bars().last() {
            self.cursor = (b.region, b.first);
        }
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
        } else if self.track == Track::NonCoding && has_coding(self.hist) {
            self.track = Track::Coding;
        }
        self.cursor_to_start();
    }

    fn set_merge(&mut self, delta: isize) {
        let m = (self.shown_merge as isize + delta).max(1) as usize;
        self.merge = Some(m);
        self.shown_merge = m;
    }

    /// The shown track as a figure, at the shown merge factor, laid out
    /// for `sheet`'s page: a side not asked for the view's shape on screen,
    /// else a page of its own; its text and margins from the sheet's font.
    /// Alone (`--panel-only`), the bars, their axes and axis titles only.
    fn figure(&self, sheet: Sheet) -> String {
        let (w, h) = sheet.fill(|| {
            self.controls
                .page(self.drawn, LEAST)
                .unwrap_or((680.0, 300.0))
        });
        let style = sheet.style((w, h));
        let mut c = Canvas::styled(w, h, &style);
        let (w, h) = c.page();
        if style.panel_only {
            self.draw_track(&mut c, (0.0, 0.0, w, h), None, String::new());
            return c.finish();
        }
        // The title and the regions' totals within the page's width: the
        // title shortened from the front, the totals wrapped onto further
        // lines, the bars below them taking the rest.
        let room = w - 24.0;
        let title = fit_from_front(&self.title, room, c.label_size(12.0) * 0.6);
        c.bold(16.0, 22.0, &title, 12.0, Anchor::Start, INK);
        let totals: Vec<String> = self
            .track
            .regions()
            .iter()
            .map(|&r| format!("{} {}", REGION_NAMES[r], count_text(self.region_total(r))))
            .collect();
        let head = format!("{} per region:", self.unit());
        let per = c.label_size(9.0) * 0.56;
        let lines = wrap_list(&head, &totals, (room / per) as usize);
        let step = 1.2 * c.label_size(9.0);
        for (i, line) in lines.iter().enumerate() {
            c.text(
                16.0,
                38.0 + step * i as f64,
                line,
                9.0,
                Anchor::Start,
                MUTED,
            );
        }
        let below = step * lines.len().saturating_sub(1) as f64;
        let title = format!("{} bin(s) per bar", self.shown_merge);
        let track = (8.0, 48.0 + below, w - 16.0, h - 56.0 - below);
        self.draw_track(&mut c, track, None, title);
        c.finish()
    }

    /// The shown track's bars in the box `(x, y, w, h)`.
    fn draw_track(
        &self,
        c: &mut Canvas,
        (x, y, w, h): (f64, f64, f64, f64),
        pointer: Option<usize>,
        title: String,
    ) {
        let bars = self.bars();
        let values: Vec<f64> = bars.iter().map(|b| b.count).collect();
        let front: Option<Vec<f64>> = bars.iter().map(|b| b.front).collect();
        let ticks = region_middles(&bars);
        let dividers = region_starts(&bars)
            .map(|(i, _)| i)
            .filter(|&i| i > 0)
            .collect();
        let palette = c.style().palette;
        Bars {
            rotated: false,
            values: &values,
            front: front.as_deref(),
            accent: &|_| false,
            colour: &|i| region_fill(bars[i].region, palette),
            y_scale: self.y_scale,
            y_max: None,
            ticks,
            pointer,
            marks: Vec::new(),
            dividers,
            title,
            x_title: self.titles.x("metagene position (MetaPlotR scale)"),
            y_title: self.titles.y(&self.unit()),
        }
        .draw(c, x, y, w, h);
    }

    fn stats_line(&self, bars: &[Bar]) -> Line<'static> {
        let dim = |t: String| Span::styled(t, DIM);
        let Some(b) = bars.get(self.cursor_bar(bars)) else {
            return Line::from(dim(" no bins".into()));
        };
        let total: f64 = bars.iter().map(|b| b.count).sum();
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
            Span::raw(match b.front {
                Some(f) => format!("{}/{} {}", count_text(f), count_text(b.count), self.unit()),
                None => format!("{} {}", count_text(b.count), self.unit()),
            }),
            dim(format!(
                " ({:.2}% of the track)",
                if total > 0.0 {
                    100.0 * b.count / total
                } else {
                    0.0
                }
            )),
        ])
    }

    fn summary_line(&self) -> Line<'static> {
        let dim = |t: String| Span::styled(t, DIM);
        let mut spans = vec![dim(format!(" {} per region:", self.unit()))];
        for &r in self.track.regions() {
            let n = self.region_total(r);
            spans.push(dim(format!("  {} ", REGION_NAMES[r])));
            spans.push(Span::raw(count_text(n)));
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

impl MetageneView<'_> {
    /// The footer's keys.
    fn footer_pairs(&self) -> Vec<(&'static str, &'static str)> {
        let mut keys = vec![
            ("←/→", "bar"),
            ("Home/End", "first/last"),
            ("+/-", "finer/coarser"),
            ("[/]", "coarser/finer"),
            ("0", "fit"),
            ("y", "scale"),
        ];
        if self.hist.has_reads() {
            keys.push(("c", "reads"));
        }
        if self.two_tracks() {
            keys.push(("Tab", "track"));
        }
        self.controls.help_keys(&mut keys);
        keys.push(("q", "quit"));
        keys.push(("?", "keys"));
        keys
    }
}

impl crate::tui::View for MetageneView<'_> {
    /// Past the save prompt or its outcome, which name none.
    #[cfg(test)]
    fn footer_keys(&self) -> Option<crate::tui::KeyList> {
        (self.controls.footer().is_none()).then(|| crate::tui::keys(&self.footer_pairs()))
    }

    /// The apply key submits the save name being typed.
    fn takes_apply(&self) -> bool {
        self.controls.prompting()
    }

    fn help(&self) -> Option<crate::tui::KeyList> {
        if self.controls.prompting() {
            return None;
        }
        let mut keys = vec![
            ("←/→ h/l", "previous/next bar"),
            ("Home", "first bar"),
            ("End", "last bar"),
            ("+/= ]", "finer: fewer bins a bar"),
            ("- [", "coarser: more bins a bar"),
            ("0", "fit the bars to the width"),
            ("y", "y scale: linear, log10, sqrt"),
        ];
        if self.hist.has_reads() {
            keys.push(("c", "bars count: converted, unconverted, both, sites"));
        }
        if self.two_tracks() {
            keys.push(("Tab/Shift+Tab", "coding or ncRNA track"));
        }
        self.controls.help_keys(&mut keys);
        crate::tui::help::say_save(&mut keys, "save the track as a figure (.pdf)");
        keys.extend([("q", "quit"), ("?", "these keys")]);
        let mut list = crate::tui::keys(&keys);
        list.push((
            "▲".into(),
            "the cursor's bar, read out under the plot".into(),
        ));
        list.push((
            "█".into(),
            "bar colour: its region, keyed under the plot".into(),
        ));
        if self.show == Show::Both {
            list.push((
                "░".into(),
                "lit: converted reads; dim behind: all reads".into(),
            ));
        }
        Some(list)
    }
}

impl Screen for MetageneView<'_> {
    /// For `faba sessions` and the session picker.
    fn status(&self) -> Option<String> {
        Some(format!("metagene of {}", self.title))
    }

    fn done(&self) -> bool {
        self.done
    }

    fn interrupt(&mut self) {
        self.done = true;
    }

    fn handle_key(&mut self, key: KeyEvent) {
        match self.controls.key(key) {
            Key::Pass => self.plot.invalidate(),
            Key::Used => return,
            Key::Save(prefix) => return self.controls.save(&self.figure(Sheet::OWN), &prefix),
        }
        match key.code {
            KeyCode::Left | KeyCode::Char('h') => self.move_cursor(-1),
            KeyCode::Right | KeyCode::Char('l') => self.move_cursor(1),
            KeyCode::Home => self.cursor_to_start(),
            KeyCode::End => self.cursor_to_end(),
            KeyCode::Tab | KeyCode::BackTab if self.two_tracks() => self.switch_track(),
            // As pileup's zoom: `+` finer (fewer bins a bar), `-` coarser.
            KeyCode::Char(']' | '+' | '=') => self.set_merge(-1),
            KeyCode::Char('[' | '-') => self.set_merge(1),
            KeyCode::Char('0') => self.merge = None,
            KeyCode::Char('y') => self.y_scale = self.y_scale.next(),
            KeyCode::Char('c') if self.hist.has_reads() => self.show = self.show.next(),
            KeyCode::Char('q') => self.done = true,
            _ => {}
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        self.controls.refresh();
        self.drawn = frame.area();
        let keys = self.footer_pairs();
        let help = self.controls.footer().unwrap_or_else(|| help_line(&keys));
        let help = crate::tui::wrap_help(help, frame.area().width);
        let [top, body, stats, summary, footer] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Fill(1),
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Length(help.len() as u16),
        ])
        .areas(frame.area());

        let track_name = match self.track {
            Track::Coding => "coding",
            Track::NonCoding => "ncRNA",
        };
        let block = panel(format!(" {track_name} track "), true);
        let [inner, key] =
            Layout::vertical([Constraint::Fill(1), Constraint::Length(1)]).areas(block.inner(body));
        let width = inner.width.saturating_sub(GUTTER) as usize;
        self.shown_merge = self
            .merge
            .unwrap_or_else(|| fitting_merge(self.hist, self.track, width));
        let bars = self.bars();

        let extra = format!(
            "{} bins/bar{} · y {}{}",
            self.shown_merge,
            if self.merge.is_none() { " (fit)" } else { "" },
            self.y_scale.name(),
            self.controls.tag()
        );
        frame.render_widget(header("metagene", &self.title, &extra), top);
        frame.render_widget(block, body);
        let next = self
            .hist
            .has_reads()
            .then(|| ('c', self.show.next().short()));
        let key_line = region_key(&self.unit(), self.track.regions(), next);
        frame.render_widget(Paragraph::new(key_line), key);

        let mut plot = std::mem::take(&mut self.plot);
        let drawn = self.controls.images().is_some_and(|picker| {
            plot.render(frame, inner, picker, |w, h| {
                let pointer = (!bars.is_empty()).then(|| self.cursor_bar(&bars));
                let bbox = (0.0, 0.0, w, h);
                figure::svg(w, h, |c| self.draw_track(c, bbox, pointer, String::new()))
            })
        });
        self.plot = plot;

        if !drawn {
            let values: Vec<f64> = bars.iter().map(|b| b.count).collect();
            let front: Option<Vec<f64>> = bars.iter().map(|b| b.front).collect();
            RegionGlyphs {
                values: &values,
                front: front.as_deref(),
                region: &|i| bars[i].region,
                names: region_middles(&bars),
                pointer: (!bars.is_empty()).then(|| self.cursor_bar(&bars)),
                y_scale: self.y_scale,
            }
            .render(frame.buffer_mut(), inner);
        }

        frame.render_widget(Paragraph::new(self.stats_line(&bars)), stats);
        frame.render_widget(Paragraph::new(self.summary_line()), summary);
        frame.render_widget(Paragraph::new(help), footer);

        self.controls.render_save(frame, frame.area());
    }
}

/// The least page a profile saved from the view grows to, in points, lest a
/// small terminal squeeze its title, totals and axes; a page asked for is
/// taken as it is, its text fit to it.
pub const LEAST: (f64, f64) = (360.0, 200.0);

/// The profile as the view first shows it, drawn without a terminal, laid
/// out for `sheet` (see [`crate::figure::SizeArgs`] and
/// [`crate::figure::style::StyleArgs`]), its axes titled as `titles` asks.
pub fn figure(
    title: &str,
    hist: &GeneFeatureHistogram,
    sheet: Sheet,
    titles: &AxisTitles,
) -> String {
    let mut view = MetageneView::new(title, hist);
    view.titles = titles.clone();
    view.figure(sheet)
}

/// `s` in at most `room` points, `per` points a character: as it is, or
/// shortened from the front behind an ellipsis.
fn fit_from_front(s: &str, room: f64, per: f64) -> String {
    let (n, fits) = (s.chars().count(), (room / per) as usize);
    match n > fits && fits > 1 {
        true => std::iter::once('…')
            .chain(s.chars().skip(n + 1 - fits))
            .collect(),
        false => s.to_string(),
    }
}

/// `head` and then `items`, comma-separated, as lines of at most `fits`
/// characters, broken only between items; an item longer than a line has
/// one of its own.
fn wrap_list(head: &str, items: &[String], fits: usize) -> Vec<String> {
    let mut lines = vec![head.to_string()];
    for (i, item) in items.iter().enumerate() {
        let last = lines.last_mut().expect("a line");
        let comma = if i + 1 < items.len() { "," } else { "" };
        let n = last.chars().count() + 1 + item.chars().count() + comma.len();
        if n <= fits {
            last.push(' ');
            last.push_str(item);
            last.push_str(comma);
        } else {
            lines.push(format!("{item}{comma}"));
        }
    }
    lines
}

/// Show the profile full screen until the user quits.
pub fn show_metagene(title: &str, hist: &GeneFeatureHistogram) -> anyhow::Result<()> {
    let mut view = MetageneView::new(title, hist);
    view.controls = Controls::new("metagene").detect();
    crate::tui::run_view(&mut view)
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
