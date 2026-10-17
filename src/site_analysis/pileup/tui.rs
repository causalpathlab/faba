//! Full-screen gene-body browser for `faba pileup`. Tracks share one axis,
//! re-binned each frame (a bar per column, via [`BinEdges`]); the cursor is a genomic position.
//! Moving about is in [`nav`], the saved figure in [`figure`].
use super::{distinct_positions, fmt_thousands, BinEdges};
use super::{genes_to_draw, SharedModels};
use crate::figure::term::PlotImage;
use crate::figure::{status_line, Controls, Edit, Key, LineInput};
use crate::site_analysis::miami::genemodel::GeneModel;
use data_beans::interactive::ui::{
    compact, header, help_line, input_line, panel, Binning, HistPlot, MirrorPlot, MirrorSide,
    Scale, Screen, ACCENTED, DIM, HIGHLIGHT, PLAIN,
};
use genomic_data::coordinates::chr_eq;
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

mod figure;
mod nav;
mod stack;
mod track;

pub use stack::show_stack;
pub use track::*;

use crate::figure::GUTTER;

/// The footer while the gene models load.
const LOADING: &str = "gene models are loading";

/// Columns between labelled ticks, at least.
const TICK_SPACING: usize = 14;

/// A row of the browser.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Row {
    Track(usize),
    Contrast(Contrast),
    /// One track above a zero line and another below it.
    Mirror(usize, usize),
    /// The annotation's genes in view.
    Genes,
}

/// What one browsing session shows besides its tracks.
pub struct View<'a> {
    pub title: &'a str,
    pub chr: &'a str,
    pub extent: (i64, i64),
    /// The converted channel's name, e.g. `methylated`.
    pub on: &'static str,
    /// The other channel's name, e.g. `unmethylated`.
    pub off: &'static str,
    /// The genes whose models to draw; `None` (searched locus) draws every gene in view.
    pub keys: Option<&'a [Box<str>]>,
    /// The annotation's gene models, possibly still loading.
    pub genes: Option<SharedModels>,
    /// Where the cursor starts; the first site when `None`.
    pub cursor: Option<i64>,
    /// Whether `M` switches to the other modality's sites.
    pub switchable: bool,
    /// Tracks drawn in pairs, each a mirrored row: the first above (in the
    /// accent, as converted reads), the second below. Empty: the first two
    /// tracks share a row when they are alike.
    pub pairs: &'a [(usize, usize)],
}

/// Minimum gap between genes in one lane, as a share of the window, so labels stay apart.
const LANE_GAP: f64 = 0.08;

/// Most gene lanes drawn.
const MAX_LANES: usize = 4;

/// Exon height and lane pitch of the gene models in figures.
const GENE_BAND: f64 = 12.0;

const GENE_LANE: f64 = 22.0;

/// How the user left the browser.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Exit {
    Quit,
    /// `g`: back to the gene list.
    Genes,
    /// `/`: a gene or locus the view cannot show itself.
    Search(String),
    /// `-` past the whole view: this chromosome's `lo..=hi`, wider, the
    /// cursor kept `at`.
    Locus {
        lo: i64,
        hi: i64,
        at: i64,
    },
    /// `M`: the other modality's sites here, the cursor kept `at`.
    Modality {
        at: i64,
    },
}

/// The widest span `-` reloads.
const MAX_SPAN: i64 = 10_000_000;

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
    controls: Controls,
    /// The `/` search being typed.
    search: LineInput,
    /// A message for the footer, until the next key.
    status: Option<String>,
    plots: Vec<PlotImage>,
    /// How the first two tracks are compared, when both carry totals.
    contrast: Option<Contrast>,
    /// The converted channel's name, e.g. `methylated`.
    on: &'static str,
    /// The other channel's name, e.g. `unmethylated`.
    off: &'static str,
    /// What the read tracks with a total draw (`c` cycles it).
    show: Show,
    /// The genes whose models to draw; `None` (searched locus) draws every gene in view.
    keys: Option<Vec<Box<str>>>,
    /// The first two tracks share one mirrored row (`m` splits them).
    mirror: bool,
    /// Tracks drawn in pairs, each a mirrored row (see [`View::pairs`]).
    pairs: &'a [(usize, usize)],
    /// The genes to draw, once the annotation has arrived.
    genes: Vec<GeneModel>,
    /// The annotation while it is still loading.
    pending_genes: Option<SharedModels>,
    /// Whether `M` switches to the other modality's sites.
    switchable: bool,
    /// The y-axis maximum `Y` locked, kept for every gene opened after.
    y_lock: Option<f64>,
    /// The header's badge: `pileup`, or in a stack, which gene has the keys.
    badge: &'static str,
    /// Where the last frame drew the view, and the column its bars start
    /// at, for clicks.
    drawn: Rect,
    chart_x: u16,
    /// The screen rows the plots (not the genes) take, for grabbing.
    chart_rows: (u16, u16),
    /// The column a press on the bars came at, until it is let go: a drag
    /// pans, a click moves the cursor.
    grabbed: Option<u16>,
    exit: Option<Exit>,
}

/// What the browser keeps from one gene to the next: how the bars are
/// drawn, not where. `None` leaves a view's own default.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Settings {
    pub y_scale: Option<Scale>,
    pub show: Option<Show>,
    pub mirror: Option<bool>,
    pub measure: Option<Measure>,
    pub y_lock: Option<f64>,
    pub images: Option<bool>,
}

impl<'a> PileupView<'a> {
    pub fn new(title: &str, chr: &str, tracks: Vec<Track<'a>>, extent: (i64, i64)) -> Self {
        let mut sites: Vec<i64> = tracks
            .iter()
            .flat_map(|t| t.front().iter().map(|p| p.0))
            .collect();
        sites.sort_unstable();
        sites.dedup();
        let (lo, hi) = (extent.0.min(extent.1), extent.0.max(extent.1));
        let paired = tracks
            .iter()
            .take(2)
            .filter(|t| t.behind().is_some())
            .count()
            == 2;
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
            controls: Controls::new(&format!("pileup_{title}")),
            search: LineInput::new(128),
            status: None,
            plots: Vec::new(),
            contrast: paired.then_some(Contrast {
                a: 0,
                b: 1,
                measure: Measure::Difference,
            }),
            on: "methylated",
            off: "unmethylated",
            show: Show::default(),
            keys: Some(vec![title.into()]),
            mirror: true,
            pairs: &[],
            genes: Vec::new(),
            pending_genes: None,
            switchable: false,
            y_lock: None,
            badge: "pileup",
            drawn: Rect::default(),
            chart_x: 0,
            chart_rows: (0, 0),
            grabbed: None,
            exit: None,
        }
    }

    /// Draw as `s` keeps from the gene before (or another in a stack); a
    /// plot drawn as an image is redrawn when that changes it.
    pub fn apply(&mut self, s: &Settings) {
        let before = self.settings();
        self.y_scale = s.y_scale.unwrap_or(self.y_scale);
        self.show = s.show.unwrap_or(self.show);
        self.mirror = s.mirror.unwrap_or(self.mirror);
        if let (Some(c), Some(measure)) = (&mut self.contrast, s.measure) {
            c.measure = measure;
        }
        self.y_lock = s.y_lock;
        if let Some(on) = s.images {
            self.controls.set_images(on);
        }
        if self.settings() != before {
            self.plots.iter_mut().for_each(PlotImage::invalidate);
        }
    }

    /// A click at screen column `x`: the cursor to the bar there, when it
    /// is over the bars; whether it moved.
    pub(super) fn click_at(&mut self, x: u16) -> bool {
        let Some(col) = self.bar_at(x) else {
            return false;
        };
        let (start, stop) = self.bar_range(col);
        self.cursor = (start + (stop - start) / 2).clamp(self.extent.0, self.extent.1);
        self.plots.iter_mut().for_each(PlotImage::invalidate);
        true
    }

    /// How this view draws, to keep for the next gene.
    pub fn settings(&self) -> Settings {
        Settings {
            y_scale: Some(self.y_scale),
            show: Some(self.show),
            mirror: Some(self.mirror),
            measure: self.contrast.map(|c| c.measure),
            y_lock: self.y_lock,
            images: Some(self.controls.images_on()),
        }
    }

    /// Take this chromosome's genes within the extent once the annotation has arrived.
    fn take_genes(&mut self) {
        let Some(read) = self.pending_genes.as_ref().and_then(|m| m.get()) else {
            return;
        };
        match read {
            Ok(models) => {
                if self.status.as_deref() == Some(LOADING) {
                    self.status = None;
                }
                let drawn = genes_to_draw(models, self.keys.as_deref());
                // The view's own genes widen it to their TSS and TES, so
                // zooming out reaches the whole gene; a searched locus
                // (no keys) keeps its bounds.
                if self.keys.is_some() {
                    let (lo, hi) = self.extent;
                    // Not zoomed yet: the window follows to the whole gene.
                    let whole = self.window == self.extent;
                    let here = drawn.iter().filter(|g| chr_eq(&g.chr, &self.chr));
                    self.extent = here
                        .filter(|g| g.hi > lo && g.lo <= hi)
                        .fold((lo, hi), |(a, b), g| (a.min(g.lo), b.max(g.hi - 1)));
                    if whole {
                        self.window = self.extent;
                    }
                }
                let (lo, hi) = self.extent;
                self.genes = drawn
                    .into_iter()
                    .filter(|g| chr_eq(&g.chr, &self.chr) && g.hi > lo && g.lo <= hi)
                    .cloned()
                    .collect();
            }
            Err(e) => self.status = Some(format!("gene models: {e}")),
        }
        self.pending_genes = None;
        self.plots.iter_mut().for_each(PlotImage::invalidate);
    }

    /// The rows top to bottom: contrast, tracks, then genes.
    fn rows(&self) -> Vec<Row> {
        let pairs = self.mirror_pairs();
        let paired = |k: usize| pairs.iter().any(|&(a, b)| k == a || k == b);
        let mut rows: Vec<Row> = self.contrast.map(Row::Contrast).into_iter().collect();
        // A pair's row where its first track is; the rest as they are.
        for k in 0..self.tracks.len() {
            match pairs.iter().find(|p| p.0 == k) {
                Some(&(a, b)) => rows.push(Row::Mirror(a, b)),
                None if !paired(k) => rows.push(Row::Track(k)),
                None => {}
            }
        }
        rows.extend((!self.genes.is_empty()).then_some(Row::Genes));
        rows
    }

    /// The tracks drawn in pairs: the view's, else the first two when they
    /// are alike; none when split (`m`).
    fn mirror_pairs(&self) -> Vec<(usize, usize)> {
        if !self.mirror {
            Vec::new()
        } else if !self.pairs.is_empty() {
            self.pairs.to_vec()
        } else if self.mirrorable() {
            vec![(0, 1)]
        } else {
            Vec::new()
        }
    }

    /// The top of `row`'s scale: the shared one, but for a track on its own
    /// scale, locked or not.
    fn row_max(&self, row: Row, bins: &Bins) -> Option<f64> {
        match row {
            Row::Mirror(k, _) | Row::Track(k) if self.tracks[k].own_scale => None,
            _ => bins.shared,
        }
    }

    /// The first two tracks are read tracks of the same measure.
    fn mirrorable(&self) -> bool {
        match &self.tracks[..] {
            [a, b, ..] => a.is_reads() && b.is_reads() && a.name == b.name,
            _ => false,
        }
    }

    fn row_title(&self, row: Row) -> String {
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
    fn gene_lanes(&self) -> Vec<(usize, &GeneModel)> {
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
    fn lanes(&self) -> usize {
        let lanes = self.gene_lanes().iter().map(|(l, _)| l + 1).max();
        lanes.unwrap_or(1)
    }

    /// Every track binned over the window, and as drawn.
    fn bins(&self) -> Bins {
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
    fn shown(&self, k: usize, edges: &BinEdges, (front, behind): &Binned) -> Shown {
        let t = &self.tracks[k];
        // The upper half of a given pair is converted reads: the accent.
        let top = self.pairs.iter().any(|&(a, _)| a == k);
        let shown = |values: Vec<f64>, front: Option<Vec<f64>>, accent: bool| Shown {
            values,
            front,
            accent: accent || top,
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
    fn y_title(&self, k: usize) -> String {
        let t = &self.tracks[k];
        let sites = self.show == Show::Sites && t.is_reads();
        if sites || t.behind().is_some() {
            self.show.label(self.on, self.off)
        } else {
            t.name.to_string()
        }
    }

    /// The contrast per bar, from raw sums (rebinned when the tracks show logs).
    fn contrast_values(&self, c: Contrast, bins: &Bins) -> Vec<Option<f64>> {
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

    fn plot_count(rows: &[Row]) -> usize {
        rows.iter().filter(|&&r| r != Row::Genes).count().max(1)
    }

    /// A press, release or wheel turn on the view: a press on the bars
    /// grabs them, and letting go elsewhere pans as far (in the same place,
    /// it moves the cursor there); the wheel, either way, pans a tenth of
    /// the view a notch. Whether it changed.
    pub(super) fn pointer(&mut self, m: ratatui::crossterm::event::MouseEvent) -> bool {
        use ratatui::crossterm::event::MouseEventKind;
        let notch = (self.columns as i64 / 10).max(1);
        match m.kind {
            // Any button: a click with any of them moves the cursor.
            MouseEventKind::Down(_) => {
                self.grabbed = self.over_bars(m.column, m.row).then_some(m.column);
                false
            }
            MouseEventKind::Up(_) => match self.grabbed.take() {
                // Dragged: the bars follow the pointer, as when grabbed.
                Some(from) if from.abs_diff(m.column) >= 2 => {
                    self.pan(i64::from(from) - i64::from(m.column))
                }
                Some(from) => self.click_at(from),
                None => false,
            },
            MouseEventKind::ScrollUp | MouseEventKind::ScrollLeft => self.pan(-notch),
            MouseEventKind::ScrollDown | MouseEventKind::ScrollRight => self.pan(notch),
            _ => false,
        }
    }

    fn readout(&self, bins: &Bins) -> Line<'static> {
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
        if let Some(c) = self.contrast {
            let (a, b) = (self.tracks[c.a].label, self.tracks[c.b].label);
            let value = self.contrast_values(c, bins).get(col).copied().flatten();
            spans.push(dim(format!("   {a} vs {b} ")));
            spans.push(Span::styled(
                value.map_or("-".into(), |v| c.measure.label(v)),
                HIGHLIGHT,
            ));
        }
        for (t, shown) in self.tracks.iter().zip(&bins.shown) {
            spans.push(dim(format!("   {} ", t.label)));
            let v = shown.values.get(col).copied().unwrap_or(0.0);
            spans.push(Span::raw(
                match shown.front.as_ref().and_then(|f| f.get(col)) {
                    Some(f) => format!("{f:.0}/{v:.0}"),
                    None => compact(v),
                },
            ));
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

impl crate::tui::View for PileupView<'_> {
    /// The apply key submits the search or the save name being typed.
    fn takes_apply(&self) -> bool {
        self.search.active() || self.controls.prompting()
    }
}

impl Screen for PileupView<'_> {
    fn takes_mouse(&self) -> bool {
        true
    }

    /// A click on the bars moves the cursor there; a drag or the wheel
    /// pans (see [`PileupView::pointer`]).
    fn mouse(&mut self, m: ratatui::crossterm::event::MouseEvent) -> bool {
        self.pointer(m)
    }

    fn done(&self) -> bool {
        self.exit.is_some()
    }

    fn interrupt(&mut self) {
        self.exit = Some(Exit::Quit);
    }

    /// Redraw once the gene models arrive in the background.
    fn tick(&mut self) -> bool {
        let waiting = self.pending_genes.is_some();
        self.take_genes();
        waiting && self.pending_genes.is_none()
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.status = None;
        if self.search.active() {
            if let Edit::Submitted(query) = self.search.handle(key) {
                self.submit(&query);
                self.plots.iter_mut().for_each(PlotImage::invalidate);
            }
            return;
        }
        match self.controls.key(key) {
            Key::Pass => self.plots.iter_mut().for_each(PlotImage::invalidate),
            Key::Used => return,
            Key::Save(prefix) => return self.controls.save(&self.figure(), &prefix),
        }
        match key.code {
            KeyCode::Char('/') => self.search.open(""),
            KeyCode::Char('c') => self.show = self.show.next(),
            KeyCode::Char('d') => {
                self.contrast = self.contrast.map(|c| Contrast {
                    measure: c.measure.next(),
                    ..c
                })
            }
            KeyCode::Char('v') => self.mirror = !self.mirror,
            KeyCode::Char('m' | 'M') if self.switchable => {
                self.exit = Some(Exit::Modality { at: self.cursor })
            }
            KeyCode::Left | KeyCode::Char('h') => self.move_cursor(-1),
            KeyCode::Right | KeyCode::Char('l') => self.move_cursor(1),
            KeyCode::PageUp => self.move_cursor(-(self.columns as i64) / 2),
            KeyCode::PageDown => self.move_cursor(self.columns as i64 / 2),
            KeyCode::Char('n') => self.jump_site(true),
            KeyCode::Char('N') => self.jump_site(false),
            KeyCode::Char('+' | '=') => self.zoom(0.5),
            KeyCode::Char('-') => self.zoom_out(),
            KeyCode::Char('0') => self.window = self.extent,
            KeyCode::Char('y') => self.y_scale = self.y_scale.next(),
            // Keep this y-axis maximum for every gene opened after, so they
            // compare; again, back to each its own.
            KeyCode::Char('Y') => {
                self.y_lock = match self.y_lock {
                    Some(_) => None,
                    None => self.bins().shared,
                }
            }
            // Esc goes back one level: the gene list it came from.
            KeyCode::Char('g') | KeyCode::Esc => self.exit = Some(Exit::Genes),
            KeyCode::Char('q') => self.exit = Some(Exit::Quit),
            _ => {}
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        self.render_in(frame, frame.area(), true);
    }
}

impl PileupView<'_> {
    /// Draw in `area`: the header, the rows and the readout, and with
    /// `footer` the keys and the save pop-up (a stack draws those once).
    pub(super) fn render_in(&mut self, frame: &mut Frame, area: Rect, footer: bool) {
        self.take_genes();
        let rows = self.rows();
        let plot_rows = Self::plot_count(&rows);
        let n_plots = plot_rows as u32;
        let genes_height = self.lanes() as u16 + 2;
        let mut layout = vec![Constraint::Length(1)];
        layout.extend(rows.iter().map(|&r| match r {
            Row::Genes => Constraint::Length(genes_height),
            _ => Constraint::Ratio(1, n_plots),
        }));
        // The footer takes as many lines as its keys need at this width.
        let help = match footer {
            true => crate::tui::wrap_help(self.footer_line(), area.width),
            false => Vec::new(),
        };
        layout.extend([Constraint::Length(1), Constraint::Length(help.len() as u16)]);
        let areas = Layout::vertical(layout).split(area);
        let (top, readout, footer_area) =
            (areas[0], areas[areas.len() - 2], areas[areas.len() - 1]);

        let inner_width = areas[1].width.saturating_sub(2 + GUTTER);
        self.columns = (inner_width as usize).max(1);
        // Inside the panel's border, past the y axis.
        (self.drawn, self.chart_x) = (area, areas[1].x + 1 + GUTTER);
        self.chart_rows = (areas[1].y, areas[plot_rows].bottom());
        self.clamp_window();
        let (lo, hi) = self.window;
        let extra = format!(
            "{}:{}-{}  {} bp/bar · {} · y {}{}{}",
            self.chr,
            fmt_thousands(lo),
            fmt_thousands(hi),
            self.bin_width(),
            self.show.label(self.on, self.off),
            self.y_scale.name(),
            match self.y_lock {
                Some(max) => format!(" ≤ {max:.0} locked"),
                None => String::new(),
            },
            self.controls.tag()
        );
        frame.render_widget(header(self.badge, &self.title, &extra), top);

        let every = TICK_SPACING.max(self.columns / 5);
        let ticks: Vec<Option<String>> = (0..self.columns)
            .map(|k| {
                (k.is_multiple_of(every) && k + 10 < self.columns)
                    .then(|| fmt_thousands(self.bar_range(k).0))
            })
            .collect();
        let column_label = |k: usize| ticks.get(k).cloned().flatten();
        let label = |k: i32| column_label(k as usize);
        let bins = self.bins();
        let pointer = self.cursor_col();
        let mut plots = std::mem::take(&mut self.plots);
        plots.resize_with(rows.len(), PlotImage::default);
        for (k, (&row, &area)) in rows.iter().zip(&areas[1..areas.len() - 2]).enumerate() {
            let block = panel(format!(" {} ", self.row_title(row)), true);
            let inner = block.inner(area);
            frame.render_widget(block, area);
            let drawn = self.controls.images().is_some_and(|picker| {
                plots[k].render(frame, inner, picker, |w, h| {
                    let bbox = (0.0, 0.0, w, h);
                    crate::figure::svg(w, h, |c| {
                        self.draw_row(row, &bins, c, bbox, Some(pointer), String::new())
                    })
                })
            });
            if drawn {
                continue;
            }
            let buf = frame.buffer_mut();
            match row {
                Row::Track(i) => {
                    let marks = self.tracks[i].sites_in(lo, hi).into_iter();
                    let shown = &bins.shown[i];
                    HistPlot {
                        bins: Binning::with_width(Scale::Linear, 1.0),
                        kmin: 0,
                        counts: &shown.values,
                        style: &|_| {
                            if shown.accent {
                                ACCENTED
                            } else {
                                PLAIN
                            }
                        },
                        subset: shown.front.as_deref(),
                        y_scale: self.y_scale,
                        y_max: self.row_max(row, &bins),
                        pointer: Some(pointer as i32),
                        marks: marks
                            .map(|s| (bins.edges.col_of(s) as i32, "+", DIM))
                            .collect(),
                        x_label: Some(&label),
                        tick_every: Some(1),
                    }
                    .render(buf, inner);
                }
                Row::Contrast(contrast) => {
                    let values = self.contrast_values(contrast, &bins);
                    let (up, down) = crate::figure::split_signed(&values);
                    let max = up.iter().chain(&down).fold(0.0, |m: f64, &v| m.max(v));
                    let side = |counts, style| MirrorSide {
                        counts,
                        subset: None,
                        style,
                        name: "",
                    };
                    let measure = contrast.measure;
                    MirrorPlot {
                        up: side(&up, ACCENTED),
                        down: side(&down, PLAIN),
                        y_scale: Scale::Linear,
                        y_max: Some(max),
                        y_labels: Some([max, 0.0, -max].map(|v| measure.label(v))),
                        pointer: Some(pointer),
                        x_label: Some(&column_label),
                    }
                    .render(buf, inner);
                }
                Row::Mirror(a, b) => {
                    let side = |k: usize| {
                        let shown = &bins.shown[k];
                        MirrorSide {
                            counts: &shown.values,
                            subset: shown.front.as_deref(),
                            style: if shown.accent { ACCENTED } else { PLAIN },
                            name: self.tracks[k].label,
                        }
                    };
                    MirrorPlot {
                        up: side(a),
                        down: side(b),
                        y_scale: self.y_scale,
                        y_max: self.row_max(row, &bins),
                        y_labels: None,
                        pointer: Some(pointer),
                        x_label: Some(&column_label),
                    }
                    .render(buf, inner);
                }
                Row::Genes => self.render_genes(buf, inner),
            }
        }
        self.plots = plots;

        frame.render_widget(Paragraph::new(self.readout(&bins)), readout);
        if footer {
            frame.render_widget(Paragraph::new(help), footer_area);
            self.controls.render_save(frame, frame.area());
        }
    }
}

impl PileupView<'_> {
    /// The genes row as text, a lane per line: heavy exons on a thin arrowed intron line.
    fn render_genes(&self, buf: &mut ratatui::buffer::Buffer, area: ratatui::layout::Rect) {
        let chart_x = area.x + GUTTER;
        if area.width <= GUTTER || area.height == 0 {
            return;
        }
        let right = area.right();
        let edges = self.edges();
        let col = |pos: i64| (chart_x + edges.col_of(pos) as u16).min(right - 1);
        for (lane, g) in self.gene_lanes() {
            let y = area.y + lane as u16;
            if y >= area.bottom() {
                break;
            }
            let (c0, c1) = draw_gene_line(buf, y, col, g, self.window, PLAIN);
            // The symbol after the gene if it fits, else before it.
            let label = g.symbol.as_ref();
            let len = label.chars().count() as u16;
            let lx = if c1 + 2 + len <= right {
                c1 + 2
            } else {
                c0.saturating_sub(len + 1).max(area.x)
            };
            buf.set_string(lx, y, label, DIM);
        }
    }

    /// The footer: the search being typed, a status, the save prompt, or
    /// the keys.
    fn footer_line(&self) -> Line<'static> {
        if let Some(buf) = self.search.text() {
            input_line(
                "search gene or chr:start-end: ",
                buf,
                &[("Enter", "go"), ("Esc", "back")],
            )
        } else if let Some(msg) = &self.status {
            status_line(msg)
        } else {
            self.controls.footer().unwrap_or_else(|| self.help())
        }
    }

    fn help(&self) -> Line<'static> {
        let mut keys = vec![
            ("←/→", "bar"),
            ("PgUp/PgDn", "half a view"),
            ("n/N", "next/prev site"),
            ("+/-", "zoom"),
            ("0", "whole"),
            ("/", "search"),
            ("y", "scale"),
            (
                "Y",
                if self.y_lock.is_some() {
                    "unlock y"
                } else {
                    "lock y"
                },
            ),
            ("g/Esc", "gene list"),
        ];
        if self.tracks.iter().any(|t| t.is_reads()) {
            keys.push(("c", self.show.next().short()));
        }
        if self.contrast.is_some() {
            keys.push(("d", "difference/fold"));
        }
        if self.mirrorable() {
            keys.push(("v", if self.mirror { "split" } else { "mirror" }));
        }
        if self.switchable {
            keys.push(("m", "m6A/A-to-I"));
        }
        self.controls.help_keys(&mut keys);
        keys.push(("q", "quit"));
        help_line(&keys)
    }
}

/// Gene model `g` on row `y`, columns from `col`: a thin intron line with
/// a strand arrow every 4th column, then full blocks for the exons that
/// reach into `shown` (a block reads thick against the line), in
/// `exon_style`. The line's first and last columns.
pub(crate) fn draw_gene_line(
    buf: &mut ratatui::buffer::Buffer,
    y: u16,
    col: impl Fn(i64) -> u16,
    g: &GeneModel,
    shown: (i64, i64),
    exon_style: ratatui::style::Style,
) -> (u16, u16) {
    let (c0, c1) = (col(g.lo), col(g.hi - 1));
    for x in c0..=c1 {
        let sym = match ((x - c0) % 4 == 2, g.forward) {
            (true, true) => "›",
            (true, false) => "‹",
            _ => "─",
        };
        buf.set_string(x, y, sym, DIM);
    }
    for &(es, ee) in &g.exons {
        if ee <= shown.0 || es > shown.1 {
            continue;
        }
        for x in col(es)..=col(ee - 1) {
            buf.set_string(x, y, "█", exon_style);
        }
    }
    (c0, c1)
}

/// The browser for `view`'s `tracks`, `status` on its footer first.
fn browser<'a>(view: View<'a>, tracks: Vec<Track<'a>>, status: Option<String>) -> PileupView<'a> {
    let mut browser = PileupView::new(view.title, view.chr, tracks, view.extent);
    browser.on = view.on;
    browser.off = view.off;
    browser.keys = view.keys.map(<[_]>::to_vec);
    browser.pending_genes = view.genes;
    browser.take_genes();
    browser.switchable = view.switchable;
    browser.pairs = view.pairs;
    if let Some(at) = view.cursor {
        browser.cursor = at.clamp(browser.extent.0, browser.extent.1);
    }
    browser.status = status.or_else(|| {
        let loading = browser.pending_genes.is_some();
        loading.then(|| LOADING.to_string())
    });
    // Text plots on screen; only the saved figure is drawn as a picture.
    browser.controls = Controls::new(&format!("pileup_{}", view.title));
    browser
}

/// Browse until the user quits, asks for the gene list, or searches outside the view.
/// `status` shows on the footer first.
/// `settings` holds how the gene before was drawn: this one starts so, and
/// leaves it as this one was.
pub fn show_pileup(
    view: View,
    tracks: Vec<Track>,
    status: Option<String>,
    settings: &mut Settings,
) -> anyhow::Result<Exit> {
    let mut browser = browser(view, tracks, status);
    browser.apply(settings);
    crate::tui::run_view(&mut browser)?;
    *settings = browser.settings();
    Ok(browser.exit.unwrap_or(Exit::Quit))
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
