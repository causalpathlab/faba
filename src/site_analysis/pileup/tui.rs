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
mod render;
mod rows;
mod scroll;
mod stack;
mod track;

use rows::Row;
pub use stack::show_stack;
pub use track::*;

use crate::figure::GUTTER;

/// The footer while the gene models load.
const LOADING: &str = "gene models are loading";

/// Columns between labelled ticks, at least.
const TICK_SPACING: usize = 14;

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
    /// The track whose row to show first under those on top when not all
    /// fit (`[` and `]` scroll them); 0 for the first.
    pub track_from: usize,
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
    /// Tracks flagged [`Track::mirror_with_next`] share a row (`v` splits them).
    mirror: bool,
    /// The track whose row is the first shown under those on top (see
    /// [`scroll`]).
    track_from: usize,
    /// Where the last frame left the rows under those on top.
    fit: scroll::Fit,
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
    /// The browser over `tracks`, each sharing a row with the next as it
    /// is flagged to; the last has none to share one with.
    pub fn new(title: &str, chr: &str, mut tracks: Vec<Track<'a>>, extent: (i64, i64)) -> Self {
        if let Some(last) = tracks.last_mut() {
            last.mirror_with_next = false;
        }
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
            track_from: 0,
            fit: scroll::Fit::default(),
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
        self.set_mirror(s.mirror.unwrap_or(self.mirror));
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
            KeyCode::Char('v') => self.set_mirror(!self.mirror),
            KeyCode::Char('[') => self.scroll_rows(-1),
            KeyCode::Char(']') => self.scroll_rows(1),
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
    browser.track_from = view.track_from;
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
/// leaves it as this one was. How the user left, and the track scrolled
/// to, for a view opened again in its place.
pub fn show_pileup(
    view: View,
    tracks: Vec<Track>,
    status: Option<String>,
    settings: &mut Settings,
) -> anyhow::Result<(Exit, usize)> {
    let mut browser = browser(view, tracks, status);
    browser.apply(settings);
    crate::tui::run_view(&mut browser)?;
    *settings = browser.settings();
    let exit = browser.exit.unwrap_or(Exit::Quit);
    Ok((exit, browser.track_from))
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
