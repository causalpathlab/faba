//! Full-screen gene-body browser for `faba pileup --interactive`. Tracks share one axis,
//! re-binned each frame (a bar per column, via [`BinEdges`]); the cursor is a genomic position.
use super::{distinct_positions, fmt_thousands, BinEdges};
use super::{genes_to_draw, SharedModels};
use crate::figure::term::PlotImage;
use crate::figure::{
    self, status_line, Anchor, Bars, Canvas, Controls, Edit, Key, LineInput, INK, MUTED,
};
use crate::site_analysis::miami::genemodel::{draw_gene_model, GeneModel};
use data_beans::interactive::ui::{
    compact, header, help_line, input_line, panel, Binning, HistPlot, MirrorPlot, MirrorSide,
    Scale, Screen, ACCENTED, DIM, HIGHLIGHT, PLAIN,
};
use genomic_data::coordinates::chr_eq;
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

mod track;

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
}

/// Minimum gap between genes in one lane, as a share of the window, so labels stay apart.
const LANE_GAP: f64 = 0.08;

/// Most gene lanes drawn.
const MAX_LANES: usize = 4;

/// Exon height and lane pitch of the gene models in figures.
const GENE_BAND: f64 = 8.0;

const GENE_LANE: f64 = 18.0;

/// How the user left the browser.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Exit {
    Quit,
    /// `g`: back to the gene list.
    Genes,
    /// `/`: a gene or locus the view cannot show itself.
    Search(String),
    /// `-` past the whole view: this chromosome's `lo..=hi`, wider.
    Locus(i64, i64),
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
    /// The genes to draw, once the annotation has arrived.
    genes: Vec<GeneModel>,
    /// The annotation while it is still loading.
    pending_genes: Option<SharedModels>,
    exit: Option<Exit>,
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
            genes: Vec::new(),
            pending_genes: None,
            exit: None,
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
                    let here = drawn.iter().filter(|g| chr_eq(&g.chr, &self.chr));
                    self.extent = here
                        .filter(|g| g.hi > lo && g.lo <= hi)
                        .fold((lo, hi), |(a, b), g| (a.min(g.lo), b.max(g.hi - 1)));
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
        let contrast = self.contrast.map(Row::Contrast);
        let genes = (!self.genes.is_empty()).then_some(Row::Genes);
        let (mirror, first) = if self.mirrored() {
            (Some(Row::Mirror(0, 1)), 2)
        } else {
            (None, 0)
        };
        contrast
            .into_iter()
            .chain(mirror)
            .chain((first..self.tracks.len()).map(Row::Track))
            .chain(genes)
            .collect()
    }

    /// The first two tracks are drawn as one mirrored row (asked for and possible).
    fn mirrored(&self) -> bool {
        self.mirror && self.mirrorable()
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
            .filter(|(t, _)| t.is_reads())
            .map(|(_, s)| s.values.iter().fold(0.0, |m: f64, &v| m.max(v)))
            .fold(0.0, f64::max);
        Bins {
            edges,
            tracks,
            shown,
            shared: (tallest > 0.0).then_some(tallest),
        }
    }

    /// Track `k` as [`Show`] asks, from its binned reads. A track with no
    /// total (or a depth track) draws as it is, but for its sites; what a
    /// mode needs beyond the binned reads is binned only in that mode.
    fn shown(&self, k: usize, edges: &BinEdges, (front, behind): &Binned) -> Shown {
        let t = &self.tracks[k];
        let shown = |values: Vec<f64>, front: Option<Vec<f64>>, accent: bool| Shown {
            values,
            front,
            accent,
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

    fn edges(&self) -> BinEdges {
        BinEdges::new(self.window.0, self.window.1, self.columns)
    }

    /// Bases per bar at the current zoom.
    fn bin_width(&self) -> i64 {
        (self.edges().span().div_ceil(self.columns.max(1) as u64) as i64).max(1)
    }

    /// Genomic range `[start, stop)` of bar `col`.
    fn bar_range(&self, col: usize) -> (i64, i64) {
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

    /// Zoom by `factor` (< 1 in, > 1 out) around the cursor; a bar spans at least one base.
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

    /// Zoom out twice as wide; with all of the view's own span in view,
    /// ask the caller for twice that span around it, flanks and their
    /// genes included, up to [`MAX_SPAN`].
    fn zoom_out(&mut self) {
        let (lo, hi) = self.extent;
        if self.window != self.extent {
            return self.zoom(2.0);
        }
        if hi - lo >= MAX_SPAN {
            self.status = Some("zoomed out as far as it goes".into());
            return;
        }
        let pad = ((hi - lo) / 2).max(self.columns as i64);
        self.exit = Some(Exit::Locus((lo - pad).max(1), hi + pad));
    }

    /// A `/` search: a locus in this view moves there; anything else exits to the caller.
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

    fn cursor_col(&self) -> usize {
        self.edges().col_of(self.cursor)
    }

    /// The visible window as a figure, one panel per track.
    fn figure(&self) -> String {
        let (lo, hi) = self.window;
        let panel_h = 170.0;
        let rows = self.rows();
        let genes_h = 24.0 + GENE_LANE * self.lanes() as f64;
        let heights: Vec<f64> = rows
            .iter()
            .map(|&r| if r == Row::Genes { genes_h } else { panel_h })
            .collect();
        let total: f64 = heights.iter().sum();
        let mut c = Canvas::new(720.0, 50.0 + total);
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
            self.draw_row(row, &bins, &mut c, (8.0, y, 704.0, h), None, title);
            y += h;
        }

        c.finish()
    }

    /// Draw one row (track, contrast, mirror or genes) in the box.
    fn draw_row(
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
                    values: &shown.values,
                    front: shown.front.as_deref(),
                    accent: &|_| shown.accent,
                    colour: &|_| figure::BAR,
                    y_scale: self.y_scale,
                    y_max: bins.shared,
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
                    up: half(a),
                    down: half(b),
                    y_scale: self.y_scale,
                    y_max: bins.shared,
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

impl crate::tui::View for PileupView<'_> {}

impl Screen for PileupView<'_> {
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
            KeyCode::Char('m') => self.mirror = !self.mirror,
            KeyCode::Left | KeyCode::Char('h') => self.move_cursor(-1),
            KeyCode::Right | KeyCode::Char('l') => self.move_cursor(1),
            KeyCode::PageUp => self.move_cursor(-(self.columns as i64) / 2),
            KeyCode::PageDown => self.move_cursor(self.columns as i64 / 2),
            KeyCode::Char('n') => self.jump_site(true),
            KeyCode::Char('p') => self.jump_site(false),
            KeyCode::Char('+' | '=') => self.zoom(0.5),
            KeyCode::Char('-') => self.zoom_out(),
            KeyCode::Char('0') => self.window = self.extent,
            KeyCode::Char('y') => self.y_scale = self.y_scale.next(),
            KeyCode::Char('g') => self.exit = Some(Exit::Genes),
            KeyCode::Char('q') | KeyCode::Esc | KeyCode::Enter => self.exit = Some(Exit::Quit),
            _ => {}
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        self.take_genes();
        let rows = self.rows();
        let n_plots = rows.iter().filter(|&&r| r != Row::Genes).count().max(1) as u32;
        let genes_height = self.lanes() as u16 + 2;
        let mut layout = vec![Constraint::Length(1)];
        layout.extend(rows.iter().map(|&r| match r {
            Row::Genes => Constraint::Length(genes_height),
            _ => Constraint::Ratio(1, n_plots),
        }));
        layout.extend([Constraint::Length(1), Constraint::Length(1)]);
        let areas = Layout::vertical(layout).split(frame.area());
        let (top, readout, footer) = (areas[0], areas[areas.len() - 2], areas[areas.len() - 1]);

        let inner_width = areas[1].width.saturating_sub(2 + GUTTER);
        self.columns = (inner_width as usize).max(1);
        self.clamp_window();
        let (lo, hi) = self.window;
        let extra = format!(
            "{}:{}-{}  {} bp/bar · {} · y {}{}",
            self.chr,
            fmt_thousands(lo),
            fmt_thousands(hi),
            self.bin_width(),
            self.show.label(self.on, self.off),
            self.y_scale.name(),
            self.controls.tag()
        );
        frame.render_widget(header("pileup", &self.title, &extra), top);

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
                    figure::svg(w, h, |c| {
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
                        style: &|_| if shown.accent { ACCENTED } else { PLAIN },
                        subset: shown.front.as_deref(),
                        y_scale: self.y_scale,
                        y_max: bins.shared,
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
                    let (up, down) = figure::split_signed(&values);
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
                        y_max: bins.shared,
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
        let help = if let Some(buf) = self.search.text() {
            input_line(
                "search gene or chr:start-end: ",
                buf,
                &[("Enter", "go"), ("Esc", "back")],
            )
        } else if let Some(msg) = &self.status {
            status_line(msg)
        } else {
            self.controls.footer().unwrap_or_else(|| self.help())
        };
        frame.render_widget(help, footer);
        self.controls.render_save(frame, frame.area());
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
        let col = |pos: i64| chart_x + edges.col_of(pos) as u16;
        for (lane, g) in self.gene_lanes() {
            let y = area.y + lane as u16;
            if y >= area.bottom() {
                break;
            }
            let (c0, c1) = (col(g.lo), col(g.hi - 1).min(right - 1));
            for x in c0..=c1 {
                let sym = match ((x - c0) % 4 == 2, g.forward) {
                    (true, true) => "›",
                    (true, false) => "‹",
                    _ => "─",
                };
                buf[(x, y)].set_symbol(sym).set_style(DIM);
            }
            for &(es, ee) in &g.exons {
                if ee <= self.window.0 || es > self.window.1 {
                    continue;
                }
                for x in col(es)..=col(ee - 1).min(right - 1) {
                    buf[(x, y)].set_symbol("━").set_style(PLAIN);
                }
            }
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

    fn help(&self) -> Line<'static> {
        let mut keys = vec![
            ("←/→", "bar"),
            ("n/p", "next/prev site"),
            ("+/-", "zoom"),
            ("0", "whole"),
            ("/", "search"),
            ("y", "scale"),
            ("g", "genes"),
        ];
        if self.tracks.iter().any(|t| t.is_reads()) {
            keys.push(("c", self.show.next().short()));
        }
        if self.contrast.is_some() {
            keys.push(("d", "difference/fold"));
        }
        if self.mirrorable() {
            keys.push(("m", if self.mirror { "split" } else { "mirror" }));
        }
        self.controls.help_keys(&mut keys);
        keys.push(("q", "quit"));
        help_line(&keys)
    }
}

/// Browse until the user quits, asks for the gene list, or searches outside the view.
/// `status` shows on the footer first.
pub fn show_pileup(view: View, tracks: Vec<Track>, status: Option<String>) -> anyhow::Result<Exit> {
    let mut browser = PileupView::new(view.title, view.chr, tracks, view.extent);
    browser.on = view.on;
    browser.off = view.off;
    browser.keys = view.keys.map(<[_]>::to_vec);
    browser.pending_genes = view.genes;
    browser.take_genes();
    browser.status = status.or_else(|| {
        let loading = browser.pending_genes.is_some();
        loading.then(|| LOADING.to_string())
    });
    browser.controls = Controls::new(&format!("pileup_{}", view.title)).detect();
    crate::tui::run_view(&mut browser)?;
    Ok(browser.exit.unwrap_or(Exit::Quit))
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
