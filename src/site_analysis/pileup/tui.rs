//! Full-screen gene-body browser for `faba pileup --interactive`.
//!
//! The matrix track, and the site track when a site table was given, share
//! one genomic axis. Every frame re-bins the raw positions of the visible
//! window, one bar per terminal column, with the same rule as the printed
//! pileup ([`BinEdges`]), so zooming in resolves single sites. The cursor is a
//! genomic coordinate, so it stays put through zooms and resizes.

use data_beans::interactive::ui::{
    compact, header, help_line, input_line, panel, Binning, HistPlot, Scale, Screen, ACCENTED, DIM,
    HIGHLIGHT, PLAIN,
};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;

use crate::figure::term::PlotImage;
use crate::figure::{
    self, status_line, Anchor, Bars, Canvas, Controls, Edit, Key, LineInput, INK, MUTED,
};

use genomic_data::coordinates::chr_eq;

use super::{distinct_positions, fmt_thousands, BinEdges};
use super::{genes_to_draw, SharedModels};
use crate::site_analysis::miami::genemodel::{gene_model_svg, GeneModel};

/// Width of HistPlot's y gutter.
const GUTTER: u16 = 6;

/// Columns between labelled ticks, at least.
const TICK_SPACING: usize = 14;

/// One track: sorted `(position, value)` pairs, optionally drawn in front
/// of a total per position (converted reads over all reads).
pub struct Track<'a> {
    pub label: &'a str,
    /// What the values are, for the panel title.
    pub name: &'a str,
    pub front: &'a [(i64, f64)],
    pub behind: Option<&'a [(i64, f64)]>,
    /// Genomic bins `(start, end, value)` instead of positions: each column
    /// shows the bin covering it (read depth).
    pub ranges: Option<&'a [(i64, i64, f64)]>,
    /// Bins become `log10(1 + sum)`, as the printed pileup does.
    pub log: bool,
}

impl<'a> Track<'a> {
    /// A track with no total.
    #[cfg(test)]
    pub fn single(label: &'a str, name: &'a str, front: &'a [(i64, f64)], log: bool) -> Self {
        Track {
            label,
            name,
            front,
            behind: None,
            ranges: None,
            log,
        }
    }

    /// A read-depth track over genomic bins.
    pub fn depth(label: &'a str, ranges: &'a [(i64, i64, f64)]) -> Self {
        Track {
            label,
            name: "reads per depth bin",
            front: &[],
            behind: None,
            ranges: Some(ranges),
            log: false,
        }
    }

    /// Binned over `edges`: front, and the total when there is one.
    fn bin(&self, edges: &BinEdges) -> (Vec<f64>, Option<Vec<f64>>) {
        if let Some(ranges) = self.ranges {
            return (ranges_per_column(ranges, edges), None);
        }
        let front = edges.bin(self.front, self.log);
        (front, self.behind.map(|b| edges.bin(b, self.log)))
    }

    /// Distinct positions inside `lo..=hi`.
    fn sites_in(&self, lo: i64, hi: i64) -> Vec<i64> {
        let a = self.front.partition_point(|p| p.0 < lo);
        let b = self.front.partition_point(|p| p.0 <= hi);
        distinct_positions(&self.front[a..b])
    }
}

/// Per column, the value of the genomic bin covering the column's middle.
fn ranges_per_column(ranges: &[(i64, i64, f64)], edges: &BinEdges) -> Vec<f64> {
    let (n, span) = (edges.num_bins as i64, edges.span() as i64);
    (0..n)
        .map(|k| {
            let mid = edges.min_pos + (2 * k + 1) * span / (2 * n);
            let i = ranges.partition_point(|r| r.1 <= mid);
            ranges.get(i).filter(|r| r.0 <= mid).map_or(0.0, |r| r.2)
        })
        .collect()
}

/// How the contrast row compares two tracks' converted fractions per bar.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Measure {
    /// Fraction A minus fraction B, in percentage points.
    Difference,
    /// log2 of fraction A over fraction B, half a read of pseudocount.
    Log2Fold,
}

impl Measure {
    fn name(self, on: &str) -> String {
        match self {
            Measure::Difference => format!("{on} fraction difference (pp)"),
            Measure::Log2Fold => format!("log2 fold of {on} fraction"),
        }
    }

    fn next(self) -> Self {
        match self {
            Measure::Difference => Measure::Log2Fold,
            Measure::Log2Fold => Measure::Difference,
        }
    }

    /// The measure for converted `(ma, na)` of A and `(mb, nb)` of B;
    /// `None` when either has no reads.
    pub fn of(self, (ma, na): (f64, f64), (mb, nb): (f64, f64)) -> Option<f64> {
        if na <= 0.0 || nb <= 0.0 {
            return None;
        }
        Some(match self {
            Measure::Difference => 100.0 * (ma / na - mb / nb),
            Measure::Log2Fold => {
                ((ma + 0.5) / (na + 1.0)).log2() - ((mb + 0.5) / (nb + 1.0)).log2()
            }
        })
    }

    fn label(self, v: f64) -> String {
        match self {
            Measure::Difference => format!("{v:+.1}"),
            Measure::Log2Fold => format!("{v:+.2}"),
        }
    }
}

/// A row of the browser: a track, or the contrast of the first two.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum Row {
    Track(usize),
    Contrast,
    /// The annotation's genes in view.
    Genes,
}

/// Signed bars around a middle zero line, as text: positive up in the
/// accent, negative down; `None` draws nothing.
struct TextDiverging<'a> {
    values: &'a [Option<f64>],
    pointer: usize,
    label: &'a dyn Fn(f64) -> String,
    /// Tick label per column.
    ticks: &'a [Option<String>],
}

impl TextDiverging<'_> {
    fn render(&self, buf: &mut ratatui::buffer::Buffer, area: ratatui::layout::Rect) {
        let [plot, axis] =
            Layout::vertical([Constraint::Min(1), Constraint::Length(1)]).areas(area);
        let [gutter, chart] =
            Layout::horizontal([Constraint::Length(GUTTER), Constraint::Min(1)]).areas(plot);
        if chart.width == 0 || chart.height < 3 {
            return;
        }
        let half = (chart.height - 1) / 2;
        let zero = chart.top() + half;
        let max = self
            .values
            .iter()
            .flatten()
            .fold(0.0f64, |m, v| m.max(v.abs()));
        for x in chart.left()..chart.right() {
            buf[(x, zero)].set_symbol("─").set_style(DIM);
        }
        for (i, v) in self.values.iter().enumerate() {
            let x = chart.x + i as u16;
            if x >= chart.right() {
                break;
            }
            let Some(v) = *v else { continue };
            if max <= 0.0 {
                continue;
            }
            let cells = ((v.abs() / max) * half as f64).round().max(1.0) as u16;
            for k in 0..cells.min(half) {
                let (y, style) = if v >= 0.0 {
                    (zero - 1 - k, ACCENTED)
                } else {
                    (zero + 1 + k, PLAIN)
                };
                buf[(x, y)].set_symbol("█").set_style(style);
            }
        }
        let gx = gutter.right() - 1;
        for y in gutter.top()..gutter.bottom() {
            buf[(gx, y)].set_symbol("│").set_style(DIM);
        }
        for (y, v) in [
            (chart.top(), max),
            (zero, 0.0),
            (chart.top() + 2 * half, -max),
        ] {
            let s = (self.label)(v);
            buf.set_string(
                gx.saturating_sub(s.len() as u16 + 1).max(gutter.x),
                y,
                &s,
                DIM,
            );
        }
        for x in axis.left()..axis.right() {
            buf[(x, axis.y)].set_symbol(" ");
        }
        let mut next_free = axis.x;
        for (i, t) in self.ticks.iter().enumerate() {
            let x = chart.x + i as u16;
            if let Some(t) = t {
                if x >= next_free && x + (t.len() as u16) <= axis.right() {
                    buf.set_string(x, axis.y, t, DIM);
                    next_free = x + t.len() as u16 + 1;
                }
            }
        }
        let px = chart.x + self.pointer as u16;
        if px < chart.right() {
            buf[(px, axis.y)].set_symbol("▲").set_style(HIGHLIGHT);
        }
    }
}

/// What one browsing session shows besides its tracks.
pub struct View<'a> {
    pub title: &'a str,
    pub chr: &'a str,
    pub extent: (i64, i64),
    /// The two channels' names, e.g. `("methylated", "unmethylated")`.
    pub channels: (&'static str, &'static str),
    /// The annotation's gene models, possibly still loading.
    pub genes: Option<SharedModels>,
}

/// Genes sharing a lane leave at least this share of the window between
/// them, so their labels stay apart.
const LANE_GAP: f64 = 0.08;

/// Most gene lanes drawn.
const MAX_LANES: usize = 4;

/// Exon height and lane pitch of the gene models in figures.
const GENE_BAND: f32 = 8.0;
const GENE_LANE: f64 = 18.0;

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
    controls: Controls,
    /// The `/` search being typed.
    search: LineInput,
    /// A message for the footer, until the next key.
    status: Option<String>,
    plots: Vec<PlotImage>,
    /// How the first two tracks are compared, when both carry totals.
    contrast: Option<Measure>,
    /// The two channels' names, e.g. `("methylated", "unmethylated")`.
    channels: (&'static str, &'static str),
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
            .flat_map(|t| t.front.iter().map(|p| p.0))
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
            controls: Controls::new(&format!("pileup_{title}")),
            search: LineInput::new(128),
            status: None,
            plots: Vec::new(),
            contrast: None,
            channels: ("methylated", "unmethylated"),
            genes: Vec::new(),
            pending_genes: None,
            exit: None,
        }
        .with_contrast()
    }

    /// Compare the first two tracks when both have totals.
    fn with_contrast(mut self) -> Self {
        let totals = self.tracks.iter().take(2).filter(|t| t.behind.is_some());
        if totals.count() == 2 {
            self.contrast = Some(Measure::Difference);
        }
        self
    }

    /// Pick this view's genes out of the annotation once it has arrived:
    /// the opened gene alone, else those on this chromosome.
    fn take_genes(&mut self) {
        let Some(models) = self.pending_genes.as_ref().and_then(|m| m.get()) else {
            return;
        };
        self.genes = genes_to_draw(models, &self.title)
            .into_iter()
            .filter(|g| chr_eq(&g.chr, &self.chr))
            .cloned()
            .collect();
        self.pending_genes = None;
        self.plots.iter_mut().for_each(PlotImage::invalidate);
    }

    /// The rows top to bottom: the contrast first, as the main row, then
    /// the tracks.
    fn rows(&self) -> Vec<Row> {
        let contrast = self.contrast.map(|_| Row::Contrast);
        let genes = (!self.genes.is_empty()).then_some(Row::Genes);
        contrast
            .into_iter()
            .chain((0..self.tracks.len()).map(Row::Track))
            .chain(genes)
            .collect()
    }

    fn row_title(&self, row: Row) -> String {
        match (row, self.contrast) {
            (Row::Track(i), _) => {
                let t = &self.tracks[i];
                format!("{} · {}", t.label, t.name)
            }
            (Row::Contrast, Some(m)) => {
                let (a, b) = (self.tracks[0].label, self.tracks[1].label);
                format!("{a} vs {b} · {}", m.name(self.channels.0))
            }
            (Row::Contrast, None) => String::new(),
            (Row::Genes, _) => "genes".into(),
        }
    }

    /// Genes overlapping the window, each with a lane so that genes in one
    /// lane (and their labels) do not collide.
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
            let (px, pw) = ((x + left) as f32, (w - left - right) as f32);
            c.raw(&gene_model_svg(g, &edges, px, pw, mid as f32, GENE_BAND));
        }
    }

    /// Gene lanes in view (at least one).
    fn lanes(&self) -> usize {
        let lanes = self.gene_lanes().iter().map(|(l, _)| l + 1).max();
        lanes.unwrap_or(1)
    }

    /// The contrast per bar over the window (raw sums, never log).
    fn contrast_values(&self, edges: &BinEdges) -> Vec<Option<f64>> {
        let Some(measure) = self.contrast else {
            return Vec::new();
        };
        let counts = |t: &Track| {
            let front = edges.bin(t.front, false);
            let total = edges.bin(t.behind.unwrap_or(t.front), false);
            (front, total)
        };
        let ((ma, na), (mb, nb)) = (counts(&self.tracks[0]), counts(&self.tracks[1]));
        (0..ma.len())
            .map(|k| measure.of((ma[k], na[k]), (mb[k], nb[k])))
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

    fn cursor_col(&self) -> usize {
        self.edges().col_of(self.cursor)
    }

    /// The visible window as a figure, one panel per track.
    fn figure(&self) -> String {
        let (lo, hi) = self.window;
        let panel_h = 170.0;
        let rows = self.rows();
        let height = |row: Row| match row {
            Row::Genes => 24.0 + GENE_LANE * self.lanes() as f64,
            _ => panel_h,
        };
        let total: f64 = rows.iter().map(|&r| height(r)).sum();
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
        let mut y = 44.0;
        for &row in &rows {
            let panel = (8.0, y, 704.0, height(row));
            self.draw_row(row, &mut c, panel, None, self.row_title(row));
            y += height(row);
        }
        c.finish()
    }

    /// A row in the box: a track's bars, or the contrast's signed bars.
    fn draw_row(
        &self,
        row: Row,
        c: &mut Canvas,
        bbox: (f64, f64, f64, f64),
        pointer: Option<usize>,
        title: String,
    ) {
        match (row, self.contrast) {
            (Row::Track(i), _) => self.draw_track(i, c, bbox, pointer, title),
            (Row::Contrast, Some(measure)) => {
                let (x, y, w, h) = bbox;
                let values = self.contrast_values(&self.edges());
                figure::Diverging {
                    values: &values,
                    ticks: self.tick_list(),
                    pointer,
                    marks: Vec::new(),
                    title,
                    x_title: format!("{} position", self.chr),
                    y_title: measure.name(self.channels.0),
                    label: &|v| measure.label(v),
                }
                .draw(c, x, y, w, h);
            }
            (Row::Contrast, None) => {}
            (Row::Genes, _) => self.draw_genes(c, bbox, title),
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
        let ticks = self.tick_list();
        let (front, behind) = t.bin(&edges);
        let marks = t
            .sites_in(lo, hi)
            .into_iter()
            .map(|s| edges.col_of(s))
            .collect();
        let stacked = behind.is_some();
        Bars {
            values: behind.as_ref().unwrap_or(&front),
            front: stacked.then_some(front.as_slice()),
            accent: &|_| stacked,
            y_scale: self.y_scale,
            ticks,
            pointer,
            marks,
            title,
            x_title: format!("{} position", self.chr),
            y_title: t.name.into(),
        }
        .draw(c, x, y, w, h);
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
        if let Some(measure) = self.contrast {
            let (a, b) = (self.tracks[0].label, self.tracks[1].label);
            let value = self.contrast_values(&edges)[col];
            spans.push(dim(format!("   {a} vs {b} ")));
            spans.push(Span::styled(
                value.map_or("-".into(), |v| measure.label(v)),
                HIGHLIGHT,
            ));
        }
        for t in &self.tracks {
            spans.push(dim(format!("   {} ", t.label)));
            spans.push(Span::raw(match t.bin(&edges) {
                (front, Some(behind)) => format!("{:.0}/{:.0}", front[col], behind[col]),
                (front, None) => compact(front[col]),
            }));
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
            KeyCode::Char('c') => self.contrast = self.contrast.map(Measure::next),
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

        // One bar per chart column, re-binned from the raw positions.
        let inner_width = areas[1].width.saturating_sub(2 + GUTTER);
        self.columns = (inner_width as usize).max(1);
        self.clamp_window();
        let edges = self.edges();
        let (lo, hi) = self.window;
        let extra = format!(
            "{}:{}-{}  {} bp/bar · y {}{}",
            self.chr,
            fmt_thousands(lo),
            fmt_thousands(hi),
            self.bin_width(),
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
        let label = |k: i32| ticks.get(k as usize).cloned().flatten();
        let pointer = Some(self.cursor_col() as i32);
        let mut plots = std::mem::take(&mut self.plots);
        plots.resize_with(rows.len(), PlotImage::default);
        for (k, (&row, &area)) in rows.iter().zip(&areas[1..areas.len() - 2]).enumerate() {
            let block = panel(format!(" {} ", self.row_title(row)), true);
            let inner = block.inner(area);
            frame.render_widget(block, area);
            let drawn = self.controls.images().is_some_and(|picker| {
                plots[k].render(frame, inner, picker, |w, h| {
                    let bbox = (0.0, 0.0, w, h);
                    let pointer = Some(self.cursor_col());
                    figure::svg(w, h, |c| {
                        self.draw_row(row, c, bbox, pointer, String::new())
                    })
                })
            });
            if drawn {
                continue;
            }
            if row == Row::Genes {
                self.render_genes(frame.buffer_mut(), inner);
                continue;
            }
            let Row::Track(i) = row else {
                let values = self.contrast_values(&edges);
                let measure = self.contrast.unwrap_or(Measure::Difference);
                let diverging = TextDiverging {
                    values: &values,
                    pointer: self.cursor_col(),
                    label: &|v| measure.label(v),
                    ticks: &ticks,
                };
                diverging.render(frame.buffer_mut(), inner);
                continue;
            };
            let t = &self.tracks[i];
            let (front, behind) = t.bin(&edges);
            let marks = t
                .sites_in(lo, hi)
                .into_iter()
                .map(|s| (edges.col_of(s) as i32, "+", DIM))
                .collect();
            let stacked = behind.is_some();
            HistPlot {
                bins: Binning::with_width(Scale::Linear, 1.0),
                kmin: 0,
                counts: behind.as_ref().unwrap_or(&front),
                style: &|_| if stacked { ACCENTED } else { PLAIN },
                subset: stacked.then_some(front.as_slice()),
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
    }
}

impl PileupView<'_> {
    /// The genes row as text, a line per lane: exons as a heavy line on a
    /// thin intron line with strand arrows, the symbol beside the gene.
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
        if self.contrast.is_some() {
            keys.push(("c", "difference/fold"));
        }
        self.controls.help_keys(&mut keys);
        keys.push(("q", "quit"));
        help_line(&keys)
    }
}

/// Browse full screen until the user quits, asks for the gene list, or
/// searches for something outside the view. `status` shows on the footer
/// first.
pub fn show_pileup(view: View, tracks: Vec<Track>, status: Option<String>) -> anyhow::Result<Exit> {
    let mut browser = PileupView::new(view.title, view.chr, tracks, view.extent);
    browser.channels = view.channels;
    browser.pending_genes = view.genes;
    browser.take_genes();
    browser.status = status.or_else(|| {
        let loading = browser.pending_genes.is_some();
        loading.then(|| "gene models are loading; they appear with the next key".to_string())
    });
    browser.controls = Controls::new(&format!("pileup_{}", view.title)).detect();
    crate::figure::term::run(&mut browser)?;
    Ok(browser.exit.unwrap_or(Exit::Quit))
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
