//! Full-screen site-threshold picker for `faba qc`.
//!
//! One row per [`Criterion`] beside a histogram of the column it cuts, for
//! one editing modality at a time. The thresholds are shared across
//! modalities, as they are on the command line. Every count on screen comes
//! from [`Criterion::fails`], the checks [`SiteFilterArgs::reason`] walks, so
//! the view and the written fileset never disagree. The histogram draws the
//! sites that pass every other threshold over all sites, so it shows
//! what the focused threshold decides among sites that would otherwise be
//! kept.

use data_beans::interactive::ui::{
    header, help_line, input_line, panel, Binned, Binning, HistPlot, Scale, Screen, ACCENTED, DIM,
    HIGHLIGHT, PLAIN,
};
use data_beans::qc::pct;
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::Style;
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;
use rustc_hash::FxHashMap;

use crate::figure::term::PlotImage;
use crate::figure::{self, Anchor, Bars, Canvas, Controls, Key, INK, MUTED};

use crate::site_analysis::metagene::{MetaLayout, MetaModels, REGION_NAMES};
use crate::site_analysis::miami::genemodel::{
    gene_models_from_records, read_records_of, GeneModel,
};
use arrow::record_batch::RecordBatch;
use rustc_hash::FxHashSet;

use super::args::SiteFilterArgs;
use super::layout::{file_name, SITE_MODALITIES};
use super::progress::Progress;
use super::sites::{genomic_sites, Criterion, GeneSites, SiteTable};
use crate::tui::{first_visible, is_go, popup, ShiftEnter, GO_KEYS};

mod annotation;
mod column;
mod draw;
mod genes;

use annotation::*;
use column::*;
use genes::*;

/// Width of HistPlot's y gutter.
use crate::figure::GUTTER;

/// Upper bound on histogram bins, whatever the scale and data.
const MAX_BINS: i32 = 400;

/// Row order on screen.
const SHOWN: [Criterion; 8] = [
    Criterion::MaxPv,
    Criterion::MinLogOdds,
    Criterion::MinFold,
    Criterion::MinCoverage,
    Criterion::MinConverted,
    Criterion::MinEditRatio,
    Criterion::MaxEditRatio,
    Criterion::MinCells,
];

/// The `faba qc` flags that reproduce `f`, as `--flag=value` so negative and
/// infinite values parse.
pub fn qc_flags(f: &SiteFilterArgs) -> String {
    SHOWN
        .iter()
        .map(|c| format!("{}={}", c.flag(), c.fmt(c.get(f))))
        .collect::<Vec<_>>()
        .join(" ")
}

/// One modality's sites, with each site's failed checks as bits.
struct SiteView<'a> {
    table: &'a SiteTable,
    /// Kept cells per site, when a `_site` matrix gave them.
    n_cells: Option<Vec<usize>>,
    /// The criteria that apply, in screen order.
    criteria: Vec<Criterion>,
    /// Per site, a [`Criterion::bit`] for every check it fails.
    fails: Vec<u8>,
    /// Sites by gene, for the gene list; `None` for a table without the
    /// columns.
    genes: Option<GeneSites>,
}

impl<'a> SiteView<'a> {
    fn new(table: &'a SiteTable, n_cells: Option<Vec<usize>>, f: &SiteFilterArgs) -> Self {
        let criteria = SHOWN
            .into_iter()
            .filter(|c| c.applies(table, n_cells.is_some()))
            .collect();
        let mut view = Self {
            table,
            n_cells,
            criteria,
            fails: vec![0; table.len()],
            genes: GeneSites::new(table).ok(),
        };
        for c in view.criteria.clone() {
            view.update(c, f);
        }
        view
    }

    fn cells(&self, i: usize) -> Option<usize> {
        self.n_cells.as_ref().map(|c| c[i])
    }

    /// Recompute `c`'s bit for every site under `f`.
    fn update(&mut self, c: Criterion, f: &SiteFilterArgs) {
        if !self.criteria.contains(&c) {
            return;
        }
        let bit = c.bit();
        for i in 0..self.fails.len() {
            let failed = c.fails(f, self.table, i, self.cells(i));
            self.fails[i] = if failed {
                self.fails[i] | bit
            } else {
                self.fails[i] & !bit
            };
        }
    }
}

/// Counts under the current thresholds, all read off the fail bits.
#[derive(Default)]
struct Tally {
    kept: usize,
    genes: usize,
    /// Sites that fail each criterion, whatever the others say. By
    /// `c as usize`.
    fail: [usize; 8],
    /// Sites that fail the focused criterion and no other: the dropped part
    /// of the histogram's front, and what turning it off would keep.
    only: usize,
    /// Bin counts of the sites that pass every other criterion.
    subset: Vec<usize>,
}

impl Tally {
    fn new(view: &SiteView, focus: Criterion, column: &Column) -> Self {
        let mut t = Tally {
            subset: vec![0; column.hist.counts.len()],
            ..Tally::default()
        };
        let mut gene_seen = vec![false; view.table.n_genes];
        let others = !focus.bit();
        for (i, &m) in view.fails.iter().enumerate() {
            if m == 0 {
                t.kept += 1;
                let g = view.table.gene_id[i] as usize;
                t.genes += usize::from(!gene_seen[g]);
                gene_seen[g] = true;
            } else {
                t.only += usize::from(m == focus.bit());
                for (b, n) in t.fail.iter_mut().enumerate() {
                    *n += usize::from(m >> b & 1 == 1);
                }
            }
            if m & others == 0 {
                t.subset[column.slot[i] as usize] += 1;
            }
        }
        t
    }
}

enum Mode {
    Browse,
    /// Typing an exact raw threshold for the focused criterion.
    Edit(String),
    /// Typing into the gene list's filter.
    Find,
    /// Asking before the thresholds are applied and the fileset written.
    Confirm,
    /// Writing the fileset; keys wait until it is done.
    Writing,
}

/// The panel Tab moves to, which the arrow keys then drive.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Panel {
    Thresholds,
    Genes,
}

/// Writes the fileset once the thresholds are confirmed: `start` spawns the
/// writer with them, and `progress` is how far it has got.
pub struct Writer<'a> {
    /// The output directory, for the confirmation and the progress.
    pub output: &'a str,
    pub progress: &'a Progress,
    pub start: Option<StartWrite<'a>>,
}

/// Starts the write with the applied thresholds and the view's figures.
pub type StartWrite<'a> = Box<dyn FnOnce(SiteFilterArgs, Vec<Figure>) + 'a>;

/// A figure of the view to save with the fileset: a file stem and its SVG.
pub struct Figure {
    pub stem: String,
    pub svg: String,
}

/// Progress that is already over, for a writer with nothing to write.
static NOTHING_TO_WRITE: Progress = Progress::finished();

impl Writer<'_> {
    /// A writer that writes nothing: applying only decides.
    fn idle() -> Self {
        Writer {
            output: "",
            progress: &NOTHING_TO_WRITE,
            start: None,
        }
    }
}

/// State of the picker, independent of the terminal so it can be tested.
struct SitePicker<'a> {
    title: String,
    views: Vec<SiteView<'a>>,
    filter: SiteFilterArgs,
    initial: SiteFilterArgs,
    modality: usize,
    /// Position of the focused criterion in the view's `criteria`.
    focus: usize,
    /// x scale per criterion (`c as usize`), an index into its scales.
    x_scale: [usize; 8],
    y_scale: Scale,
    column: Column,
    tally: Tally,
    mode: Mode,
    panel: Panel,
    /// Shift+Enter reporting, asked for at the first draw.
    shift_enter: ShiftEnter,
    /// Plain Enter was pressed: say how to apply instead.
    enter_hint: bool,
    /// Steps done the last time the write was drawn.
    drawn_steps: usize,
    controls: Controls,
    plot: PlotImage,
    meta: Meta,
    /// What the gene and metagene bars add up.
    weight: Weight,
    /// The focused view's metagene bars: all sites, then the kept ones.
    meta_bars: (Vec<usize>, Vec<usize>),
    meta_plot: PlotImage,
    list: GeneList,
    /// Writes the fileset in the view once the thresholds are applied.
    writer: Writer<'a>,
    decision: Option<Picked>,
}

impl<'a> SitePicker<'a> {
    fn new(title: &str, views: Vec<SiteView<'a>>, filter: SiteFilterArgs, meta: Meta) -> Self {
        assert!(!views.is_empty(), "no site table to pick thresholds on");
        let c = views[0].criteria[0];
        let column = Column::new(&views[0], c, c.scales()[0]);
        let tally = Tally::new(&views[0], c, &column);
        Self {
            title: title.to_string(),
            views,
            initial: filter.clone(),
            filter,
            modality: 0,
            focus: 0,
            x_scale: [0; 8],
            y_scale: Scale::Log,
            column,
            tally,
            mode: Mode::Browse,
            panel: Panel::Thresholds,
            shift_enter: ShiftEnter::default(),
            enter_hint: false,
            drawn_steps: 0,
            controls: Controls::new("qc_sites"),
            plot: PlotImage::default(),
            meta,
            weight: Weight::Sites,
            meta_bars: Default::default(),
            meta_plot: PlotImage::default(),
            list: GeneList::default(),
            writer: Writer::idle(),
            decision: None,
        }
        .pin_genes(&[])
        .with_refreshed_meta()
    }

    fn view(&self) -> &SiteView<'a> {
        &self.views[self.modality]
    }

    fn criterion(&self) -> Criterion {
        self.view().criteria[self.focus]
    }

    fn scale(&self) -> Scale {
        let c = self.criterion();
        c.scales()[self.x_scale[c as usize]]
    }

    fn retally(&mut self) {
        self.tally = Tally::new(self.view(), self.criterion(), &self.column);
        self.refresh_meta();
        self.refresh_genes(false);
    }

    /// Rebuild the focused column, after a modality or focus change.
    fn rebuild_column(&mut self) {
        self.column = Column::new(self.view(), self.criterion(), self.scale());
        self.retally();
    }

    fn set_threshold(&mut self, raw: f64) {
        let c = self.criterion();
        c.set(&mut self.filter, raw);
        for view in &mut self.views {
            view.update(c, &self.filter);
        }
        self.retally();
    }

    fn set_off(&mut self) {
        self.set_threshold(self.criterion().permissive());
    }

    /// Move the threshold one bar right (`dir > 0`) or left on the displayed
    /// axis. Walking out through the kept end turns it off.
    fn step_bin(&mut self, dir: i32) {
        let c = self.criterion();
        let off = c.is_off(&self.filter);
        // An off threshold sits beyond the kept end of the axis.
        let cur = match (off, c.keeps_high()) {
            (false, _) => self.column.key_of(c, c.get(&self.filter)),
            (true, true) => i32::MIN,
            (true, false) => i32::MAX,
        };
        let stops = &self.column.stops;
        let target = if dir > 0 {
            stops.iter().find(|s| s.0 > cur)
        } else {
            stops.iter().rev().find(|s| s.0 < cur)
        };
        match target {
            Some(&(_, raw)) => self.set_threshold(raw),
            None if !off && (dir > 0) != c.keeps_high() => self.set_off(),
            None => {}
        }
    }

    /// Nudge a whole-count threshold by one; others step a bar, `+` always
    /// tightening.
    fn nudge(&mut self, delta: i64) {
        let c = self.criterion();
        if c.is_integer() {
            let v = c.get(&self.filter) as i64 + delta;
            self.set_threshold(v.max(0) as f64);
        } else {
            let dir = if c.keeps_high() { delta } else { -delta };
            self.step_bin(dir as i32);
        }
    }

    fn reset(&mut self) {
        let c = self.criterion();
        self.set_threshold(c.get(&self.initial));
    }

    fn switch_modality(&mut self, delta: isize) {
        let c = self.criterion();
        let n = self.views.len() as isize;
        self.modality = (self.modality as isize + delta).rem_euclid(n) as usize;
        self.refresh_genes(true);
        self.focus = self
            .view()
            .criteria
            .iter()
            .position(|&x| x == c)
            .unwrap_or(0);
        self.rebuild_column();
    }

    fn move_focus(&mut self, delta: isize) {
        let n = self.view().criteria.len() as isize;
        self.focus = (self.focus as isize + delta).rem_euclid(n) as usize;
        self.rebuild_column();
    }

    fn cycle_x_scale(&mut self) {
        let c = self.criterion();
        let i = c as usize;
        self.x_scale[i] = (self.x_scale[i] + 1) % c.scales().len();
        let scale = self.scale();
        self.column.rebin(&self.views[self.modality], c, scale);
        self.retally();
    }

    /// End the session, handing the terminal's keyboard back as it was.
    fn decide(&mut self, picked: Picked) {
        self.shift_enter.release();
        self.decision = Some(picked);
    }

    /// Apply the confirmed thresholds: start the write, which [`Screen::tick`]
    /// then follows to the end.
    fn apply(&mut self) {
        self.mode = Mode::Writing;
    }

    /// Draw the figures and start the writer, once the writing pop-up is up.
    fn start_writing(&mut self) {
        if let Some(start) = self.writer.start.take() {
            let figures = self.figures();
            start(self.filter.clone(), figures);
        }
    }

    /// The view's figure for every modality and every knob, under the
    /// applied thresholds and the scales chosen per knob. The modality on
    /// screen keeps its selected gene; the others show their top gene, the
    /// gene filter aside.
    fn figures(&mut self) -> Vec<Figure> {
        let (modality, focus, at) = (self.modality, self.focus, self.list.at);
        let find = std::mem::take(&mut self.list.find);
        let mut out = Vec::new();
        for m in 0..self.views.len() {
            self.modality = m;
            if m == modality {
                self.list.find.clone_from(&find);
            }
            self.refresh_genes(true);
            self.list.find.clear();
            if m == modality {
                self.list.at = at;
            }
            for f in 0..self.views[m].criteria.len() {
                self.focus = f;
                self.rebuild_column();
                let knob = self.criterion().flag().trim_start_matches("--site-");
                out.push(Figure {
                    stem: format!("{}_{knob}", self.view().table.modality),
                    svg: self.figure(),
                });
            }
        }
        self.modality = modality;
        self.list.find = find;
        self.refresh_genes(true);
        self.list.at = at;
        self.focus = focus;
        self.rebuild_column();
        out
    }

    /// The histogram bin holding the focused threshold, `None` when off.
    fn pointer_key(&self) -> Option<i32> {
        let c = self.criterion();
        (!c.is_off(&self.filter)).then(|| self.column.key_of(c, c.get(&self.filter)))
    }
}

impl Screen for SitePicker<'_> {
    fn done(&self) -> bool {
        self.decision.is_some()
    }

    fn interrupt(&mut self) {
        // A write under way finishes: stopping it would leave half a fileset.
        if !matches!(self.mode, Mode::Writing) {
            self.decide(Picked::Cancelled);
        }
    }

    fn tick(&mut self) -> bool {
        if matches!(self.mode, Mode::Writing) {
            if self.writer.start.is_some() {
                self.start_writing();
                return true;
            }
            let progress = self.writer.progress;
            if progress.is_finished() {
                self.decide(Picked::Apply(self.filter.clone()));
                return false;
            }
            // Redraw only when a step finished.
            let done = progress.done();
            return std::mem::replace(&mut self.drawn_steps, done) != done;
        }
        let arrived = self.meta.poll();
        if arrived {
            self.refresh_meta();
        }
        arrived
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.enter_hint = false;
        if matches!(self.mode, Mode::Browse) {
            match self.controls.key(key) {
                Key::Pass => {}
                Key::Used => return,
                Key::Save(prefix) => return self.controls.save(&self.figure(), &prefix),
            }
        }
        self.plot.invalidate();
        match &mut self.mode {
            Mode::Edit(buf) => match key.code {
                KeyCode::Char(ch)
                    if (ch.is_ascii_digit() || "-+.eEinf".contains(ch)) && buf.len() < 24 =>
                {
                    buf.push(ch)
                }
                KeyCode::Backspace => {
                    buf.pop();
                }
                KeyCode::Enter => {
                    let typed = buf.parse::<f64>().ok();
                    self.mode = Mode::Browse;
                    if let Some(v) = typed {
                        self.set_threshold(v);
                    }
                }
                KeyCode::Esc => self.mode = Mode::Browse,
                _ => {}
            },
            Mode::Find => match key.code {
                KeyCode::Char(ch) if self.list.find.len() < 32 => self.set_find(|f| f.push(ch)),
                KeyCode::Backspace => self.set_find(|f| {
                    f.pop();
                }),
                KeyCode::Enter => self.mode = Mode::Browse,
                KeyCode::Esc => {
                    self.set_find(String::clear);
                    self.mode = Mode::Browse;
                }
                _ => {}
            },
            Mode::Confirm => match key.code {
                _ if is_go(&key) => self.apply(),
                KeyCode::Char('y') => self.apply(),
                KeyCode::Esc | KeyCode::Char('n' | 'q') => self.mode = Mode::Browse,
                _ => {}
            },
            Mode::Writing => {}
            Mode::Browse => match key.code {
                _ if is_go(&key) => self.mode = Mode::Confirm,
                KeyCode::Enter => self.enter_hint = true,
                KeyCode::Tab | KeyCode::BackTab => {
                    self.panel = match self.panel {
                        Panel::Thresholds => Panel::Genes,
                        Panel::Genes => Panel::Thresholds,
                    }
                }
                KeyCode::Char('m') => self.switch_modality(1),
                KeyCode::Char('M') => self.switch_modality(-1),
                KeyCode::Char('[') => self.step_gene(-1),
                KeyCode::Char(']') => self.step_gene(1),
                KeyCode::Char('/') => self.mode = Mode::Find,
                KeyCode::Char('c') => self.switch_weight(),
                KeyCode::Up | KeyCode::Char('k') if self.panel == Panel::Genes => {
                    self.step_gene(-1)
                }
                KeyCode::Down | KeyCode::Char('j') if self.panel == Panel::Genes => {
                    self.step_gene(1)
                }
                KeyCode::PageUp if self.panel == Panel::Genes => self.step_gene(-10),
                KeyCode::PageDown if self.panel == Panel::Genes => self.step_gene(10),
                KeyCode::Up | KeyCode::Char('k') => self.move_focus(-1),
                KeyCode::Down | KeyCode::Char('j') => self.move_focus(1),
                KeyCode::Left | KeyCode::Char('h') => self.step_bin(-1),
                KeyCode::Right | KeyCode::Char('l') => self.step_bin(1),
                KeyCode::Char('-' | ',') => self.nudge(-1),
                KeyCode::Char('+' | '=' | '.') => self.nudge(1),
                KeyCode::Char('o') => self.set_off(),
                KeyCode::Char('r') => self.reset(),
                KeyCode::Char('x') => self.cycle_x_scale(),
                KeyCode::Char('y') => self.y_scale = self.y_scale.next(),
                KeyCode::Char('e') => self.mode = Mode::Edit(String::new()),
                KeyCode::Char(ch) if ch.is_ascii_digit() => self.mode = Mode::Edit(ch.to_string()),
                KeyCode::Char('p') => self.decide(Picked::PrintOnly(self.filter.clone())),
                KeyCode::Char('q') | KeyCode::Esc => self.decide(Picked::Cancelled),
                _ => {}
            },
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        self.shift_enter.arm();
        let [top, tabs, body, footer] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Fill(1),
            Constraint::Length(1),
        ])
        .areas(frame.area());

        let scales = format!(
            "x {} · y {}{}",
            self.scale().name(),
            self.y_scale.name(),
            self.controls.tag()
        );
        frame.render_widget(header("qc", &self.title, &scales), top);

        let mut spans = vec![Span::raw(" ")];
        for (i, v) in self.views.iter().enumerate() {
            let style = if i == self.modality { HIGHLIGHT } else { DIM };
            spans.push(Span::styled(format!("{}  ", v.table.modality), style));
        }
        let tabs_line = Line::from(spans);
        // Keep the mark clear of the header's text and the modality tabs.
        let reserve = (self.title.chars().count() + scales.chars().count() + 12)
            .max(tabs_line.width()) as u16;
        frame.render_widget(tabs_line, tabs);
        let corner = Rect::new(top.x, top.y, top.width, 2);
        crate::figure::logo::draw_mini_logo(frame.buffer_mut(), corner, reserve);

        let [left, right] =
            Layout::horizontal([Constraint::Length(56), Constraint::Fill(1)]).areas(body);
        let table = self.criteria_lines();
        let [table_area, genes] = Layout::vertical([
            Constraint::Length(table.len() as u16 + 2),
            Constraint::Fill(1),
        ])
        .areas(left);
        let block = panel(" site thresholds ".into(), self.panel == Panel::Thresholds);
        let inner = block.inner(table_area);
        frame.render_widget(block, table_area);
        frame.render_widget(Paragraph::new(table), inner);
        self.render_gene_list(frame, genes);
        let [hist, gene, meta] = Layout::vertical([
            Constraint::Fill(3),
            Constraint::Fill(2),
            Constraint::Fill(2),
        ])
        .areas(right);
        self.render_hist(frame, hist);
        self.render_gene(frame, gene);
        self.render_meta(frame, meta);

        let help = match (&self.mode, self.controls.footer()) {
            (_, Some(line)) => line,
            (Mode::Edit(buf), None) => input_line(
                &format!("{} {}: ", self.criterion().label(), self.criterion().flag()),
                buf,
                &[("Enter", "set"), ("Esc", "back")],
            ),
            (Mode::Confirm, None) => help_line(&[
                (GO_KEYS, "apply and write"),
                ("y", "apply and write"),
                ("Esc/n", "back"),
            ]),
            (Mode::Writing, None) => help_line(&[("", "writing; the view closes when done")]),
            (Mode::Find, None) => input_line(
                "gene: ",
                &self.list.find,
                &[("Enter", "keep"), ("Esc", "clear")],
            ),
            (Mode::Browse, None) => {
                let apply = (GO_KEYS, "apply");
                if self.enter_hint {
                    let mut line = help_line(&[apply]);
                    let hint = Span::styled("Enter does nothing here; ", DIM);
                    line.spans.insert(1, hint);
                    return frame.render_widget(line, footer);
                }
                let mut keys = vec![("Tab", "panel")];
                keys.extend(match self.panel {
                    Panel::Thresholds => [("↑/↓", "knob"), ("←/→", "bin")],
                    Panel::Genes => [("↑/↓", "gene"), ("/", "find")],
                });
                // What matters most first: a narrow footer cuts from the end.
                keys.extend([apply, ("p", "print flags"), ("q", "cancel")]);
                self.controls.help_keys(&mut keys);
                keys.extend([
                    ("m", "modality"),
                    ("-/+", "±1"),
                    ("0-9", "type"),
                    ("o", "off"),
                    ("r", "reset"),
                    ("x/y", "scale"),
                ]);
                help_line(&keys)
            }
        };
        frame.render_widget(help, footer);
        match self.mode {
            Mode::Confirm => self.render_confirm(frame, body),
            Mode::Writing => self.render_writing(frame, body),
            _ => {}
        }
    }
}

/// How a picker session ended.
#[derive(Debug, Clone)]
pub enum Picked {
    Cancelled,
    /// Cut and written with these thresholds.
    Apply(SiteFilterArgs),
    /// Print the flags for these thresholds; write nothing.
    PrintOnly(SiteFilterArgs),
}

/// Open the picker over the non-empty site tables, in [`SITE_MODALITIES`]
/// order, with each modality's kept cells per site where `_site` matrices
/// gave them. Needs a terminal and at least one non-empty site table.
pub fn run_site_picker(
    input_dir: &str,
    tables: &FxHashMap<Box<str>, SiteTable>,
    site_cells: &FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>,
    filter: &SiteFilterArgs,
    gff: Option<&str>,
    pinned: &[Box<str>],
    writer: Writer<'_>,
) -> anyhow::Result<Picked> {
    let views: Vec<SiteView> = SITE_MODALITIES
        .iter()
        .filter_map(|m| tables.get(*m))
        .filter(|t| !t.is_empty())
        .map(|t| {
            let n_cells = site_cells.get(&t.modality).map(|acc| t.cells_per_site(acc));
            SiteView::new(t, n_cells, filter)
        })
        .collect();
    anyhow::ensure!(
        !views.is_empty(),
        "no editing site table to pick thresholds on"
    );
    let batches = views.iter().map(|v| v.table.batch.clone()).collect();
    let keys = views
        .iter()
        .filter_map(|v| v.genes.as_ref())
        .flat_map(|g| g.keys.iter().cloned())
        .collect();
    let meta = Meta::start(gff, batches, keys);
    let mut picker =
        SitePicker::new(&file_name(input_dir), views, filter.clone(), meta).pin_genes(pinned);
    picker.writer = writer;
    picker.shift_enter = ShiftEnter::wanted();
    picker.controls = Controls::new("qc_sites").detect();
    data_beans::interactive::ui::run_screen(&mut picker)?;
    Ok(picker.decision.unwrap_or(Picked::Cancelled))
}

#[cfg(test)]
#[path = "tests/site_tui.rs"]
mod tests;
