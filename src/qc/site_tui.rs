//! Full-screen site-threshold picker for `faba qc --interactive`.
//!
//! One row per [`Criterion`] beside a histogram of the column it cuts, for
//! one editing modality at a time. The thresholds are shared across
//! modalities, as they are on the command line. Every count on screen comes
//! from [`Criterion::fails`], the checks [`SiteFilterArgs::reason`] walks, so
//! the view and the written fileset never disagree. The histogram draws the
//! sites that pass every other threshold over all sites, so it shows
//! what the focused threshold decides among sites that would otherwise be
//! kept.

use data_beans::interactive::tui_available;
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
use crate::site_analysis::miami::genemodel::{gene_models_from_records, GeneModel};
use arrow::record_batch::RecordBatch;
use genomic_data::gff::read_gff_record_vec;
use rustc_hash::FxHashSet;

use super::args::SiteFilterArgs;
use super::layout::{file_name, SITE_MODALITIES};
use super::sites::{genomic_sites, Criterion, GeneSites, SiteTable};

mod annotation;
mod column;
mod draw;
mod genes;

use annotation::*;
use column::*;
use genes::*;

/// Width of HistPlot's y gutter.
const GUTTER: u16 = 6;

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
    controls: Controls,
    plot: PlotImage,
    meta: Meta,
    /// What the gene and metagene bars add up.
    weight: Weight,
    /// The focused view's metagene over all sites, which no threshold moves.
    meta_all: Vec<usize>,
    /// The focused view's metagene over the kept sites, refreshed with the
    /// tally.
    meta_kept: Vec<usize>,
    meta_plot: PlotImage,
    /// Per view, its genes (dense ids) in list order: pinned first, then
    /// by number of putative sites.
    gene_order: Vec<Vec<u32>>,
    /// The gene list's filter, matched against symbols.
    gene_find: String,
    /// Position of the selected gene in the filtered list.
    gene_at: usize,
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
            controls: Controls::new("qc_sites"),
            plot: PlotImage::default(),
            meta: Meta::Unavailable(String::new()),
            weight: Weight::Sites,
            meta_all: Vec::new(),
            meta_kept: Vec::new(),
            meta_plot: PlotImage::default(),
            gene_order: Vec::new(),
            gene_find: String::new(),
            gene_at: 0,
            decision: None,
        }
        .with_meta(meta)
        .pin_genes(&[])
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
        self.update_meta();
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
        self.gene_at = 0;
        self.refresh_meta();
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
        self.decision = Some(Picked::Cancelled);
    }

    fn tick(&mut self) -> bool {
        let arrived = self.meta.poll();
        if arrived {
            self.refresh_meta();
        }
        arrived
    }

    fn handle_key(&mut self, key: KeyEvent) {
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
                KeyCode::Char(ch) if self.gene_find.len() < 32 => {
                    self.gene_find.push(ch);
                    self.gene_at = 0;
                }
                KeyCode::Backspace => {
                    self.gene_find.pop();
                    self.gene_at = 0;
                }
                KeyCode::Enter => self.mode = Mode::Browse,
                KeyCode::Esc => {
                    self.gene_find.clear();
                    self.gene_at = 0;
                    self.mode = Mode::Browse;
                }
                _ => {}
            },
            Mode::Browse => match key.code {
                KeyCode::Char('[') => self.step_gene(-1),
                KeyCode::Char(']') => self.step_gene(1),
                KeyCode::Char('/') => self.mode = Mode::Find,
                KeyCode::Char('c') => self.switch_weight(),
                KeyCode::Up | KeyCode::Char('k') => self.move_focus(-1),
                KeyCode::Down | KeyCode::Char('j') => self.move_focus(1),
                KeyCode::Tab => self.switch_modality(1),
                KeyCode::BackTab => self.switch_modality(-1),
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
                KeyCode::Enter => self.decision = Some(Picked::Apply(self.filter.clone())),
                KeyCode::Char('p') => self.decision = Some(Picked::PrintOnly(self.filter.clone())),
                KeyCode::Char('q') | KeyCode::Esc => self.decision = Some(Picked::Cancelled),
                _ => {}
            },
        }
    }

    fn render(&mut self, frame: &mut Frame) {
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
        frame.render_widget(Line::from(spans), tabs);

        let [left, right] =
            Layout::horizontal([Constraint::Length(56), Constraint::Fill(1)]).areas(body);
        let table = self.criteria_lines();
        let [table_area, genes] = Layout::vertical([
            Constraint::Length(table.len() as u16 + 2),
            Constraint::Fill(1),
        ])
        .areas(left);
        let block = panel(" site thresholds ".into(), false);
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
            (Mode::Find, None) => input_line(
                "gene: ",
                &self.gene_find,
                &[("Enter", "keep"), ("Esc", "clear")],
            ),
            (Mode::Browse, None) => {
                let mut keys = vec![
                    ("↑/↓", "knob"),
                    ("←/→", "bin"),
                    ("-/+", "±1"),
                    ("0-9", "type"),
                    ("o", "off"),
                    ("r", "reset"),
                    ("x/y", "scale"),
                    ("Tab", "modality"),
                ];
                self.controls.help_keys(&mut keys);
                keys.extend([("Enter", "apply"), ("p", "print flags"), ("q", "cancel")]);
                help_line(&keys)
            }
        };
        frame.render_widget(help, footer);
    }
}

/// How a picker session ended.
#[derive(Debug, Clone)]
pub enum Picked {
    /// No terminal, or no site table: the thresholds stand as given.
    Skipped,
    Cancelled,
    /// Cut with these thresholds.
    Apply(SiteFilterArgs),
    /// Print the flags for these thresholds; write nothing.
    PrintOnly(SiteFilterArgs),
}

/// Open the picker over the non-empty site tables, in [`SITE_MODALITIES`]
/// order, with each modality's kept cells per site where `_site` matrices
/// gave them.
pub fn run_site_picker(
    input_dir: &str,
    tables: &FxHashMap<Box<str>, SiteTable>,
    site_cells: &FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>,
    filter: &SiteFilterArgs,
    gff: Option<&str>,
    pinned: &[Box<str>],
) -> anyhow::Result<Picked> {
    if !tui_available() {
        log::warn!("--interactive needs stdin and stdout on a terminal; skipping the view");
        return Ok(Picked::Skipped);
    }
    let views: Vec<SiteView> = SITE_MODALITIES
        .iter()
        .filter_map(|m| tables.get(*m))
        .filter(|t| !t.is_empty())
        .map(|t| {
            let n_cells = site_cells.get(&t.modality).map(|acc| t.cells_per_site(acc));
            SiteView::new(t, n_cells, filter)
        })
        .collect();
    if views.is_empty() {
        log::warn!("--interactive: no editing site table to pick thresholds on");
        return Ok(Picked::Skipped);
    }
    let batches = views.iter().map(|v| v.table.batch.clone()).collect();
    let keys = views
        .iter()
        .filter_map(|v| v.genes.as_ref())
        .flat_map(|g| g.keys.iter().cloned())
        .collect();
    let meta = Meta::start(gff, batches, keys);
    let mut picker =
        SitePicker::new(&file_name(input_dir), views, filter.clone(), meta).pin_genes(pinned);
    picker.controls = Controls::new("qc_sites").detect();
    data_beans::interactive::ui::run_screen(&mut picker)?;
    Ok(picker.decision.unwrap_or(Picked::Cancelled))
}

#[cfg(test)]
#[path = "tests/site_tui.rs"]
mod tests;
