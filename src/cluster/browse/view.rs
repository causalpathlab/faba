//! The modules browser, the first three of four levels, each a margin of
//! the next: every gene module against every cell module (left), the sites
//! of the module in focus (right), and the genes the user pins, side by side
//! (below), each on its own axis under its gene model, with an "all cells"
//! track (the margin over cell modules) above one track per cell module.
//! Each gene is in its main gene module, with its posterior-mean reads
//! there. Enter on a pinned gene opens the fourth, its full pileup over
//! every gene module (`super::gene`).
//!
//! Reads are drawn as the pileup draws them, mirror bars (see [`bars`]). A
//! grid can show its blocks as bars, a % heatmap or a z heatmap (see
//! [`grid`]); the strip is in [`strip`]. The sites list is in [`list`], the
//! pins in [`pins`], the cell modules shown in [`columns`], the `/` search
//! in [`find`], keys and the mouse in [`input`], the panes' titles and the
//! footer in [`chrome`], and the fit's trace, full pane on `t`, in
//! [`trace`].

use super::Modules;
use crate::site_analysis::miami::genemodel::GeneModel;
use crate::site_analysis::pileup::SharedModels;
use crate::tui::Hits;
use data_beans::interactive::ui::{header, Screen};
use ratatui::crossterm::event::{KeyEvent, MouseEvent};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::widgets::Paragraph;
use ratatui::Frame;
use rustc_hash::FxHashMap;
use std::path::PathBuf;
use std::sync::Arc;

mod bars;
mod chrome;
mod columns;
mod find;
mod grid;
mod input;
mod list;
mod pins;
mod strip;
mod trace;

#[cfg(test)]
pub use bars::{eighths, hanging};
#[cfg(test)]
pub use grid::{box_size, Stretch};
pub use grid::{Look, Ramp};
#[cfg(test)]
pub use strip::{columns, pack, segment_width};

use find::{Entry, Find};
use grid::{Kept, CELL};
use input::Held;

/// The footer while the gene models load.
const LOADING: &str = "gene models are loading";

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Pane {
    Modules,
    Sites,
    Genes,
}

/// The panes in Tab's order, each a level finer than the last.
const PANES: [Pane; 3] = [Pane::Modules, Pane::Sites, Pane::Genes];

impl Pane {
    /// Its place among [`PANES`], its level less one.
    fn index(self) -> usize {
        self as usize
    }

    /// The pane `by` places along in Tab's order, round from the end.
    fn cycled(self, by: usize) -> Self {
        PANES[(self.index() + by) % PANES.len()]
    }
}

/// What the sites list is sorted by; `o` cycles it, `O` reverses it.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Sort {
    /// The module's order: deepest genes first, each by position.
    Module,
    /// A cell module's methylated %, by its place among the columns.
    Column(usize),
    Reads,
    Cells,
}

/// A pinned gene: its gene role and its place among the role's genes.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct Pin {
    pub role: usize,
    pub gene: usize,
}

/// What a click lands on.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Spot {
    /// A gene role.
    Module(usize),
    /// A member of the module in focus, in the sites grid.
    Site(usize),
    /// A pinned gene's segment of the strip.
    Gene(usize),
    /// A column of pinned gene `k`'s segment, at a member of its role.
    Strip(usize, usize),
    Pane(Pane),
}

/// State of the browser, independent of the terminal so it can be tested.
pub struct ModulesView<'a> {
    m: &'a Modules,
    title: String,
    pub(super) pane: Pane,
    /// The gene role in focus.
    pub(super) module: usize,
    /// The site in focus, an index into the module's members.
    pub(super) site: usize,
    /// The module's members in the sites list's order.
    list: Vec<usize>,
    /// The first site row shown, a place in the list.
    pub(super) from: usize,
    /// Site rows the last frame had room for.
    rows: usize,
    pub(super) pins: Vec<Pin>,
    /// The pinned gene in focus, and its site in focus (a member of its role).
    pub(super) pin: usize,
    pub(super) pin_site: usize,
    /// The first row of the strip shown.
    row_from: usize,
    /// The first cell module shown, a place among the columns, the same in
    /// every pane so their columns line up.
    pub(super) col_from: usize,
    /// Columns per cell module in the grids this frame: [`CELL`], or
    /// [`COMPACT`] when not all fit at that.
    cell_w: u16,
    /// The fewest cell modules a grid had room for in the last frame.
    col_fit: usize,
    /// How the modules grid and the sites grid draw.
    pub(super) looks: [Look; 2],
    pub(super) ramp: Ramp,
    pub(super) sort: Sort,
    /// Sorted high to low (the module's order: as it is).
    pub(super) descending: bool,
    /// The annotation while it is still loading.
    pending: Option<SharedModels>,
    /// The sites' genes' models, by gene, shared with a report being written.
    models: Arc<FxHashMap<Box<str>, GeneModel>>,
    catalog: Vec<Entry>,
    find: Option<Find>,
    held: Option<Held>,
    /// A message for the footer, until the next key.
    status: Option<String>,
    /// The gene whose full pileup Enter asked for.
    open: Option<String>,
    /// Where `e` writes the report, and the modules for its thread.
    report: Option<(PathBuf, Arc<Modules>)>,
    /// The report being written on a thread of its own, to `report`'s
    /// path; `e` again meanwhile writes it anew after, with what has
    /// arrived since.
    writing: crate::report::Writing,
    /// The sites grid's colour ranges, for the gene role and scale they
    /// were worked out for.
    site_stretch: Kept<(usize, Ramp)>,
    /// The modules grid's colour ranges, for the scale they were worked
    /// out for: over every module's blocks, whatever the rows' order.
    pub(super) module_stretch: Kept<Ramp>,
    /// Each gene's widest site label for the gene role it was worked out
    /// for (see [`ModulesView::site_label_w`]).
    pub(super) site_labels: Kept<usize, Vec<usize>>,
    hits: Hits<Spot>,
    /// The fit's trace is open, full pane.
    tracing: bool,
    done: bool,
}

impl<'a> ModulesView<'a> {
    pub fn new(title: &str, m: &'a Modules) -> Self {
        let mut view = Self {
            m,
            title: title.to_string(),
            pane: Pane::Modules,
            module: m.order.first().copied().unwrap_or(0),
            site: 0,
            list: Vec::new(),
            from: 0,
            rows: 1,
            pins: Vec::new(),
            pin: 0,
            pin_site: 0,
            row_from: 0,
            col_from: 0,
            cell_w: CELL,
            col_fit: m.cell_roles,
            looks: [Look::Bars; 2],
            ramp: Ramp::Sqrt,
            sort: Sort::Module,
            descending: false,
            pending: None,
            models: Arc::default(),
            catalog: find::catalog(m),
            find: None,
            held: None,
            status: None,
            open: None,
            report: None,
            writing: Default::default(),
            site_stretch: Kept::default(),
            module_stretch: Kept::default(),
            site_labels: Kept::default(),
            hits: Hits::default(),
            tracing: false,
            done: false,
        };
        view.relist();
        view
    }

    /// Draw gene models from `genes` once they have loaded.
    pub fn with_genes(mut self, genes: Option<SharedModels>) -> Self {
        self.status = genes.is_some().then(|| LOADING.to_string());
        self.pending = genes;
        self.take_genes();
        self
    }

    /// The gene Enter opened, the view ready to run again.
    pub fn take_open(&mut self) -> Option<String> {
        self.done = false;
        self.open.take()
    }

    /// The first cell module shown, a place among the columns.
    pub(super) fn first_column(&self) -> usize {
        self.col_from
    }

    /// Say `msg` on the footer until the next key.
    pub fn set_status(&mut self, msg: String) {
        self.status = Some(msg);
    }

    /// Keep the models of the sites' genes once the annotation has arrived.
    fn take_genes(&mut self) {
        let Some(read) = self.pending.as_ref().and_then(|m| m.get()) else {
            return;
        };
        match read {
            Ok(all) => {
                self.models = Arc::new(super::gene_models(self.m, all));
                if self.status.as_deref() == Some(LOADING) {
                    self.status = None;
                }
            }
            Err(e) => self.status = Some(format!("gene models: {e}")),
        }
        self.pending = None;
    }

    /// Where `e` writes the report of `modules` (those the view shows).
    pub fn with_report(mut self, path: PathBuf, modules: Arc<Modules>) -> Self {
        self.report = Some((path, modules));
        self
    }

    /// `e`: gather and write the report on a thread of its own, saying
    /// where; [`Self::take_written`] says when it is done. Asked again
    /// while one is being written, it is written anew after.
    fn export(&mut self) {
        if let Some(said) = self.start_report() {
            self.status = Some(said);
        }
    }

    /// Start the report, or queue one more while one is under way: what
    /// to say; none without a report to write.
    fn start_report(&mut self) -> Option<String> {
        let (path, modules) = self.report.clone()?;
        if let Some(said) = self.writing.busy() {
            return Some(said);
        }
        let (models, to) = (self.models.clone(), path.clone());
        Some(self.writing.spawn(path, move || {
            super::report::save(&super::report::data(&modules, &models), &to)
        }))
    }

    /// Say where the report was written (or why not) once it is, then
    /// write it again if that was asked meanwhile; true when it says so.
    fn take_written(&mut self) -> bool {
        let Some(done) = self.writing.poll() else {
            return false;
        };
        self.status = Some(done.then(|| self.start_report().unwrap_or_default()));
        true
    }

    /// The gene role in focus; none when no module has sites.
    fn role(&self) -> Option<usize> {
        (!self.m.order.is_empty()).then_some(self.module)
    }

    fn sites(&self) -> usize {
        self.role().map_or(0, |g| self.m.members[g].len())
    }

    /// The gene role and member in focus: the strip's when it has the
    /// focus and genes, else the sites list's.
    fn focus(&self) -> Option<(usize, usize)> {
        match self.pins.get(self.pin) {
            Some(p) if self.pane == Pane::Genes => Some((p.role, self.pin_site)),
            _ => self.role().map(|g| (g, self.site)),
        }
    }

    /// The modules grid's rows: by shift range as bars, diagonal as a heatmap.
    fn module_order(&self) -> &[usize] {
        match self.looks[0] {
            Look::Bars => &self.m.order,
            _ => &self.m.diagonal,
        }
    }

    /// Up or down in the pane in focus.
    fn step(&mut self, by: isize) {
        match self.pane {
            Pane::Modules => {
                let order = self.module_order();
                let Some(last) = order.len().checked_sub(1) else {
                    return;
                };
                let at = order.iter().position(|&g| g == self.module).unwrap_or(0);
                self.select_module(order[crate::tui::moved(at, by, last)]);
            }
            Pane::Sites => self.step_listed(by),
            Pane::Genes => {
                if let Some(last) = self.pins.len().checked_sub(1) {
                    self.focus_pin(crate::tui::moved(self.pin, by, last), None);
                }
            }
        }
    }
}

impl crate::tui::View for ModulesView<'_> {}

impl Screen for ModulesView<'_> {
    fn done(&self) -> bool {
        self.done
    }

    fn interrupt(&mut self) {
        self.done = true;
    }

    fn tick(&mut self) -> bool {
        let waiting = self.pending.is_some();
        self.take_genes();
        let wrote = self.take_written();
        (waiting && self.pending.is_none()) || wrote
    }

    fn takes_mouse(&self) -> bool {
        true
    }

    fn mouse(&mut self, event: MouseEvent) -> bool {
        self.on_mouse(event)
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.on_key(key)
    }

    fn render(&mut self, frame: &mut Frame) {
        self.hits.clear();
        let area = frame.area();
        let m = self.m;
        let mut extra = format!(
            "{} gene modules × {} cell modules",
            m.order.len(),
            m.cell_roles
        );
        if let Some(badge) = self.badge() {
            extra.push_str(&format!(" · {badge}"));
        }
        if self.tracing {
            let [head, body, foot] = Layout::vertical([
                Constraint::Length(1),
                Constraint::Fill(1),
                Constraint::Length(1),
            ])
            .areas(area);
            frame.render_widget(header("modules", &self.title, &extra), head);
            return self.draw_trace(frame, body, foot);
        }
        self.fit_columns(area.width);
        let help = self.footer(area.width);
        // Half the screen for the grids, more when the modules need it,
        // and room left for the strip.
        let modules_h = 2 * m.order.len() as u16 + 3;
        let top_h = (area.height / 2)
            .max(modules_h.min(area.height.saturating_sub(14)))
            .max(4.min(area.height));
        let [head, levels, top, genes, stats, footer] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Length(1),
            Constraint::Length(top_h),
            Constraint::Fill(1),
            Constraint::Length(1),
            Constraint::Length(help.len() as u16),
        ])
        .areas(area);
        let [modules, sites] = Layout::horizontal([
            Constraint::Length(self.split(top.width)),
            Constraint::Fill(1),
        ])
        .areas(top);

        frame.render_widget(header("modules", &self.title, &extra), head);
        frame.render_widget(Paragraph::new(self.levels()), levels);
        // Keep the mark clear of the header's text.
        let reserve = (self.title.chars().count() + extra.chars().count() + 16) as u16;
        let corner = Rect::new(head.x, head.y, head.width, 2);
        crate::figure::logo::draw_mini_logo(frame.buffer_mut(), corner, reserve);

        let inner = self.panels(frame, [modules, sites, genes]);
        let buf = frame.buffer_mut();
        self.draw_modules(buf, inner[0]);
        self.draw_sites(buf, inner[1]);
        self.draw_genes(buf, inner[2]);
        frame.render_widget(Paragraph::new(self.stats_line(stats.width)), stats);
        frame.render_widget(Paragraph::new(help), footer);
    }
}
