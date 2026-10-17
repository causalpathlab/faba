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
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseEvent, MouseEventKind};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::style::Style;
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;
use rustc_hash::FxHashMap;

use crate::figure::term::PlotImage;
use crate::figure::{self, Anchor, Bars, Canvas, Controls, Key, INK, MUTED};

use crate::site_analysis::metagene::{
    region_key, MetaLayout, MetaModels, RegionGlyphs, REGION_COLOURS, REGION_NAMES,
};
use crate::site_analysis::miami::genemodel::{
    gene_models_from_records, read_records_of, GeneModel,
};
use arrow::record_batch::RecordBatch;
use rustc_hash::FxHashSet;

use super::args::SiteFilterArgs;
use super::layout::{file_name, SITE_MODALITIES};
use super::progress::Progress;
use super::sites::{genomic_sites, Criterion, GeneSites, RawColumns, SiteTable};
use crate::tui::browser::{is_annotation, Browser, Nav};
use crate::tui::{
    apply_key, button_popup, enter_hint, first_visible, is_apply, popup, popup_frame, run_view,
    tab_bar, wheel_key, Hits, View, APPLY_KEYS,
};

mod annotation;
mod column;
mod draw;
mod frame;
mod genes;
mod keys;
mod mouse;
mod report;
mod writer;

use annotation::*;
use column::*;
use draw::{render_gff, render_leave, LEGEND};
use genes::*;
use keys::*;
use mouse::*;
pub use writer::{Figure, Writer};

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
    /// Sites by gene, for the gene list, shared with a report being
    /// written; `None` for a table without the columns.
    genes: Option<std::sync::Arc<GeneSites>>,
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
            genes: GeneSites::new(table).ok().map(Into::into),
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
    /// Browsing for an annotation to read.
    Gff(Browser),
    /// Asking before leaving without writing changed thresholds; the key
    /// that asked ([`leave_key`]) leaves as the [`Picked`] says.
    Leave(Picked),
}

/// The panel Tab moves to, which the arrow keys then drive.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Panel {
    Thresholds,
    Genes,
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
    /// Plain Enter was pressed: say how to apply instead.
    enter_hint: bool,
    /// Where the last draw put what, for the mouse.
    hits: Hits<Hit>,
    /// Steps done the last time the write was drawn.
    drawn_steps: usize,
    /// The screen the view last took, which a saved figure takes the shape of.
    drawn: Rect,
    controls: Controls,
    plot: PlotImage,
    meta: Meta,
    /// Where `meta` comes from, to read another annotation in its place.
    annotation: AnnotationSource,
    /// What the gene and metagene bars add up.
    weight: Weight,
    /// The focused view's metagene bars: all sites, then the kept ones.
    meta_bars: (Vec<f64>, Vec<f64>),
    meta_plot: PlotImage,
    list: GeneList,
    /// Writes the fileset in the view once the thresholds are applied.
    writer: Writer<'a>,
    decision: Option<Picked>,
    /// Where `Ctrl+E` writes the report: the working directory unless told.
    report_dir: std::path::PathBuf,
    /// The report being written, on a thread of its own; dropped with
    /// the picker, however its run ends, it finishes first.
    report: crate::report::Writing,
    /// A message for the footer, until the next key.
    status: Option<String>,
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
            enter_hint: false,
            hits: Hits::default(),
            drawn_steps: 0,
            drawn: Rect::default(),
            controls: Controls::new("qc_sites"),
            plot: PlotImage::default(),
            meta,
            annotation: AnnotationSource::default(),
            weight: Weight::Coverage,
            meta_bars: Default::default(),
            meta_plot: PlotImage::default(),
            list: GeneList::default(),
            writer: Writer::idle(),
            decision: None,
            report_dir: std::env::current_dir().unwrap_or_default(),
            report: Default::default(),
            status: None,
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
    /// Leave as `picked` says; first ask, when the thresholds have changed.
    fn leave(&mut self, picked: Picked) {
        if qc_flags(&self.filter) == qc_flags(&self.initial) {
            self.decide(picked);
        } else {
            self.mode = Mode::Leave(picked);
        }
    }

    fn decide(&mut self, picked: Picked) {
        self.decision = Some(picked);
    }

    /// The histogram bin holding the focused threshold, `None` when off.
    fn pointer_key(&self) -> Option<i32> {
        let c = self.criterion();
        (!c.is_off(&self.filter)).then(|| self.column.key_of(c, c.get(&self.filter)))
    }
}

impl View for SitePicker<'_> {
    fn takes_apply(&self) -> bool {
        true
    }

    fn stray_enter(&mut self) {
        self.enter_hint = matches!(self.mode, Mode::Browse);
    }

    /// The report's chord, only while browsing: elsewhere it is dropped,
    /// so it types nothing.
    fn takes_chord(&self, key: &KeyEvent) -> bool {
        is_report(key) && self.browsing()
    }
}

impl Screen for SitePicker<'_> {
    fn takes_mouse(&self) -> bool {
        true
    }

    fn mouse(&mut self, m: MouseEvent) -> bool {
        let acts = crate::tui::is_click_or_wheel(&m);
        if acts {
            self.click(m);
        }
        acts
    }

    fn reports_chords(&self) -> bool {
        true
    }

    fn done(&self) -> bool {
        self.decision.is_some()
    }

    /// Ctrl-C asks first when the thresholds have changed, as `q` does;
    /// a second one leaves. A write under way finishes: stopping it would
    /// leave half a fileset.
    fn interrupt(&mut self) {
        match self.mode {
            Mode::Writing => {}
            Mode::Leave(_) => self.decide(Picked::Cancelled),
            _ => self.leave(Picked::Cancelled),
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
        // Said once browsing again: the other footers would hide it.
        let written = self.browsing() && self.end_report();
        arrived || written
    }

    fn handle_key(&mut self, key: KeyEvent) {
        self.enter_hint = false;
        self.status = None;
        // The chord comes only while browsing (see `takes_chord`).
        if is_report(&key) {
            self.status = Some(self.export());
            return;
        }
        if matches!(self.mode, Mode::Browse) {
            match self.controls.key(key) {
                Key::Pass => {}
                Key::Used => return,
                Key::Save(prefix) => return self.controls.save(&self.figure(), &prefix),
            }
        }
        self.plot.invalidate();
        match self.mode {
            Mode::Gff(_) => self.gff_key(key),
            Mode::Edit(_) => self.edit_key(key),
            Mode::Find => self.find_key(key),
            Mode::Confirm => self.confirm_key(key),
            Mode::Leave(_) => self.leave_key_press(key),
            Mode::Writing => {}
            Mode::Browse => self.browse_key(key),
        }
    }

    fn render(&mut self, frame: &mut Frame) {
        self.render_frame(frame);
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
/// Returns how it ended and the annotation it ended on, which `g` may have
/// changed from `gff`.
pub fn run_site_picker(
    input_dir: &str,
    tables: &FxHashMap<Box<str>, SiteTable>,
    site_cells: &FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>,
    filter: &SiteFilterArgs,
    gff: Option<&str>,
    pinned: &[Box<str>],
    writer: Writer<'_>,
) -> anyhow::Result<(Picked, Option<Box<str>>)> {
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
    let annotation = AnnotationSource {
        gff: gff.map(Box::from),
        given: gff.map(Box::from),
        dir: std::path::PathBuf::from(input_dir),
        batches: views.iter().map(|v| v.table.batch.clone()).collect(),
        keys: views
            .iter()
            .filter_map(|v| v.genes.as_ref())
            .flat_map(|g| g.keys.iter().cloned())
            .collect(),
    };
    let meta = Meta::start(&annotation);
    let mut picker =
        SitePicker::new(&file_name(input_dir), views, filter.clone(), meta).pin_genes(pinned);
    picker.annotation = annotation;
    picker.writer = writer;
    picker.controls = Controls::new("qc_sites").detect();
    run_view(&mut picker)?;
    let gff = picker.gff().map(Box::from);
    Ok((picker.decision.unwrap_or(Picked::Cancelled), gff))
}

#[cfg(test)]
#[path = "tests/site_tui/mod.rs"]
mod tests;
