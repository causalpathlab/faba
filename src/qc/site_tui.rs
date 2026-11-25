//! Full-screen site-threshold picker for `faba qc`, and the same view read
//! only for `faba qc-view` (see [`run_site_viewer`]).
//!
//! One row per [`Criterion`] beside a histogram of the column it cuts, for
//! one editing modality at a time. The thresholds are shared across
//! modalities, as they are on the command line. Every count on screen comes
//! from [`Criterion::fails`], the checks [`SiteFilterArgs::reason`] walks, so
//! the view and the written fileset never disagree. The histogram draws the
//! sites that pass every other threshold over all sites, so it shows
//! what the focused threshold decides among sites that would otherwise be
//! kept.

use crate::figure::style::Sheet;
use crate::tui::ui::{
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
use super::pairs::{self, Axis, Sites};
use super::progress::Progress;
use super::sites::{genomic_sites, Criterion, GeneSites, RawColumns, SiteTable};
use crate::tui::browser::{self, is_annotation, Browser, Outcome};
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
mod pairs_view;
mod report;
mod writer;

use annotation::*;
use column::*;
pub(crate) use draw::{knob_least, knob_room};
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
    /// The 2D view's panels: each pair of axes, and why it is shown.
    pairs: Vec<(Axis, Axis, String)>,
    /// Each axis' values, read once.
    cache: pairs::Cache,
    /// Bumped whenever a fail bit may have changed, so a panel binned at
    /// an earlier one is binned again.
    generation: u64,
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
            pairs: Vec::new(),
            cache: pairs::Cache::default(),
            generation: 0,
        };
        for c in view.criteria.clone() {
            view.update(c, f);
        }
        view
    }

    /// The panels shown: those chosen, else Auto's at the thresholds now
    /// (picked only when something shows them).
    fn shown_pairs(&self) -> Vec<(Axis, Axis, String)> {
        match self.pairs.is_empty() {
            true => pairs::auto_panels(&self.sites(), pairs::PANELS),
            false => self.pairs.clone(),
        }
    }

    /// Its sites as the 2D panels read them.
    fn sites(&self) -> Sites<'_> {
        Sites {
            table: self.table,
            n_cells: self.n_cells.as_deref(),
            fails: &self.fails,
            criteria: &self.criteria,
            cache: &self.cache,
        }
    }

    /// The panels Auto picks at the thresholds as they are.
    fn auto_pairs(&mut self) {
        self.pairs = pairs::auto_panels(&self.sites(), pairs::PANELS);
    }

    fn cells(&self, i: usize) -> Option<usize> {
        self.n_cells.as_ref().map(|c| c[i])
    }

    /// Recompute `c`'s bit for every site under `f`, in parallel; one pass,
    /// whatever the other criteria.
    fn update(&mut self, c: Criterion, f: &SiteFilterArgs) {
        use rayon::prelude::*;
        if !self.criteria.contains(&c) {
            return;
        }
        let bit = c.bit();
        let (table, n_cells) = (self.table, self.n_cells.as_deref());
        let changed = (self.fails.par_iter_mut().enumerate())
            .with_min_len(1 << 15)
            .map(|(i, m)| {
                let failed = c.fails(f, table, i, n_cells.map(|n| n[i]));
                let was = *m;
                *m = if failed { *m | bit } else { *m & !bit };
                *m != was
            })
            .reduce(|| false, |a, b| a || b);
        // A panel binned at an earlier generation is binned again.
        self.generation += u64::from(changed);
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
    /// Sites that fail each criterion and no other, by `c as usize`: for
    /// the focused one, the dropped part of the histogram's front, and what
    /// turning it off would keep.
    alone: [usize; 8],
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
                if m.is_power_of_two() {
                    t.alone[m.trailing_zeros() as usize] += 1;
                }
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
    Gff(Box<Browser>),
    /// Asking before leaving without writing changed thresholds; the key
    /// that asked ([`leave_key`]) leaves as the [`Picked`] says.
    Leave(Picked),
}

/// What a 2D panel is binned from: the modality, its axes, the fail bits'
/// generation and the squares a side.
type GridKey = (usize, Axis, Axis, u64, usize);

/// The panel Tab moves to, which the arrow keys then drive.
#[derive(Clone, Copy, PartialEq, Eq, Debug)]
enum Panel {
    Thresholds,
    Genes,
    /// The 2D view's panels, when it is shown.
    Pairs,
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
    /// The 2D view in place of the histogram, gene and metagene (`v`), its
    /// focused panel, and each panel's image.
    two_d: bool,
    pair_at: usize,
    pair_plots: [PlotImage; pairs::PANELS],
    /// What each panel's image was last drawn of, and the panels binned
    /// lately, so a draw that changes nothing bins nothing.
    pair_drawn: [Option<GridKey>; pairs::PANELS],
    pair_grids: Vec<(GridKey, std::sync::Arc<pairs::Grid>)>,
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
    /// Each batch's cell QC, for the report.
    cells: CellsShared,
    /// Read only (`faba qc-view`): nothing to apply, `q` leaves at once,
    /// and `s` saves the panel on screen as a paper panel, as these panels
    /// and this sheet lay it out (see [`Self::screen_panel`]); `None` in
    /// `faba qc`.
    viewing: Option<(crate::qc::site_plot::Plan, Sheet)>,
}

impl<'a> SitePicker<'a> {
    fn new(title: &str, views: Vec<SiteView<'a>>, filter: SiteFilterArgs, meta: Meta) -> Self {
        assert!(!views.is_empty(), "no site table to pick thresholds on");
        let c = views[0].criteria[0];
        let x_scale = knob_scales(&views[0].cache.scales);
        let column = Column::new(&views[0], c, c.scales()[x_scale[c as usize]]);
        let tally = Tally::new(&views[0], c, &column);
        Self {
            title: title.to_string(),
            views,
            initial: filter.clone(),
            filter,
            modality: 0,
            focus: 0,
            x_scale,
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
            two_d: false,
            pair_at: 0,
            pair_plots: Default::default(),
            pair_drawn: Default::default(),
            pair_grids: Vec::new(),
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
            cells: CellsShared::default(),
            viewing: None,
        }
        .pin_genes(&[])
        .with_refreshed_meta()
    }

    fn view(&self) -> &SiteView<'a> {
        &self.views[self.modality]
    }

    /// Whether this is `faba qc-view`'s view, which writes no fileset.
    fn read_only(&self) -> bool {
        self.viewing.is_some()
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

    /// `s`, once a name is confirmed: the view's figure at its own size,
    /// or read only the panel on screen (see [`Self::screen_panel`]),
    /// saved under `prefix` (never over a file: the prompt refuses one).
    fn save(&mut self, prefix: &str) {
        let svg = match &self.viewing {
            None => Ok(self.figure(Sheet::OWN).expect("its own size fits")),
            Some((plan, sheet)) => self.screen_panel(plan, *sheet),
        };
        match svg {
            Ok(svg) => self.controls.save(&svg, prefix),
            Err(e) => self.status = Some(format!("not saved: {e}")),
        }
    }

    /// The panel on screen as the paper panel `faba site-plot` draws of
    /// it, laid out for `sheet` with `plan`'s marginals: the focused
    /// pair in 2D, else the focused knob's histogram, each on the view's
    /// scales and bins, at the thresholds now.
    pub(super) fn screen_panel(
        &self,
        plan: &crate::qc::site_plot::Plan,
        sheet: Sheet,
    ) -> anyhow::Result<String> {
        use crate::qc::site_plot::{Choice, Panels, Plan};
        let v = self.view();
        let pair = (self.two_d)
            .then(|| v.shown_pairs().get(self.pair_at).map(|p| (p.0, p.1)))
            .flatten();
        let (hist, pairs) = match pair {
            Some(p) => (Vec::new(), vec![p]),
            None => (vec![Axis::of(self.criterion())], Vec::new()),
        };
        let plan = Plan {
            hist: Choice::Given(hist),
            pairs: Choice::Given(pairs),
            marginals: plan.marginals,
            scales: v.cache.scales,
            bins: Some(v.cache.bins),
        };
        let panels = Panels::of_view(&v.sites(), &self.filter, &plan, sheet)?;
        (panels.figures().pop())
            .map(|f| f.svg)
            .ok_or_else(|| anyhow::anyhow!("{} has no such panel", v.table.modality))
    }
}

impl View for SitePicker<'_> {
    /// Read only, there is nothing to apply.
    fn takes_apply(&self) -> bool {
        !self.read_only()
    }

    fn stray_enter(&mut self) {
        self.enter_hint = matches!(self.mode, Mode::Browse) && !self.read_only();
    }

    /// The report's chord, only while browsing: elsewhere it is dropped,
    /// so it types nothing.
    fn takes_chord(&self, key: &KeyEvent) -> bool {
        is_report(key) && self.browsing()
    }

    /// None while a threshold, a gene or a name is typed, an annotation is
    /// browsed for, or the fileset is written.
    fn help(&self) -> Option<crate::tui::KeyList> {
        if self.controls.prompting() {
            return None;
        }
        let keys: Vec<(&str, &str)> = match &self.mode {
            Mode::Confirm => vec![(APPLY_KEYS, "apply and write"), ("Esc/n", "back")],
            Mode::Leave(picked) => vec![
                match picked {
                    Picked::PrintOnly(_) => ("p", "print the flags and leave"),
                    _ => ("q", "leave without writing"),
                },
                ("Esc/n", "back"),
            ],
            Mode::Browse => self.browse_help(),
            Mode::Edit(_) | Mode::Find | Mode::Gff(_) | Mode::Writing => return None,
        };
        Some(crate::tui::keys(&keys))
    }
}

impl Screen for SitePicker<'_> {
    /// `faba kill`: leave without asking, changed thresholds dropped and
    /// nothing written; a write under way finishes first, the view done
    /// once it is.
    fn quit(&mut self) {
        if !matches!(self.mode, Mode::Writing) {
            self.decide(Picked::Cancelled);
        }
    }

    /// For `faba sessions` and the session picker.
    fn status(&self) -> Option<String> {
        Some(match self.mode {
            Mode::Writing => format!("writing {}", self.writer.output),
            _ if self.read_only() => format!("viewing the site thresholds of {}", self.title),
            _ => format!("picking site thresholds for {}", self.title),
        })
    }

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
    /// leave half a fileset. Read only, it leaves at once.
    fn interrupt(&mut self) {
        match self.mode {
            Mode::Writing => {}
            _ if self.read_only() => self.decide(Picked::Cancelled),
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
        let listed = match &mut self.mode {
            Mode::Gff(b) => crate::tui::browser::tick(b),
            _ => false,
        };
        let arrived = self.meta.poll();
        if arrived {
            self.refresh_meta();
        }
        // Said once browsing again: the other footers would hide it.
        let written = self.browsing() && self.end_report();
        listed || arrived || written
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
                Key::Save(prefix) => return self.save(&prefix),
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
        self.controls.refresh();
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

/// The non-empty site tables as views, in [`SITE_MODALITIES`] order, with
/// each modality's kept cells per site where `_site` matrices gave them.
/// Their 2D axes are read on `scales`' scales, their rulers in about its
/// bins (see [`pairs::Cache`]); on the defaults without.
fn site_views<'a>(
    tables: &'a FxHashMap<Box<str>, SiteTable>,
    site_cells: &FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>,
    filter: &SiteFilterArgs,
    scales: Option<(crate::qc::scale::Scales, usize)>,
) -> Vec<SiteView<'a>> {
    SITE_MODALITIES
        .iter()
        .filter_map(|m| tables.get(*m))
        .filter(|t| !t.is_empty())
        .map(|t| {
            let n_cells = site_cells.get(&t.modality).map(|acc| t.cells_per_site(acc));
            let mut v = SiteView::new(t, n_cells, filter);
            if let Some((scales, bins)) = scales {
                v.cache = pairs::Cache::new(scales, bins);
            }
            v
        })
        .collect()
}

/// The picker over the site tables, its annotation being read; its knobs'
/// histograms and 2D axes on `scales` (see [`site_views`]).
fn picker<'a>(
    input_dir: &str,
    tables: &'a FxHashMap<Box<str>, SiteTable>,
    site_cells: &FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>,
    filter: &SiteFilterArgs,
    (gff, pinned): (Option<&str>, &[Box<str>]),
    scales: Option<(crate::qc::scale::Scales, usize)>,
) -> anyhow::Result<SitePicker<'a>> {
    let views = site_views(tables, site_cells, filter, scales);
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
    Ok(picker)
}

#[cfg(test)]
impl SitePicker<'_> {
    /// Every view's 2D axes read again on `scales`, their rulers in about
    /// `bins` bins (see [`pairs::Cache`]); each knob's histogram back on
    /// its quantity's scale.
    pub(super) fn set_scales(&mut self, scales: crate::qc::scale::Scales, bins: usize) {
        for v in &mut self.views {
            v.cache = pairs::Cache::new(scales, bins);
            v.generation += 1;
        }
        self.x_scale = knob_scales(&scales);
        self.rebuild_column();
    }
}

/// Each batch's cell QC, shared by the view, its report and its figures.
pub type CellsShared = std::sync::Arc<[crate::qc::cells_plot::Cells]>;

/// Each batch's cell QC as a figure, laid out for `sheet`: its own size
/// grown to the font, a side asked for checked as a knob figure's is
/// ([`knob_least`]).
pub(crate) fn cell_figures(cells: &[crate::qc::cells_plot::Cells], sheet: Sheet) -> Vec<Figure> {
    let size = sheet.fill_own(crate::qc::cells_plot::SIZE);
    let style = sheet.style(size);
    (cells.iter())
        .map(|c| Figure {
            stem: c.stem(),
            svg: c.figure(size, &style),
        })
        .collect()
}

/// What a run without the view takes from the site tables and the cells,
/// built once: the dry run's counts and Auto's pairs, the figures, the
/// page. Its picker is there when some site table is.
pub struct Batch<'a> {
    picker: Option<SitePicker<'a>>,
    cells: CellsShared,
}

impl<'a> Batch<'a> {
    /// Over `tables` (and their kept cells per site) at `filter`, the
    /// annotation `gff` read for the metagene, `pinned` genes first, with
    /// each batch's `cells`; its knobs' histograms and 2D panels' axes on
    /// `scales` (`--scale`, its rulers in about `bins` bins), else the
    /// defaults.
    pub fn new(
        input_dir: &str,
        tables: &'a FxHashMap<Box<str>, SiteTable>,
        site_cells: &FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>,
        filter: &SiteFilterArgs,
        genes: (Option<&str>, &[Box<str>]),
        (cells, scales): (CellsShared, Option<(crate::qc::scale::Scales, usize)>),
    ) -> anyhow::Result<Self> {
        let has_sites = tables.values().any(|t| !t.is_empty());
        let picker = match has_sites {
            true => {
                let mut p = picker(input_dir, tables, site_cells, filter, genes, scales)?;
                p.cells = cells.clone();
                Some(p)
            }
            false => None,
        };
        Ok(Self { picker, cells })
    }

    /// Every applicable threshold's counts, and Auto's pairs, modality by
    /// modality.
    pub fn counts(&self) -> (Vec<ThresholdCount>, Vec<(Box<str>, pairs::Pair)>) {
        match &self.picker {
            Some(p) => counts_of(&p.views, &p.filter),
            None => Default::default(),
        }
    }

    /// The figures the view writes under `qc_plots/` on apply, laid out for
    /// `sheet`: every modality and knob, the metagene once the annotation is
    /// read, the pairs, and each batch's cells; then the site panels `plan`
    /// asks for, if any (see [`SitePicker::written_figures`]).
    pub fn figures(
        &mut self,
        plan: Option<&crate::qc::site_plot::Plan>,
        sheet: Sheet,
    ) -> anyhow::Result<Vec<Figure>> {
        match &mut self.picker {
            Some(p) => {
                p.meta.wait();
                if let Meta::Unavailable(why) = &p.meta {
                    log::warn!("qc_plots without a metagene: {why}");
                }
                p.written_figures(plan, sheet)
            }
            None => Ok(cell_figures(&self.cells, sheet)),
        }
    }

    /// The page `Ctrl+E` writes, to `path` (never over a file); false, and
    /// nothing written, without a site table.
    pub fn write_page(&mut self, path: &std::path::Path) -> anyhow::Result<bool> {
        let Some(p) = &mut self.picker else {
            log::warn!("no editing site table: --html writes no page");
            return Ok(false);
        };
        p.meta.wait();
        let title = format!("faba qc · {}", p.title);
        crate::report::write_new(path, &title, report::FRAGMENT, &p.report_source().data())?;
        Ok(true)
    }
}

/// Each step of the 2D view timed on real tables (see the ignored
/// `qc_timing` test), in seconds: building each modality's fail masks and
/// Auto's pairs, the pairs alone, binning four panels, the view's first
/// draw, a threshold's keypress, an axis' step, the figures and the page,
/// written into `out`.
#[cfg(test)]
pub(crate) fn time_steps(
    tables: &FxHashMap<Box<str>, SiteTable>,
    site_cells: &FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>,
    filter: &SiteFilterArgs,
    cells: CellsShared,
    out: &std::path::Path,
) -> Vec<(String, f64)> {
    use ratatui::crossterm::event::KeyModifiers;
    use std::time::Instant;
    let mut said = Vec::new();
    let mut time = |what: String, t: Instant| said.push((what, t.elapsed().as_secs_f64()));
    let t = Instant::now();
    let views = site_views(tables, site_cells, filter, None);
    let sizes: Vec<String> = views
        .iter()
        .map(|v| format!("{} {}", v.table.modality, v.table.len()))
        .collect();
    time(
        format!("fail masks + Auto, every modality ({})", sizes.join(", ")),
        t,
    );
    let big = (0..views.len())
        .max_by_key(|&m| views[m].table.len())
        .unwrap_or(0);
    let s = views[big].sites();
    let t = Instant::now();
    let auto = pairs::auto_pairs(&s, pairs::PANELS);
    time(
        format!("Auto's pairs alone, {}", views[big].table.modality),
        t,
    );
    let t = Instant::now();
    for p in &auto {
        pairs::grid(&s, (p.x, p.y), filter, pairs::BINS);
    }
    time("binning four panels".into(), t);
    let mut p = SitePicker::new(
        "timing",
        views,
        filter.clone(),
        Meta::Unavailable("timing".into()),
    );
    p.modality = big;
    p.cells = cells;
    let mut term = ratatui::Terminal::new(ratatui::backend::TestBackend::new(200, 60)).unwrap();
    let key = |c| KeyEvent::new(c, KeyModifiers::NONE);
    let t = Instant::now();
    crate::tui::feed(&mut p, key(KeyCode::Char('d')));
    term.draw(|f| p.render(f)).unwrap();
    time("opening the 2D view (d, one draw)".into(), t);
    while p.criterion() != Criterion::MinCoverage {
        crate::tui::feed(&mut p, key(KeyCode::Down));
    }
    let t = Instant::now();
    crate::tui::feed(&mut p, key(KeyCode::Char('+')));
    term.draw(|f| p.render(f)).unwrap();
    time("a threshold's keypress (+, one draw)".into(), t);
    crate::tui::feed(&mut p, key(KeyCode::Tab));
    crate::tui::feed(&mut p, key(KeyCode::Tab));
    let t = Instant::now();
    crate::tui::feed(&mut p, key(KeyCode::Char('x')));
    term.draw(|f| p.render(f)).unwrap();
    time("a panel's axis (x, one draw)".into(), t);
    let t = Instant::now();
    let figures = p.figures(Sheet::OWN).unwrap();
    time(format!("the view's {} figures, drawn", figures.len()), t);
    let t = Instant::now();
    let html = out.join("qc.html");
    let title = "faba qc · timing".to_string();
    crate::report::write(&html, &title, report::FRAGMENT, &p.report_source().data()).unwrap();
    let size = std::fs::metadata(&html).map_or(0, |m| m.len());
    time(
        format!("the page, written ({:.1} MB)", size as f64 / 1e6),
        t,
    );
    said
}

/// The 2D figures, the cell figures and the page of `tables` (each with
/// its kept cells per site) and `cells`, at `filter`, written into `out`
/// as SVG and PNG: for redrawing them quickly (see the ignored `qc_redraw`
/// test).
#[cfg(test)]
pub(crate) fn redraw(
    tables: &[(SiteTable, Option<Vec<usize>>)],
    filter: &SiteFilterArgs,
    cells: Vec<crate::qc::cells_plot::Cells>,
    out: &std::path::Path,
) -> Vec<std::path::PathBuf> {
    let views: Vec<SiteView> = (tables.iter())
        .map(|(t, n)| SiteView::new(t, n.clone(), filter))
        .collect();
    let mut p = SitePicker::new(
        "redraw",
        views,
        filter.clone(),
        Meta::Unavailable("redraw".into()),
    );
    let mut written = Vec::new();
    let mut save = |stem: &str, svg: &str| {
        let prefix = out.join(stem).to_string_lossy().into_owned();
        written.extend(
            crate::figure::save_with_png(svg, &prefix)
                .unwrap()
                .into_iter()
                .map(Into::into),
        );
    };
    for m in 0..p.views.len() {
        let stem = pairs_view::pairs_stem(&p.views[m].table.modality);
        save(&stem, &p.pairs_figure(m, Sheet::OWN));
    }
    for c in &cells {
        let svg = c.figure((720.0, 320.0), &Default::default());
        save(&format!("{}_cells", c.batch), &svg);
    }
    p.cells = cells.into();
    let html = out.join("qc.html");
    let title = "faba qc · redraw".to_string();
    crate::report::write(&html, &title, report::FRAGMENT, &p.report_source().data()).unwrap();
    written.push(html);
    written
}

/// What one site threshold does to one modality, for `faba qc --dry-run`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ThresholdCount {
    pub modality: Box<str>,
    pub flag: &'static str,
    /// The threshold as given, `off` when it is not set.
    pub threshold: String,
    /// Sites failing it, whatever the others say.
    pub removes: usize,
    /// Sites failing it and no other.
    pub only: usize,
    /// Sites passing every threshold.
    pub kept: usize,
    pub sites: usize,
}

/// Every applicable threshold's counts at `filter`, modality by modality.
#[cfg(test)]
pub fn threshold_counts(
    tables: &FxHashMap<Box<str>, SiteTable>,
    site_cells: &FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>,
    filter: &SiteFilterArgs,
) -> Vec<ThresholdCount> {
    counts_of(&site_views(tables, site_cells, filter, None), filter).0
}

/// Every applicable threshold's counts at `filter` over `views`, and the
/// pairs Auto picks, modality by modality (for `faba qc --dry-run`).
fn counts_of(
    views: &[SiteView],
    filter: &SiteFilterArgs,
) -> (Vec<ThresholdCount>, Vec<(Box<str>, pairs::Pair)>) {
    let mut pairs_of = Vec::new();
    let mut out = Vec::new();
    for v in views {
        let m = &v.table.modality;
        for p in pairs::auto_pairs(&v.sites(), pairs::PANELS) {
            pairs_of.push((m.clone(), p));
        }
        let kept = v.fails.iter().filter(|&&f| f == 0).count();
        for &c in &v.criteria {
            let bit = c.bit();
            out.push(ThresholdCount {
                modality: m.clone(),
                flag: c.flag(),
                threshold: c.shown(filter).unwrap_or_else(|| "off".into()),
                removes: v.fails.iter().filter(|&&f| f & bit != 0).count(),
                only: v.fails.iter().filter(|&&f| f == bit).count(),
                kept,
                sites: v.fails.len(),
            });
        }
    }
    (out, pairs_of)
}

/// Open the picker over the non-empty site tables (their knobs and 2D
/// axes on `scales`, see [`site_views`]), in [`SITE_MODALITIES`]
/// order, with each modality's kept cells per site where `_site` matrices
/// gave them. Needs a terminal and at least one non-empty site table.
/// Returns how it ended and the annotation it ended on, which `a` may have
/// changed from `gff`.
pub fn run_site_picker(
    input_dir: &str,
    tables: &FxHashMap<Box<str>, SiteTable>,
    site_cells: &FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>,
    filter: &SiteFilterArgs,
    (gff, pinned): (Option<&str>, &[Box<str>]),
    (scales, bins): (crate::qc::scale::Scales, usize),
    writer: Writer<'_>,
) -> anyhow::Result<(Picked, Option<Box<str>>)> {
    let mut picker = picker(
        input_dir,
        tables,
        site_cells,
        filter,
        (gff, pinned),
        Some((scales, bins)),
    )?;
    picker.page_problem(&writer.sheet)?;
    picker.cells = writer.cells.clone();
    picker.writer = writer;
    picker.controls = Controls::new("qc_sites").detect();
    match run_view(&mut picker) {
        // Stopped by `faba kill` once the fileset was written: it is
        // recorded as any other.
        Err(e)
            if crate::tui::is_stopped(&e) && matches!(picker.decision, Some(Picked::Apply(_))) => {}
        shown => shown?,
    }
    let gff = picker.gff().map(Box::from);
    Ok((picker.decision.unwrap_or(Picked::Cancelled), gff))
}

/// Open the picker read only (`faba qc-view`) over the non-empty site
/// tables at `filter`, as [`run_site_picker`] opens it but with nothing to
/// apply: its knobs and 2D axes on `plan`'s scales, `s` saving the panel
/// on screen as `plan` and `sheet` lay it out, each batch's `cells` on its
/// page. Needs a terminal and at least one non-empty site table.
pub fn run_site_viewer(
    input_dir: &str,
    (tables, site_cells): (&FxHashMap<Box<str>, SiteTable>, &super::sites::SiteCells),
    filter: &SiteFilterArgs,
    genes: (Option<&str>, &[Box<str>]),
    (plan, sheet): (crate::qc::site_plot::Plan, Sheet),
    cells: CellsShared,
) -> anyhow::Result<()> {
    let scales = (plan.scales, plan.bins_on(&sheet));
    let mut picker = picker(input_dir, tables, site_cells, filter, genes, Some(scales))?;
    picker.cells = cells;
    picker.viewing = Some((plan, sheet));
    picker.controls = Controls::new("qc_view").detect();
    run_view(&mut picker)
}

#[cfg(test)]
#[path = "tests/site_tui/mod.rs"]
mod tests;
