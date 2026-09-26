//! Full-screen site-threshold picker for `faba qc --interactive` and
//! `faba qc-report --interactive`.
//!
//! One row per [`Criterion`] beside a histogram of the column it cuts, for
//! one editing modality at a time. The thresholds are shared across
//! modalities, as they are on the command line. Every count on screen comes
//! from [`Criterion::fails`], the checks [`SiteFilterArgs::reason`] walks, so
//! the view and the written fileset never disagree. The histogram draws the
//! sites that pass every other threshold in front of all sites, so it shows
//! what the focused threshold decides among sites that would otherwise be
//! kept.

use data_beans::interactive::tui_available;
use data_beans::interactive::ui::{
    header, help_line, input_line, panel, run_screen, Binned, Binning, HistPlot, Scale, Screen,
    ACCENTED, DIM, HIGHLIGHT, PLAIN,
};
use data_beans::qc::pct;
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::layout::{Constraint, Layout, Rect};
use ratatui::text::{Line, Span};
use ratatui::widgets::Paragraph;
use ratatui::Frame;
use rustc_hash::FxHashMap;

use super::args::SiteFilterArgs;
use super::layout::{file_name, SITE_MODALITIES};
use super::sites::{Criterion, SiteTable};

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

/// How each criterion is drawn.
impl Criterion {
    fn label(self) -> &'static str {
        match self {
            Criterion::MaxPv => "p-value",
            Criterion::MinLogOdds => "log odds",
            Criterion::MinFold => "fold",
            Criterion::MinCoverage => "coverage",
            Criterion::MinConverted => "converted",
            Criterion::MinEditRatio => "edit ratio ≥",
            Criterion::MaxEditRatio => "edit ratio ≤",
            Criterion::MinCells => "cells",
        }
    }

    /// What the histogram's x axis shows.
    fn axis(self) -> &'static str {
        match self {
            Criterion::MaxPv => "-log10 p",
            Criterion::MinLogOdds => "raw log odds ratio",
            Criterion::MinFold => "signal / control edit rate",
            Criterion::MinCoverage => "signal + control reads",
            Criterion::MinConverted => "converted signal reads",
            Criterion::MinEditRatio | Criterion::MaxEditRatio => "converted / coverage",
            Criterion::MinCells => "kept cells with a converted read",
        }
    }

    /// Whether the site keeps the high end of the displayed axis (all but
    /// the upper edit-ratio bound: a p-value is shown as -log10 p).
    fn keeps_high(self) -> bool {
        self != Criterion::MaxEditRatio
    }

    fn is_integer(self) -> bool {
        matches!(
            self,
            Criterion::MinCoverage | Criterion::MinConverted | Criterion::MinCells
        )
    }

    fn scales(self) -> &'static [Scale] {
        match self {
            Criterion::MinLogOdds => &[Scale::Linear],
            // Log bins resolve the p ~ 0.05 region the decision is made in.
            Criterion::MaxPv => &[Scale::Log, Scale::Sqrt, Scale::Linear],
            Criterion::MinEditRatio | Criterion::MaxEditRatio => &[Scale::Linear, Scale::Sqrt],
            _ => &[Scale::Log, Scale::Sqrt, Scale::Linear],
        }
    }

    /// Position of a raw value on the displayed axis (may be infinite).
    fn display(self, raw: f64) -> f64 {
        match self {
            Criterion::MaxPv if raw <= 0.0 => f64::INFINITY,
            Criterion::MaxPv => (-raw.log10()).max(0.0),
            _ => raw,
        }
    }

    /// Exact, for the flags.
    fn fmt(self, v: f64) -> String {
        if self.is_integer() {
            format!("{}", v as u64)
        } else {
            format!("{}", v as f32)
        }
    }

    /// Three significant digits, for the screen.
    fn fmt_short(self, v: f64) -> String {
        let a = v.abs();
        if self.is_integer() || !v.is_finite() || a == 0.0 {
            self.fmt(v)
        } else if !(1e-3..1e5).contains(&a) {
            format!("{v:.2e}")
        } else {
            let decimals = (2 - a.log10().floor() as i32).max(0) as usize;
            let s = format!("{v:.decimals$}");
            if s.contains('.') {
                s.trim_end_matches('0').trim_end_matches('.').to_string()
            } else {
                s
            }
        }
    }

    fn bit(self) -> u8 {
        1 << self as u8
    }
}

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
    /// Sites each criterion drops alone (the others off), by `c as usize`.
    alone: [usize; 8],
    /// Sites each criterion drops first, by `c as usize`.
    first: [usize; 8],
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
                t.first[m.trailing_zeros() as usize] += 1;
                for (b, n) in t.alone.iter_mut().enumerate() {
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

/// The focused criterion's column, binned.
struct Column {
    /// Display value per site, infinities clamped to the finite range.
    display: Vec<f32>,
    n_inf: usize,
    lo: f64,
    hi: f64,
    hist: Binned,
    /// Histogram slot per site.
    slot: Vec<u16>,
    /// Thresholds the bin steps visit, ascending on the displayed axis:
    /// `(bin, raw)` with `raw` the loosest threshold that keeps the whole bin
    /// and drops the bins on the far side. Stops that keep everything are
    /// left out: they read as off.
    stops: Vec<(i32, f64)>,
}

impl Column {
    fn new(view: &SiteView, c: Criterion, scale: Scale) -> Self {
        let mut display: Vec<f32> = (0..view.table.len())
            .map(|i| c.display(c.raw(view.table, i, view.cells(i))) as f32)
            .collect();
        let (mut lo, mut hi, mut n_inf) = (f64::INFINITY, f64::NEG_INFINITY, 0);
        for &d in &display {
            if d.is_finite() {
                lo = lo.min(d as f64);
                hi = hi.max(d as f64);
            } else {
                n_inf += 1;
            }
        }
        if lo > hi {
            (lo, hi) = (0.0, 0.0);
        }
        display
            .iter_mut()
            .for_each(|d| *d = d.clamp(lo as f32, hi as f32));
        let (hist, slot, stops) = bin(view, c, scale, &display, lo, hi);
        Self {
            display,
            n_inf,
            lo,
            hi,
            hist,
            slot,
            stops,
        }
    }

    fn rebin(&mut self, view: &SiteView, c: Criterion, scale: Scale) {
        (self.hist, self.slot, self.stops) = bin(view, c, scale, &self.display, self.lo, self.hi);
    }

    /// Bin key of a raw threshold, clamped to the histogram.
    fn key_of(&self, c: Criterion, raw: f64) -> i32 {
        let d = c.display(raw).clamp(self.lo, self.hi);
        let kmax = self.hist.kmax();
        self.hist.bins.key(d).clamp(self.hist.kmin, kmax)
    }
}

/// Bin `display` on `scale`: the histogram, each site's slot, and the stops.
fn bin(
    view: &SiteView,
    c: Criterion,
    scale: Scale,
    display: &[f32],
    lo: f64,
    hi: f64,
) -> (Binned, Vec<u16>, Vec<(i32, f64)>) {
    // A signed axis spans -max..max; widen so it still gets ~50 bins.
    let max = if lo < 0.0 {
        2.0 * lo.abs().max(hi.abs())
    } else {
        hi
    };
    let bins = Binning::new(scale, max, c.is_integer());
    let kmin = bins.key(lo);
    let n = (bins.key(hi) - kmin + 1).clamp(1, MAX_BINS) as usize;
    let mut counts = vec![0; n];
    let mut slot = Vec::with_capacity(display.len());
    let mut edge: Vec<Option<f64>> = vec![None; n];
    for (i, &d) in display.iter().enumerate() {
        let b = (bins.key(d as f64) - kmin).clamp(0, n as i32 - 1) as usize;
        counts[b] += 1;
        slot.push(b as u16);
        let raw = c.raw(view.table, i, view.cells(i));
        edge[b] = Some(match edge[b] {
            None => raw,
            Some(e) if c.is_max() => e.max(raw),
            Some(e) => e.min(raw),
        });
    }
    let stops = edge
        .iter()
        .enumerate()
        .filter_map(|(b, e)| e.map(|raw| (kmin + b as i32, raw)))
        .filter(|&(_, raw)| {
            let mut f = SiteFilterArgs::permissive();
            c.set(&mut f, raw);
            !c.is_off(&f)
        })
        .collect();
    (Binned { bins, kmin, counts }, slot, stops)
}

enum Mode {
    Browse,
    /// Typing an exact raw threshold for the focused criterion.
    Edit(String),
}

/// What Enter does.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Purpose {
    /// `faba qc`: apply the thresholds.
    Apply,
    /// `faba qc-report`: print the matching `faba qc` flags.
    Explore,
}

/// State of the picker, independent of the terminal so it can be tested.
struct SitePicker<'a> {
    title: String,
    views: Vec<SiteView<'a>>,
    purpose: Purpose,
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
    decision: Option<Option<SiteFilterArgs>>,
}

impl<'a> SitePicker<'a> {
    fn new(
        title: &str,
        views: Vec<SiteView<'a>>,
        filter: SiteFilterArgs,
        purpose: Purpose,
    ) -> Self {
        assert!(!views.is_empty(), "no site table to pick thresholds on");
        let c = views[0].criteria[0];
        let column = Column::new(&views[0], c, c.scales()[0]);
        let tally = Tally::new(&views[0], c, &column);
        Self {
            title: title.to_string(),
            views,
            purpose,
            initial: filter.clone(),
            filter,
            modality: 0,
            focus: 0,
            x_scale: [0; 8],
            y_scale: Scale::Log,
            column,
            tally,
            mode: Mode::Browse,
            decision: None,
        }
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

    fn criteria_lines(&self) -> Vec<Line<'static>> {
        let dim = |t: String| Span::styled(t, DIM);
        let mut lines = vec![Line::from(dim(format!(
            "  {:<13}{:>11}{:>9}{:>9}",
            "", "threshold", "alone", "first"
        )))];
        for (j, &c) in self.view().criteria.iter().enumerate() {
            let focused = j == self.focus;
            let (value, value_style) = if c.is_off(&self.filter) {
                ("off".to_string(), DIM)
            } else {
                let op = if c.is_max() { "≤" } else { "≥" };
                let value = format!("{op} {}", c.fmt_short(c.get(&self.filter)));
                (value, if focused { HIGHLIGHT } else { PLAIN })
            };
            lines.push(Line::from(vec![
                Span::styled(if focused { "▸ " } else { "  " }, HIGHLIGHT),
                Span::styled(
                    format!("{:<13}", c.label()),
                    if focused { HIGHLIGHT } else { PLAIN },
                ),
                Span::styled(format!("{value:>11}"), value_style),
                Span::styled(format!("{:>9}", self.tally.alone[c as usize]), ACCENTED),
                Span::styled(format!("{:>9}", self.tally.first[c as usize]), DIM),
            ]));
        }
        let n = self.view().table.len();
        lines.push(Line::from(""));
        lines.push(Line::from(vec![
            dim("  kept ".into()),
            Span::styled(format!("{}", self.tally.kept), HIGHLIGHT),
            dim(format!(" / {n} sites ({:.2}%)", pct(self.tally.kept, n))),
        ]));
        lines.push(Line::from(vec![
            dim("  in ".into()),
            Span::raw(format!("{}", self.tally.genes)),
            dim(format!(" / {} genes", self.view().table.n_genes)),
        ]));
        lines.push(Line::from(dim("  alone: with the others off".into())));
        lines.push(Line::from(dim("  first: the first check it fails".into())));
        lines
    }

    fn render_hist(&self, frame: &mut Frame, area: Rect) {
        let c = self.criterion();
        let block = panel(format!(" {} ", c.axis()), true);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let [stats, plot] =
            Layout::vertical([Constraint::Length(2), Constraint::Min(5)]).areas(inner);

        let col = &self.column;
        let dim = |t: String| Span::styled(t, DIM);
        let mut first = vec![
            dim("min ".into()),
            Span::raw(c.fmt_short(col.lo)),
            dim("   max ".into()),
            Span::raw(c.fmt_short(col.hi)),
        ];
        if col.n_inf > 0 {
            first.push(dim(format!(
                "   {} infinite, drawn at the end bin",
                col.n_inf
            )));
        }
        let front: usize = self.tally.subset.iter().sum();
        let second = vec![
            dim("front: ".into()),
            Span::raw(front.to_string()),
            dim(" sites passing every other threshold   ".into()),
            Span::styled("accent", ACCENTED),
            dim(" = dropped".into()),
        ];
        frame.render_widget(
            Paragraph::new(vec![Line::from(first), Line::from(second)]),
            stats,
        );

        let pointer = (!c.is_off(&self.filter)).then(|| col.key_of(c, c.get(&self.filter)));
        let keeps_high = c.keeps_high();
        let style = |k: i32| match pointer {
            Some(p) if (keeps_high && k < p) || (!keeps_high && k > p) => ACCENTED,
            _ => PLAIN,
        };
        HistPlot {
            bins: col.hist.bins,
            kmin: col.hist.kmin,
            counts: &col.hist.counts,
            style: &style,
            subset: Some(&self.tally.subset),
            y_scale: self.y_scale,
            pointer,
            marks: Vec::new(),
        }
        .render(frame.buffer_mut(), plot);
    }
}

impl Screen for SitePicker<'_> {
    fn done(&self) -> bool {
        self.decision.is_some()
    }

    fn interrupt(&mut self) {
        self.decision = Some(None);
    }

    fn handle_key(&mut self, key: KeyEvent) {
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
            Mode::Browse => match key.code {
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
                KeyCode::Enter => self.decision = Some(Some(self.filter.clone())),
                KeyCode::Char('q') | KeyCode::Esc => self.decision = Some(None),
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

        let scales = format!("x {} · y {}", self.scale().name(), self.y_scale.name());
        frame.render_widget(header("qc", &self.title, &scales), top);

        let mut spans = vec![Span::raw(" ")];
        for (i, v) in self.views.iter().enumerate() {
            let style = if i == self.modality { HIGHLIGHT } else { DIM };
            spans.push(Span::styled(format!("{}  ", v.table.modality), style));
        }
        frame.render_widget(Line::from(spans), tabs);

        let [left, right] =
            Layout::horizontal([Constraint::Length(48), Constraint::Fill(1)]).areas(body);
        let block = panel(" site thresholds ".into(), false);
        let inner = block.inner(left);
        frame.render_widget(block, left);
        frame.render_widget(Paragraph::new(self.criteria_lines()), inner);
        self.render_hist(frame, right);

        let enter = match self.purpose {
            Purpose::Apply => "apply",
            Purpose::Explore => "print flags",
        };
        let help = match &self.mode {
            Mode::Edit(buf) => input_line(
                &format!("{} {}: ", self.criterion().label(), self.criterion().flag()),
                buf,
                &[("Enter", "set"), ("Esc", "back")],
            ),
            Mode::Browse => help_line(&[
                ("↑/↓", "knob"),
                ("←/→", "bin"),
                ("-/+", "±1"),
                ("0-9", "type"),
                ("o", "off"),
                ("r", "reset"),
                ("x/y", "scale"),
                ("Tab", "modality"),
                ("Enter", enter),
                ("q", "cancel"),
            ]),
        };
        frame.render_widget(help, footer);
    }
}

/// How a picker session ended.
pub enum Picked {
    /// No terminal, or no site table: the thresholds stand as given.
    Skipped,
    Cancelled,
    Chosen(SiteFilterArgs),
}

/// Open the picker over the non-empty site tables, in [`SITE_MODALITIES`]
/// order, with each modality's kept cells per site where `_site` matrices
/// gave them.
pub fn run_site_picker(
    input_dir: &str,
    tables: &FxHashMap<Box<str>, SiteTable>,
    site_cells: &FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>,
    filter: &SiteFilterArgs,
    purpose: Purpose,
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
    let mut picker = SitePicker::new(&file_name(input_dir), views, filter.clone(), purpose);
    run_screen(&mut picker)?;
    Ok(match picker.decision.flatten() {
        Some(f) => Picked::Chosen(f),
        None => Picked::Cancelled,
    })
}

#[cfg(test)]
#[path = "tests/site_tui.rs"]
mod tests;
