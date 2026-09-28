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

/// How each criterion is drawn.
impl Criterion {
    fn label(self) -> &'static str {
        match self {
            Criterion::MaxPv => "p-value",
            Criterion::MinLogOdds => "log odds",
            Criterion::MinFold => "fold",
            Criterion::MinCoverage => "coverage",
            Criterion::MinConverted => "converted",
            Criterion::MinEditRatio => "min edit ratio",
            Criterion::MaxEditRatio => "max edit ratio",
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
            Criterion::MinEditRatio => "converted / coverage: low is weak editing",
            Criterion::MaxEditRatio => "converted / coverage: high is variant-like",
            Criterion::MinCells => "kept cells with a converted read",
        }
    }

    /// Whether the site keeps the high end of the displayed axis (all but
    /// the upper edit-ratio bound: a p-value is shown as -log10 p).
    fn keeps_high(self) -> bool {
        self != Criterion::MaxEditRatio
    }

    /// Whether bin `k` lies on the dropped side of the threshold's bin.
    fn drops(self, k: i32, threshold: Option<i32>) -> bool {
        match threshold {
            Some(p) if self.keeps_high() => k < p,
            Some(p) => k > p,
            None => false,
        }
    }

    /// The threshold as shown: `≥ v` or `≤ v`, `None` when off.
    fn shown(self, f: &SiteFilterArgs) -> Option<String> {
        let op = if self.is_max() { "≤" } else { "≥" };
        (!self.is_off(f)).then(|| format!("{op} {}", self.fmt_short(self.get(f))))
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
    /// The smallest display value in each bin, for tick labels.
    lowest: Vec<f64>,
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
        let (hist, slot, stops, lowest) = bin(view, c, scale, &display, lo, hi);
        Self {
            display,
            n_inf,
            lo,
            hi,
            hist,
            slot,
            stops,
            lowest,
        }
    }

    fn rebin(&mut self, view: &SiteView, c: Criterion, scale: Scale) {
        (self.hist, self.slot, self.stops, self.lowest) =
            bin(view, c, scale, &self.display, self.lo, self.hi);
    }

    /// Bin key of a raw threshold, clamped to the histogram.
    fn key_of(&self, c: Criterion, raw: f64) -> i32 {
        let d = c.display(raw).clamp(self.lo, self.hi);
        let kmax = self.hist.kmax();
        self.hist.bins.key(d).clamp(self.hist.kmin, kmax)
    }
}

/// A column binned: the histogram, each site's slot, the stops, and each
/// bin's smallest value.
type ColumnBins = (Binned, Vec<u16>, Vec<(i32, f64)>, Vec<f64>);

/// Bin `display` on `scale`: the histogram, each site's slot, and the stops.
fn bin(
    view: &SiteView,
    c: Criterion,
    scale: Scale,
    display: &[f32],
    lo: f64,
    hi: f64,
) -> ColumnBins {
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
    let mut lowest = vec![f64::INFINITY; n];
    for (i, &d) in display.iter().enumerate() {
        let b = (bins.key(d as f64) - kmin).clamp(0, n as i32 - 1) as usize;
        counts[b] += 1;
        lowest[b] = lowest[b].min(d as f64);
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
    (Binned { bins, kmin, counts }, slot, stops, lowest)
}

/// Total metagene bins across 5'UTR, CDS and 3'UTR: enough for the shape,
/// and at two columns a bin they fill a panel on a 150-column terminal.
const META_BINS: usize = 48;

/// Each coding region's name at its middle bin, where a short UTR's label
/// does not run into the next one.
fn meta_ticks(regions: [usize; 3]) -> Vec<(usize, String)> {
    let mut at = 0;
    let mut out = Vec::new();
    for (r, &n) in regions.iter().enumerate() {
        if n > 0 {
            out.push((at + n / 2, REGION_NAMES[r].to_string()));
        }
        at += n;
    }
    out
}

/// The first bin of the CDS and of the 3'UTR: the region boundaries.
fn meta_marks(regions: [usize; 3]) -> Vec<usize> {
    vec![regions[0], regions[0] + regions[1]]
}

/// What the gene and metagene bars add up: sites, or their converted reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Weight {
    Sites,
    Converted,
}

impl Weight {
    /// The y axis, as the titles name it.
    fn unit(self) -> &'static str {
        match self {
            Weight::Sites => "sites",
            Weight::Converted => "converted reads",
        }
    }

    fn other(self) -> Self {
        match self {
            Weight::Sites => Weight::Converted,
            Weight::Converted => Weight::Sites,
        }
    }
}

/// What the annotation gives the picker: per view, in view order, the
/// metagene layout (`None` where no site is on a coding transcript); and the
/// gene models of the genes the tables name, by `{gene_id}_{symbol}`.
struct Annotation {
    views: Vec<Option<MetaLayout>>,
    /// Why there is no metagene at all, when the annotation is to blame.
    no_metagene: Option<String>,
    models: FxHashMap<Box<str>, GeneModel>,
}

/// The annotation's state.
enum Meta {
    /// Nothing to draw; the reason, for the panels.
    Unavailable(String),
    /// The annotation being read and the sites placed, on a thread.
    Pending(std::thread::JoinHandle<Result<Annotation, String>>),
    Ready(Annotation),
}

impl Meta {
    /// Read `gff` once, on a thread: place every table's sites (`batches`,
    /// in view order) on the metagene, and keep the models of the genes in
    /// `keys`.
    fn start(gff: Option<&str>, batches: Vec<RecordBatch>, keys: FxHashSet<Box<str>>) -> Self {
        let Some(gff) = gff.map(str::to_string) else {
            return Meta::Unavailable(
                "no annotation: pass --gff, or keep the run record next to the sites".into(),
            );
        };
        Meta::Pending(std::thread::spawn(move || {
            let fail = |e: anyhow::Error| format!("{gff}: {e:#}");
            let records = read_gff_record_vec(&gff).map_err(fail)?;
            let meta = MetaModels::from_records(&records);
            let no_metagene = meta.is_empty().then(|| {
                format!("{gff}: no coding transcript (exon and CDS lines) to build a metagene")
            });
            let views = batches
                .iter()
                .map(|b| {
                    let sites = genomic_sites(b).map_err(|e| format!("{e:#}"))?;
                    Ok(meta.layout(&sites, META_BINS))
                })
                .collect::<Result<Vec<_>, String>>()?;
            let models = gene_models_from_records(&records, |k| keys.contains(k))
                .map_err(fail)?
                .into_iter()
                .map(|m| (m.key.clone(), m))
                .collect();
            Ok(Annotation {
                views,
                no_metagene,
                models,
            })
        }))
    }

    #[cfg(test)]
    fn ready(layouts: Vec<Option<MetaLayout>>, models: Vec<GeneModel>) -> Self {
        Meta::Ready(Annotation {
            views: layouts,
            no_metagene: None,
            models: models.into_iter().map(|m| (m.key.clone(), m)).collect(),
        })
    }

    /// Take the thread's result if it has arrived; true when it just did.
    fn poll(&mut self) -> bool {
        let Meta::Pending(handle) = self else {
            return false;
        };
        if !handle.is_finished() {
            return false;
        }
        let Meta::Pending(handle) = std::mem::replace(self, Meta::Unavailable(String::new()))
        else {
            unreachable!()
        };
        *self = match handle.join() {
            Ok(Ok(annotation)) => Meta::Ready(annotation),
            Ok(Err(e)) => Meta::Unavailable(e),
            Err(_) => Meta::Unavailable("reading the annotation panicked".into()),
        };
        true
    }
}

/// A colour key: a block in each bar style, then what it stands for.
fn key_line(items: &[(Style, String)]) -> Line<'static> {
    let mut spans = Vec::new();
    for (style, what) in items {
        spans.push(Span::styled("█ ", *style));
        spans.push(Span::styled(format!("{what}   "), DIM));
    }
    Line::from(spans)
}

/// The key of the gene and metagene plots, with the key that switches
/// what they add up.
fn kept_key(weight: Weight) -> Line<'static> {
    let mut line = key_line(&[
        (DIM, format!("all {}", weight.unit())),
        (PLAIN, "kept".into()),
    ]);
    line.push_span(Span::styled(format!("c: {}", weight.other().unit()), DIM));
    line
}

/// One gene's sites along its span, all and kept, with its exons.
struct GeneProfile {
    symbol: String,
    chr: Option<String>,
    lo: i64,
    hi: i64,
    forward: Option<bool>,
    all: Vec<usize>,
    kept: Vec<usize>,
    /// What `all` and `kept` add up.
    unit: &'static str,
    /// Per bin, whether it overlaps an exon; `None` without a gene model.
    exons: Option<Vec<bool>>,
}

impl GeneProfile {
    /// Base pairs per bin.
    fn bin_bp(&self) -> i64 {
        ((self.hi - self.lo) as f64 / self.all.len().max(1) as f64).round() as i64
    }

    fn title(&self) -> String {
        let strand = match self.forward {
            Some(true) => " (+)",
            Some(false) => " (-)",
            None => "",
        };
        let chr = self
            .chr
            .as_deref()
            .map_or(String::new(), |c| format!("{c}:"));
        format!("{}  {chr}{}-{}{strand}", self.symbol, self.lo, self.hi)
    }
}

/// A gene's profile in the box `(x, y, w, h)`: its sites as bars, all
/// behind and kept in front, over its exons.
fn draw_gene(canvas: &mut Canvas, p: &GeneProfile, (x, y, w, h): (f64, f64, f64, f64)) {
    let all: Vec<f64> = p.all.iter().map(|&v| v as f64).collect();
    let kept: Vec<f64> = p.kept.iter().map(|&v| v as f64).collect();
    Bars {
        values: &all,
        front: Some(&kept),
        accent: &|_| false,
        y_scale: Scale::Linear,
        y_max: None,
        ticks: Vec::new(),
        pointer: None,
        marks: Vec::new(),
        title: format!("{}: {} per {} bp", p.title(), p.unit, p.bin_bp()),
        x_title: String::new(),
        y_title: p.unit.into(),
    }
    .draw(canvas, x, y, w, h);
    // The exon track in the space Bars keeps for x labels (52 left, 12
    // right, 38 below), which a gene's bars leave empty.
    let (px, pw) = (x + 52.0, w - 64.0);
    let ty = y + h - 38.0 + 14.0;
    canvas.line(px, ty, px + pw, ty, MUTED, 0.6);
    if let Some(exons) = &p.exons {
        let bw = pw / exons.len().max(1) as f64;
        for (b, _) in exons.iter().enumerate().filter(|(_, &e)| e) {
            let x0 = px + b as f64 * bw;
            canvas.line(x0, ty, x0 + bw, ty, INK, 5.0);
        }
    }
}

/// A view's metagene: all sites, the kept ones, and the bins per region.
struct MetaCounts<'m> {
    all: &'m [usize],
    kept: &'m [usize],
    regions: [usize; 3],
    unassigned: usize,
}

enum Mode {
    Browse,
    /// Typing an exact raw threshold for the focused criterion.
    Edit(String),
    /// Typing into the gene list's filter.
    Find,
}

/// What the table's two counts mean, a line each for the screen and the
/// figure.
const LEGEND: [&str; 4] = [
    "filtered out: sites failing this threshold.",
    "  A site failing several thresholds counts in",
    "  each row, so the rows add up to more than the",
    "  sites removed.",
];

/// What the figure's bar colours mean; on screen each plot has its own key.
const FIGURE_KEY: [&str; 3] = [
    "Bars: light gray, all sites; dark gray, kept;",
    "  orange, filtered out only by the focused",
    "  threshold.",
];

/// The count column's header, a line each.
const COUNT_HEADER: (&str, &str) = ("sites", "filtered out");

/// Widths of the table's columns: marker and name, threshold, count.
const TABLE_WIDTHS: [usize; 3] = [17, 10, 14];

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

    /// List `pinned` (symbols or keys, case-insensitive) first, in the order
    /// given, then every other gene by number of putative sites.
    fn pin_genes(mut self, pinned: &[Box<str>]) -> Self {
        let rank = |key: &str, symbol: &str| {
            pinned
                .iter()
                .position(|p| p.eq_ignore_ascii_case(key) || p.eq_ignore_ascii_case(symbol))
        };
        self.gene_order = self
            .views
            .iter()
            .map(|v| {
                let Some(g) = &v.genes else {
                    return Vec::new();
                };
                let mut order: Vec<u32> = (0..g.keys.len() as u32).collect();
                order.sort_by_key(|&i| {
                    let i = i as usize;
                    let pin = rank(&g.keys[i], g.symbol(i)).unwrap_or(usize::MAX);
                    (
                        pin,
                        std::cmp::Reverse(g.rows[i].len()),
                        g.symbol(i).to_string(),
                    )
                });
                order
            })
            .collect();
        self
    }

    /// The focused view's genes that pass the filter, in list order.
    fn gene_list(&self) -> Vec<u32> {
        let Some(g) = &self.view().genes else {
            return Vec::new();
        };
        let find = self.gene_find.to_lowercase();
        self.gene_order[self.modality]
            .iter()
            .copied()
            .filter(|&i| find.is_empty() || g.symbol(i as usize).to_lowercase().contains(&find))
            .collect()
    }

    /// The selected gene's dense id.
    fn gene(&self) -> Option<usize> {
        self.gene_list().get(self.gene_at).map(|&g| g as usize)
    }

    fn step_gene(&mut self, delta: isize) {
        let n = self.gene_list().len() as isize;
        if n > 0 {
            self.gene_at = (self.gene_at as isize + delta).clamp(0, n - 1) as usize;
        }
    }

    /// Kept and all sites of gene `g` in the focused view.
    fn gene_kept(&self, g: usize) -> (usize, usize) {
        let Some(genes) = &self.view().genes else {
            return (0, 0);
        };
        let rows = &genes.rows[g];
        let fails = &self.view().fails;
        let kept = rows.iter().filter(|&&i| fails[i as usize] == 0).count();
        (kept, rows.len())
    }

    /// The selected gene's sites along its span in `n` bins: all, kept,
    /// and which bins overlap an exon (`None` without a gene model).
    fn gene_profile(&self, n: usize) -> Option<GeneProfile> {
        let g = self.gene()?;
        let genes = self.view().genes.as_ref()?;
        let rows = &genes.rows[g];
        let model = match &self.meta {
            Meta::Ready(a) => a.models.get(&genes.keys[g]),
            _ => None,
        };
        let pos = |i: &u32| genes.pos[*i as usize];
        // The model's span, widened to every site: an annotation other than
        // the one the sites were called on may not cover them all.
        let (min, max) = (rows.iter().map(pos).min()?, rows.iter().map(pos).max()? + 1);
        let (lo, hi) = model.map_or((min, max), |m| (m.lo.min(min), m.hi.max(max)));
        let n = n.max(1);
        let span = (hi - lo).max(1);
        let bin = |p: i64| (((p - lo).clamp(0, span - 1) * n as i64) / span) as usize;
        let (mut all, mut kept) = (vec![0usize; n], vec![0usize; n]);
        let fails = &self.view().fails;
        for i in rows {
            let b = bin(pos(i));
            let w = self.site_weight(*i as usize);
            all[b] += w;
            if fails[*i as usize] == 0 {
                kept[b] += w;
            }
        }
        let exons = model.map(|m| {
            (0..n as i64)
                .map(|b| {
                    let (a, z) = (lo + b * span / n as i64, lo + (b + 1) * span / n as i64);
                    m.exons.iter().any(|&(s, e)| s < z.max(a + 1) && e > a)
                })
                .collect()
        });
        Some(GeneProfile {
            symbol: genes.symbol(g).to_string(),
            chr: model.map(|m| m.chr.to_string()),
            lo,
            hi,
            forward: model.map(|m| m.forward),
            all,
            kept,
            unit: self.weight.unit(),
            exons,
        })
    }

    fn render_gene_list(&self, frame: &mut Frame, area: Rect) {
        let block = panel(" genes: kept / all sites ".into(), false);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let Some(genes) = &self.view().genes else {
            let why = "this site table has no gene column";
            frame.render_widget(Paragraph::new(Span::styled(why, DIM)), inner);
            return;
        };
        let list = self.gene_list();
        let rows = inner.height.saturating_sub(1) as usize;
        let first = self.gene_at.saturating_sub(rows.saturating_sub(1) / 2);
        let first = first.min(list.len().saturating_sub(rows));
        let width = inner.width as usize;
        let mut lines: Vec<Line> = list
            .iter()
            .enumerate()
            .skip(first)
            .take(rows)
            .map(|(j, &g)| {
                let g = g as usize;
                let (kept, all) = self.gene_kept(g);
                let counts = format!("{kept} / {all}");
                let name_w = width.saturating_sub(counts.len() + 3);
                let selected = j == self.gene_at;
                let style = if selected { HIGHLIGHT } else { PLAIN };
                Line::from(vec![
                    Span::styled(if selected { "▸ " } else { "  " }, HIGHLIGHT),
                    Span::styled(format!("{:<name_w$.name_w$}", genes.symbol(g)), style),
                    Span::styled(format!(" {counts}"), DIM),
                ])
            })
            .collect();
        if list.is_empty() {
            lines.push(Line::from(Span::styled("  no gene matches", DIM)));
        }
        let find = match self.mode {
            Mode::Find => format!("  / {}_", self.gene_find),
            _ if !self.gene_find.is_empty() => format!("  / {}", self.gene_find),
            _ => format!("  {} genes   [ ] move  / find", list.len()),
        };
        let [list_area, find_area] =
            Layout::vertical([Constraint::Fill(1), Constraint::Length(1)]).areas(inner);
        frame.render_widget(Paragraph::new(lines), list_area);
        frame.render_widget(Paragraph::new(Span::styled(find, DIM)), find_area);
    }

    fn render_gene(&self, frame: &mut Frame, area: Rect) {
        let n = area.width.saturating_sub(2 + 6).max(10) as usize;
        let profile = self.gene_profile(n);
        let title = match &profile {
            Some(p) => format!(" {} · y: {} per {} bp ", p.title(), p.unit, p.bin_bp()),
            None => " gene ".into(),
        };
        let block = panel(title, true);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let Some(p) = profile else {
            let why = "no gene selected";
            frame.render_widget(Paragraph::new(Span::styled(why, DIM)), inner);
            return;
        };
        let [plot, track, key] = Layout::vertical([
            Constraint::Min(4),
            Constraint::Length(1),
            Constraint::Length(1),
        ])
        .areas(inner);
        frame.render_widget(Paragraph::new(kept_key(self.weight)), key);
        HistPlot {
            bins: Binning::with_width(Scale::Linear, 1.0),
            kmin: 0,
            counts: &p.all,
            style: &|_| PLAIN,
            subset: Some(&p.kept),
            y_scale: Scale::Linear,
            y_max: None,
            pointer: None,
            marks: Vec::new(),
            x_label: Some(&|_| None),
            tick_every: None,
        }
        .render(frame.buffer_mut(), plot);
        let model = match &p.exons {
            Some(exons) => exons
                .iter()
                .map(|&e| if e { "▬" } else { "─" })
                .collect::<String>(),
            None => "(no gene model)".into(),
        };
        let line = Line::from(vec![
            Span::styled(format!("{:<6}", "exons"), DIM),
            Span::styled(model, DIM),
        ]);
        frame.render_widget(Paragraph::new(line), track);
    }

    fn with_meta(mut self, meta: Meta) -> Self {
        self.set_meta(meta);
        self
    }

    fn set_meta(&mut self, meta: Meta) {
        self.meta = meta;
        self.refresh_meta();
    }

    fn view_meta(&self) -> Option<&MetaLayout> {
        let Meta::Ready(a) = &self.meta else {
            return None;
        };
        a.views.get(self.modality)?.as_ref()
    }

    /// Site `i`'s bar weight in the focused view.
    fn site_weight(&self, i: usize) -> usize {
        match self.weight {
            Weight::Sites => 1,
            Weight::Converted => self.view().table.converted[i] as usize,
        }
    }

    /// Recount the focused view's metagene over all sites, after the view,
    /// the annotation or the weight changed; then the kept sites.
    fn refresh_meta(&mut self) {
        self.meta_all = self
            .view_meta()
            .map(|m| m.counts(|i| self.site_weight(i)))
            .unwrap_or_default();
        self.meta_plot.invalidate();
        self.update_meta();
    }

    /// Recount the focused view's kept sites; redraw only if they changed.
    fn update_meta(&mut self) {
        let fails = &self.views[self.modality].fails;
        let kept = self
            .view_meta()
            .map(|m| {
                m.counts(|i| {
                    if fails[i] == 0 {
                        self.site_weight(i)
                    } else {
                        0
                    }
                })
            })
            .unwrap_or_default();
        if kept != self.meta_kept {
            self.meta_kept = kept;
            self.meta_plot.invalidate();
        }
    }

    fn switch_weight(&mut self) {
        self.weight = self.weight.other();
        self.refresh_meta();
    }

    /// The focused view's metagene under the current thresholds, `None`
    /// until the layouts arrive or when no site is on a coding transcript.
    fn meta_counts(&self) -> Option<MetaCounts<'_>> {
        let m = self.view_meta()?;
        Some(MetaCounts {
            all: &self.meta_all,
            kept: &self.meta_kept,
            regions: m.region_bins(),
            unassigned: m.unassigned,
        })
    }

    /// Why there is no metagene to draw.
    fn meta_status(&self) -> String {
        match &self.meta {
            Meta::Unavailable(why) => why.clone(),
            Meta::Pending(_) => "reading gene models ...".into(),
            Meta::Ready(a) => a
                .no_metagene
                .clone()
                .unwrap_or_else(|| "no site on a coding transcript".into()),
        }
    }

    /// The metagene in the box `(x, y, w, h)`: all sites behind, the kept
    /// sites in front.
    fn draw_meta(
        &self,
        canvas: &mut Canvas,
        m: &MetaCounts,
        bbox: (f64, f64, f64, f64),
        title: String,
    ) {
        let all: Vec<f64> = m.all.iter().map(|&v| v as f64).collect();
        let kept: Vec<f64> = m.kept.iter().map(|&v| v as f64).collect();
        Bars {
            values: &all,
            front: Some(&kept),
            accent: &|_| false,
            y_scale: Scale::Linear,
            y_max: None,
            ticks: meta_ticks(m.regions),
            pointer: None,
            marks: meta_marks(m.regions),
            title,
            x_title: "metagene position (MetaPlotR scale)".into(),
            y_title: self.weight.unit().into(),
        }
        .draw(canvas, bbox.0, bbox.1, bbox.2, bbox.3);
    }

    fn render_meta(&mut self, frame: &mut Frame, area: Rect) {
        let mut image = std::mem::take(&mut self.meta_plot);
        self.draw_meta_panel(frame, area, &mut image);
        self.meta_plot = image;
    }

    /// The metagene panel, drawn as an image into `image` where the
    /// terminal shows images.
    fn draw_meta_panel(&self, frame: &mut Frame, area: Rect, image: &mut PlotImage) {
        let block = panel(
            format!(" metagene · y: {} per bin ", self.weight.unit()),
            true,
        );
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let Some(m) = self.meta_counts() else {
            frame.render_widget(
                Paragraph::new(Line::from(Span::styled(self.meta_status(), DIM))),
                inner,
            );
            return;
        };
        let [stats, plot, key] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Min(4),
            Constraint::Length(1),
        ])
        .areas(inner);
        frame.render_widget(Paragraph::new(kept_key(self.weight)), key);
        let region_sum = |c: &[usize], r: usize| {
            let start: usize = m.regions[..r].iter().sum();
            c.get(start..start + m.regions[r])
                .map_or(0, |bins| bins.iter().sum::<usize>())
        };
        let dim = |t: String| Span::styled(t, DIM);
        let mut line = vec![dim("kept / all  ".into())];
        for (r, name) in REGION_NAMES[..3].iter().enumerate() {
            line.push(dim(format!("{name} ")));
            line.push(Span::raw(format!(
                "{}/{}   ",
                region_sum(m.kept, r),
                region_sum(m.all, r)
            )));
        }
        line.push(dim(format!("off-transcript {}", m.unassigned)));
        frame.render_widget(Paragraph::new(Line::from(line)), stats);

        let drawn = self.controls.images().is_some_and(|picker| {
            image.render(frame, plot, picker, |w, h| {
                figure::svg(w, h, |c| {
                    self.draw_meta(c, &m, (0.0, 0.0, w, h), String::new())
                })
            })
        });
        if drawn {
            return;
        }
        // Stretch the bins over the whole chart, one column each: HistPlot
        // gives a bin a whole number of columns, which leaves the rest of a
        // panel empty.
        let n = m.all.len().max(1);
        let cols = (plot.width.saturating_sub(GUTTER) as usize).max(n);
        let bin_of = |x: usize| x * n / cols;
        let stretch = |v: &[usize]| (0..cols).map(|x| v[bin_of(x)]).collect::<Vec<_>>();
        let (all, kept) = (stretch(m.all), stretch(m.kept));
        let ticks: Vec<(usize, String)> = meta_ticks(m.regions)
            .into_iter()
            .map(|(b, t)| ((2 * b + 1) * cols / (2 * n), t))
            .collect();
        let label = |k: i32| {
            ticks
                .iter()
                .find(|(i, _)| *i as i32 == k)
                .map(|(_, t)| t.clone())
        };
        HistPlot {
            bins: Binning::with_width(Scale::Linear, 1.0),
            kmin: 0,
            counts: &all,
            style: &|_| PLAIN,
            subset: Some(&kept),
            y_scale: Scale::Linear,
            y_max: None,
            pointer: None,
            marks: Vec::new(),
            x_label: Some(&label),
            tick_every: Some(1),
        }
        .render(frame.buffer_mut(), plot);
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

    fn criteria_lines(&self) -> Vec<Line<'static>> {
        let dim = |t: String| Span::styled(t, DIM);
        let [name_w, value_w, count_w] = TABLE_WIDTHS;
        let header = |first: &str, second: &str, count: &str| {
            Line::from(dim(format!(
                "{first:<name_w$}{second:>value_w$} │{count:>w$}",
                w = count_w - 1
            )))
        };
        let rule = format!(
            "{}┼{}",
            "─".repeat(name_w + value_w + 1),
            "─".repeat(count_w)
        );
        let mut lines = vec![
            header("", "", COUNT_HEADER.0),
            header("  threshold", "value", COUNT_HEADER.1),
            Line::from(dim(rule)),
        ];
        for (j, &c) in self.view().criteria.iter().enumerate() {
            let focused = j == self.focus;
            let (value, value_style) = match c.shown(&self.filter) {
                None => ("off".to_string(), DIM),
                Some(v) => (v, if focused { HIGHLIGHT } else { PLAIN }),
            };
            lines.push(Line::from(vec![
                Span::styled(if focused { "▸ " } else { "  " }, HIGHLIGHT),
                Span::styled(
                    format!("{:<w$}", c.label(), w = name_w - 2),
                    if focused { HIGHLIGHT } else { PLAIN },
                ),
                Span::styled(format!("{value:>value_w$}"), value_style),
                dim(" │".into()),
                Span::styled(
                    format!("{:>w$}", self.tally.fail[c as usize], w = count_w - 1),
                    ACCENTED,
                ),
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
        lines.push(Line::from(""));
        lines.extend(LEGEND.iter().map(|l| Line::from(dim(format!("  {l}")))));
        lines
    }

    /// The current view as a figure: the threshold table beside the
    /// focused column's histogram.
    fn figure(&self) -> String {
        let c = self.criterion();
        let view = self.view();
        let n = view.table.len();
        let meta = self.meta_counts();
        let gene = self.gene_profile(80);
        let height = 340.0
            + if gene.is_some() { 220.0 } else { 0.0 }
            + if meta.is_some() { 220.0 } else { 0.0 };
        let mut canvas = Canvas::new(720.0, height);
        canvas.bold(16.0, 22.0, &self.title, 12.0, Anchor::Start, INK);
        canvas.text(
            16.0,
            38.0,
            &format!(
                "{}: kept {} of {} sites ({:.2}%) in {} of {} genes",
                view.table.modality,
                self.tally.kept,
                n,
                pct(self.tally.kept, n),
                self.tally.genes,
                view.table.n_genes
            ),
            9.0,
            Anchor::Start,
            MUTED,
        );

        let (x0, mut y) = (16.0, 60.0);
        // Right edges of the value and count columns, and the rule between.
        let cols = [x0 + 150.0, x0 + 230.0];
        let sep = cols[0] + 8.0;
        canvas.text(x0, y + 10.0, "threshold", 8.0, Anchor::Start, MUTED);
        canvas.text(cols[0], y + 10.0, "value", 8.0, Anchor::End, MUTED);
        canvas.text(cols[1], y, COUNT_HEADER.0, 8.0, Anchor::End, MUTED);
        canvas.text(cols[1], y + 10.0, COUNT_HEADER.1, 8.0, Anchor::End, MUTED);
        let top = y - 9.0;
        y += 15.0;
        canvas.line(x0, y, cols[1] + 4.0, y, MUTED, 0.6);
        let rows_end = y + 16.0 * view.criteria.len() as f64 + 5.0;
        canvas.line(sep, top, sep, rows_end, MUTED, 0.6);
        y -= 5.0;
        for &k in &view.criteria {
            y += 16.0;
            let value = k.shown(&self.filter).unwrap_or_else(|| "off".into());
            let colour = if k == c { figure::ACCENT } else { INK };
            canvas.text(x0, y, k.label(), 9.0, Anchor::Start, colour);
            canvas.text(cols[0], y, &value, 9.0, Anchor::End, colour);
            let n = self.tally.fail[k as usize].to_string();
            canvas.text(cols[1], y, &n, 9.0, Anchor::End, INK);
        }
        y += 8.0;
        for line in LEGEND.iter().chain(&FIGURE_KEY) {
            y += 10.0;
            canvas.text(x0, y, line, 7.5, Anchor::Start, MUTED);
        }

        let title = format!("{}: sites per bin", c.label());
        self.draw_hist(&mut canvas, (320.0, 50.0, 390.0, 280.0), title);
        let mut below = 340.0;
        if let Some(p) = &gene {
            draw_gene(&mut canvas, p, (320.0, below, 390.0, 210.0));
            below += 220.0;
        }
        if let Some(m) = &meta {
            let title = format!(
                "{} metagene: {} per bin",
                view.table.modality,
                self.weight.unit()
            );
            self.draw_meta(&mut canvas, m, (320.0, below, 390.0, 210.0), title);
        }
        canvas.finish()
    }

    /// The histogram bin holding the focused threshold, `None` when off.
    fn pointer_key(&self) -> Option<i32> {
        let c = self.criterion();
        (!c.is_off(&self.filter)).then(|| self.column.key_of(c, c.get(&self.filter)))
    }

    /// The focused column's histogram in the box `(x, y, w, h)`, as the
    /// export and the in-terminal image draw it.
    fn draw_hist(&self, canvas: &mut Canvas, (x, y, w, h): (f64, f64, f64, f64), title: String) {
        let c = self.criterion();
        // Histogram ticks: the smallest value in about six bins.
        let col = &self.column;
        let (nb, lowest) = (col.lowest.len(), &col.lowest);
        let every = (nb / 6).max(1);
        let ticks = (0..nb)
            .step_by(every)
            .filter_map(|b| (b..nb).find(|&j| lowest[j].is_finite()))
            .map(|b| (b, c.fmt_short(lowest[b])))
            .collect::<Vec<_>>();
        let key = self.pointer_key();
        let pointer = key.map(|k| (k - col.hist.kmin) as usize);
        let dropped = |b: usize| c.drops(col.hist.kmin + b as i32, key);
        let values: Vec<f64> = col.hist.counts.iter().map(|&v| v as f64).collect();
        let front: Vec<f64> = self.tally.subset.iter().map(|&v| v as f64).collect();
        Bars {
            values: &values,
            front: Some(&front),
            accent: &dropped,
            y_scale: self.y_scale,
            y_max: None,
            ticks,
            pointer,
            marks: Vec::new(),
            title,
            x_title: format!("{} ({} bins)", c.axis(), self.scale().name()),
            y_title: "sites".into(),
        }
        .draw(canvas, x, y, w, h);
    }

    fn render_hist(&mut self, frame: &mut Frame, area: Rect) {
        let c = self.criterion();
        let block = panel(format!(" {} · y: sites per bin ", c.axis()), true);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let front: usize = self.tally.subset.iter().sum();
        let n = self.view().table.len();
        let key_text = key_line(&[
            (DIM, format!("all sites {n}")),
            (PLAIN, format!("pass the other thresholds {front}")),
            (
                ACCENTED,
                format!("filtered out only by this {}", self.tally.only),
            ),
        ]);
        let key_rows = if key_text.width() > inner.width as usize {
            2
        } else {
            1
        };
        let [stats, plot, key] = Layout::vertical([
            Constraint::Length(1),
            Constraint::Min(5),
            Constraint::Length(key_rows),
        ])
        .areas(inner);

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
        frame.render_widget(Paragraph::new(Line::from(first)), stats);
        frame.render_widget(
            Paragraph::new(key_text).wrap(ratatui::widgets::Wrap { trim: true }),
            key,
        );

        let mut image = std::mem::take(&mut self.plot);
        let drawn = self.controls.images().is_some_and(|picker| {
            image.render(frame, plot, picker, |w, h| {
                figure::svg(w, h, |c| self.draw_hist(c, (0.0, 0.0, w, h), String::new()))
            })
        });
        self.plot = image;
        if drawn {
            return;
        }
        let col = &self.column;
        let pointer = self.pointer_key();
        let style = |k: i32| if c.drops(k, pointer) { ACCENTED } else { PLAIN };
        HistPlot {
            bins: col.hist.bins,
            kmin: col.hist.kmin,
            counts: &col.hist.counts,
            style: &style,
            subset: Some(&self.tally.subset),
            y_scale: self.y_scale,
            y_max: None,
            pointer,
            marks: Vec::new(),
            x_label: None,
            tick_every: None,
        }
        .render(frame.buffer_mut(), plot);
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
