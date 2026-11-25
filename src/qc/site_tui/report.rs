//! `Ctrl+E` in the site view: its thresholds as one page that recounts them
//! live (see [`crate::report`]). The data is what the page recounts from:
//! each applicable criterion's raw value per site, rounded as a threshold
//! set at it rounds, so every verdict the page reaches is the view's.

use super::*;
use crate::qc::sites::OwnedColumns;
use crate::report::thousandth;
use serde::{Serialize, Serializer};
use std::sync::Arc;

/// The page's own markup, styles and code.
pub(super) const FRAGMENT: &str = include_str!("report.html");

/// The file `Ctrl+E` writes, numbered before the extension when it exists.
const FILE: (&str, &str) = ("qc_report", ".html");

/// What an infinity is written as: past f32's range, so `Math.fround`
/// reads it back as one.
const INFINITE: f64 = 1e39;

/// The most sites a page carries per modality: the strongest by p-value
/// (ties at the last included). Every site left off has a larger p-value,
/// so with the p-value threshold held at or below the last one carried
/// (the page's limit), each one left off fails it and every count stays
/// exact. The view and `faba mcp` read every site.
pub(crate) const CARRIED: usize = 50_000;

/// What the page reads from `#data`.
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub(super) struct Data {
    title: String,
    /// Every criterion some modality applies, in the view's order.
    criteria: Vec<Knob>,
    modalities: Vec<Modality>,
    /// The modality, criterion (by id) and gene on screen when written.
    modality: usize,
    focus: usize,
    gene: Option<usize>,
    /// What the gene bars add up, and the histogram's y scale.
    unit: &'static str,
    y_scale: &'static str,
    legend: [&'static str; 4],
    /// The metagene's regions, by name and colour.
    regions: [&'static str; 3],
    colours: [&'static str; 3],
    /// The 2D panels' axes, by key.
    axes: Vec<AxisInfo>,
    /// Each batch's barcodes and the cell QC's verdict on each.
    cells: Vec<CellsData>,
}

/// An axis of the 2D panels: its key, names, and the criteria it carries.
#[derive(Serialize)]
struct AxisInfo {
    key: Axis,
    name: &'static str,
    criteria: Vec<usize>,
}

/// An axis' ruler as the page reads it (see `crate::qc::scale`): its
/// range along its line, whether it places a shown value there at
/// `ln(1 + v)` (else at `v`), and its ticks there.
#[derive(Serialize)]
struct RulerInfo {
    range: (f64, f64),
    log: bool,
    ticks: Vec<(f64, String)>,
}

/// A panel's axes, and why it is shown.
#[derive(Serialize)]
struct PairInfo {
    x: Axis,
    y: Axis,
    why: String,
}

/// A batch's cell QC as the page draws it, the size of a figure whatever
/// the barcodes: the barcodes with each verdict (kept, MAD outlier, near
/// empty), the knee as its thinned points (log10 rank, log10 counts,
/// verdict) and the rank axis' end, and the density's grid.
#[derive(Serialize)]
pub(super) struct CellsData {
    batch: String,
    tally: [usize; 3],
    knee: Vec<(f32, f32, u8)>,
    knee_end: f64,
    density: crate::qc::cells_plot::Grid,
}

impl CellsData {
    pub(super) fn of(c: &crate::qc::cells_plot::Cells) -> Self {
        let (points, knee_end) = c.knee_points(crate::qc::cells_plot::KNEE_STEPS);
        Self {
            batch: c.batch.to_string(),
            tally: c.tally(),
            knee: (points.iter())
                .map(|&(x, y, v)| (x as f32, y as f32, v as u8))
                .collect(),
            knee_end,
            density: c.density_grid(),
        }
    }
}

/// A criterion as the view shows and steps it.
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct Knob {
    /// `c as usize`: its bit in a site's failed checks.
    id: usize,
    label: &'static str,
    axis: &'static str,
    flag: &'static str,
    /// An upper bound on the raw value.
    max: bool,
    keeps_high: bool,
    integer: bool,
    /// Drawn as `-log10` of the raw value.
    neg_log: bool,
    /// The x scales `x` cycles through, and the one chosen.
    scales: Vec<&'static str>,
    scale: usize,
    /// The threshold as written, and where it keeps everything.
    now: Value,
    off: Value,
}

/// One modality's sites, as columns aligned to its table.
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
struct Modality {
    name: String,
    /// The sites carried (see [`CARRIED`]), and all the table's.
    sites: usize,
    total: usize,
    /// The largest p-value carried when not every site is: the p-value
    /// threshold the page holds at or below it.
    limit: Option<f32>,
    /// The 2D panels as the view had them, the pairs Auto picks at the
    /// thresholds as written, and each axis' ruler.
    pairs: Vec<PairInfo>,
    auto: Vec<PairInfo>,
    rulers: Vec<(Axis, RulerInfo)>,
    n_genes: usize,
    /// The criteria that apply, by id in the view's order, and for each
    /// its raw value per site.
    criteria: Vec<usize>,
    raw: Vec<Vec<Value>>,
    /// Per criterion, the sites whose raw value to the nearest f32, where
    /// the view bins it, is one f32 on from `raw` (up for a lower bound,
    /// down for an upper).
    near: Vec<Vec<u32>>,
    /// Dense gene id per site.
    gene: Vec<u32>,
    /// `primary_pos` per site and the genes by dense id, when the table
    /// has the columns.
    pos: Option<Vec<i64>>,
    genes: Option<Vec<Gene>>,
    /// The gene bars' weight per site; none when they count sites.
    weight: Option<Vec<u64>>,
    /// The metagene at the thresholds as written, once there is one.
    meta: Option<MetaBars>,
}

/// A gene: its symbol, its span `[lo, hi)` as the view bins it (see
/// [`gene_span`]; none for a gene with neither sites nor a model), and its
/// model's place and exons.
#[derive(Serialize)]
struct Gene {
    name: String,
    #[serde(skip_serializing_if = "Option::is_none")]
    lo: Option<i64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    hi: Option<i64>,
    #[serde(skip_serializing_if = "Option::is_none")]
    chr: Option<String>,
    #[serde(skip_serializing_if = "Option::is_none")]
    forward: Option<bool>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    exons: Vec<(i64, i64)>,
    /// Its place among the pinned genes, which list first.
    #[serde(skip_serializing_if = "Option::is_none")]
    pin: Option<usize>,
}

/// The metagene's bars, all sites and kept, per 5'UTR, CDS and 3'UTR bin.
#[derive(Serialize)]
struct MetaBars {
    all: Vec<f64>,
    kept: Vec<f64>,
    regions: [usize; 3],
    unassigned: usize,
}

/// A value as JSON carries it: a whole count as it is; a real as its
/// shortest f32 decimal, an infinity as [`INFINITE`], NaN as null.
enum Value {
    Count(u64),
    Real(f32),
}

impl Value {
    /// `raw` for criterion `c`, rounded as setting a threshold at it rounds
    /// (a lower bound down, an upper bound up): against any f32 threshold
    /// it fails exactly when `raw` does.
    fn of(c: Criterion, raw: f64) -> Self {
        let mut f = SiteFilterArgs::permissive();
        c.set(&mut f, raw);
        let v = c.get(&f);
        if c.is_integer() {
            Value::Count(v as u64)
        } else {
            Value::Real(v as f32)
        }
    }
}

impl Serialize for Value {
    fn serialize<S: Serializer>(&self, s: S) -> Result<S::Ok, S::Error> {
        match *self {
            Value::Count(n) => s.serialize_u64(n),
            Value::Real(x) if x.is_nan() => s.serialize_none(),
            Value::Real(x) if x.is_infinite() => s.serialize_f64(INFINITE.copysign(x as f64)),
            Value::Real(x) => s.serialize_f32(x),
        }
    }
}

/// What the report's thread builds the page from, taken from the picker on
/// the UI thread by copying columns and sharing what the picker shares, so
/// the work per site (each raw value and weight and its rounding, the
/// spans, the metagene) is the thread's.
pub(super) struct Source {
    /// The page's data but its modalities, which [`Source::data`] fills.
    data: Data,
    views: Vec<ViewSource>,
    /// The annotation, once read.
    annotation: Option<Arc<Annotation>>,
    /// Each batch's cell QC, summarised on the report's thread.
    cells: super::CellsShared,
}

/// One view's columns, as [`Source`] takes them.
struct ViewSource {
    name: String,
    sites: usize,
    n_genes: usize,
    criteria: Vec<Criterion>,
    /// What the raw values and weights are read from, and the kept cells
    /// per site.
    columns: OwnedColumns,
    n_cells: Option<Vec<usize>>,
    gene: Vec<u32>,
    genes: Option<Arc<GeneSites>>,
    /// Each gene's place among the pinned genes, by dense id.
    pins: Vec<Option<usize>>,
    /// What the bars add up.
    weight: Weight,
    fails: Vec<u8>,
    /// The 2D panels, Auto's at the thresholds now, and the axes' rulers.
    pairs: Vec<(Axis, Axis, String)>,
    auto: Vec<(Axis, Axis, String)>,
    rulers: Vec<(Axis, RulerInfo)>,
}

impl Source {
    /// The page's data: every modality, per site.
    pub(super) fn data(self) -> Data {
        let Source {
            mut data,
            views,
            annotation,
            cells,
        } = self;
        data.cells = cells.iter().map(CellsData::of).collect();
        let a = annotation.as_deref();
        data.modalities = views
            .into_iter()
            .enumerate()
            .map(|(m, v)| {
                let layout = a.and_then(|a| a.views.get(m)?.as_ref());
                v.modality(a.map(|a| &a.models), layout)
            })
            .collect();
        data
    }
}

impl ViewSource {
    /// Its columns as the page reads them, its genes on `models` and its
    /// metagene on `layout` where there are some.
    fn modality(
        self,
        models: Option<&FxHashMap<Box<str>, GeneModel>>,
        layout: Option<&MetaLayout>,
    ) -> Modality {
        let columns = self.columns.as_raw();
        let cells = |i: usize| self.n_cells.as_ref().map(|n| n[i]);
        let (carried, limit) = self.carried(columns);
        let (raw, near) = self
            .criteria
            .iter()
            .map(|&c| {
                let mut near = Vec::new();
                let raw = (carried.iter().enumerate())
                    .map(|(j, &i)| {
                        let x = c.raw_in(columns, i, cells(i));
                        let value = Value::of(c, x);
                        if matches!(value, Value::Real(y) if y != x as f32 && !y.is_nan()) {
                            near.push(j as u32);
                        }
                        value
                    })
                    .collect();
                (raw, near)
            })
            .unzip();
        // Every site's weight once: the metagene's, and the carried sites'.
        let all = (self.weight != Weight::Sites).then(|| {
            (0..self.sites)
                .map(|i| self.weight.of(columns, i) as u64)
                .collect::<Vec<_>>()
        });
        let weight = (all.as_ref()).map(|w| carried.iter().map(|&i| w[i]).collect::<Vec<_>>());
        let info = |p: &[(Axis, Axis, String)]| -> Vec<PairInfo> {
            (p.iter())
                .map(|(x, y, why)| PairInfo {
                    x: *x,
                    y: *y,
                    why: why.clone(),
                })
                .collect()
        };
        // The metagene is drawn here, of every site, not only those carried.
        let meta = layout.map(|layout| self.meta(layout, all.as_deref()));
        let genes = self.genes.as_deref();
        Modality {
            name: self.name,
            sites: carried.len(),
            total: self.sites,
            limit,
            pairs: info(&self.pairs),
            auto: info(&self.auto),
            rulers: self.rulers,
            n_genes: self.n_genes,
            criteria: self.criteria.iter().map(|&c| c as usize).collect(),
            raw,
            near,
            pos: genes.map(|g| carried.iter().map(|&i| g.pos.value(i)).collect()),
            genes: genes.map(|g| genes_of(g, &self.pins, models)),
            gene: carried.iter().map(|&i| self.gene[i]).collect(),
            weight,
            meta,
        }
    }

    /// The sites the page carries, in table order, and the largest p-value
    /// among them when not every site is (see [`CARRIED`]).
    /// A site without a p-value (NaN) is always carried: it passes any
    /// p-value threshold, so only sites with a larger one are left off.
    fn carried(&self, columns: RawColumns) -> (Vec<usize>, Option<f32>) {
        let mut by_p: Vec<f32> = columns.pv.iter().copied().filter(|p| !p.is_nan()).collect();
        if by_p.len() <= CARRIED || !self.criteria.contains(&Criterion::MaxPv) {
            return ((0..self.sites).collect(), None);
        }
        let (_, limit, _) = by_p.select_nth_unstable_by(CARRIED - 1, f32::total_cmp);
        let limit = *limit;
        let carried = (0..self.sites)
            .filter(|&i| columns.pv[i].is_nan() || columns.pv[i] <= limit)
            .collect();
        (carried, Some(limit))
    }

    /// Its metagene on `layout` under the thresholds as written, each site
    /// as `weight` weighs it (one when none).
    fn meta(&self, layout: &MetaLayout, weight: Option<&[u64]>) -> MetaBars {
        let weight = |i: usize| weight.map_or(1.0, |w| w[i] as f64);
        let bars = |kept_only: bool| {
            layout
                .counts(|i| {
                    if kept_only && self.fails[i] != 0 {
                        0.0
                    } else {
                        weight(i)
                    }
                })
                .into_iter()
                .map(thousandth)
                .collect()
        };
        MetaBars {
            all: bars(false),
            kept: bars(true),
            regions: layout.region_bins(),
            unassigned: layout.unassigned,
        }
    }
}

/// `genes` by dense id, each placed among the pinned as `pins` says, on its
/// model in `models` where it has one.
fn genes_of(
    genes: &GeneSites,
    pins: &[Option<usize>],
    models: Option<&FxHashMap<Box<str>, GeneModel>>,
) -> Vec<Gene> {
    (0..genes.keys.len())
        .map(|g| {
            let model = models.and_then(|all| all.get(&genes.keys[g]));
            let pos = genes.rows[g].iter().map(|&i| genes.pos.value(i as usize));
            let span = gene_span(model, pos);
            Gene {
                name: genes.symbol(g).to_string(),
                lo: span.map(|(lo, _)| lo),
                hi: span.map(|(_, hi)| hi),
                chr: model.map(|md| md.chr.to_string()),
                forward: model.map(|md| md.forward),
                exons: model.map(|md| md.exons.clone()).unwrap_or_default(),
                pin: pins[g],
            }
        })
        .collect()
}

/// The three coding regions' entries of a per-region table.
fn coding(per_region: [&'static str; 4]) -> [&'static str; 3] {
    let [utr5, cds, utr3, _] = per_region;
    [utr5, cds, utr3]
}

impl SitePicker<'_> {
    /// What the report's thread needs, taken without work per site but
    /// copying the columns.
    pub(super) fn report_source(&self) -> Source {
        let criteria = SHOWN
            .into_iter()
            .filter(|c| self.views.iter().any(|v| v.criteria.contains(c)))
            .map(|c| Knob {
                id: c as usize,
                label: c.label(),
                axis: c.axis(),
                flag: c.flag(),
                max: c.is_max(),
                keeps_high: c.keeps_high(),
                integer: c.is_integer(),
                neg_log: c.neg_log(),
                scales: c.scales().iter().map(|s| s.name()).collect(),
                scale: self.x_scale[c as usize],
                now: Value::of(c, c.get(&self.filter)),
                off: Value::of(c, c.permissive()),
            })
            .collect();
        let data = Data {
            title: self.title.clone(),
            criteria,
            modalities: Vec::new(),
            modality: self.modality,
            focus: self.criterion() as usize,
            gene: self.gene(),
            unit: self.weight.unit(),
            y_scale: self.y_scale().name(),
            legend: LEGEND,
            regions: coding(REGION_NAMES),
            colours: coding(REGION_COLOURS),
            axes: (Axis::ALL.iter())
                .map(|&a| AxisInfo {
                    key: a,
                    name: a.name(),
                    criteria: a.criteria().iter().map(|&c| c as usize).collect(),
                })
                .collect(),
            cells: Vec::new(),
        };
        Source {
            data,
            cells: self.cells.clone(),
            views: (0..self.views.len()).map(|m| self.view_source(m)).collect(),
            annotation: match &self.meta {
                Meta::Ready(a) => Some(a.clone()),
                _ => None,
            },
        }
    }

    /// View `m`'s columns, copied, and its genes, shared.
    fn view_source(&self, m: usize) -> ViewSource {
        let v = &self.views[m];
        let t = v.table;
        // Auto's pairs once: the panels too, when none were picked.
        let auto = pairs::auto_panels(&v.sites(), pairs::PANELS);
        let axes = pairs::axes(&v.sites());
        let rulers = (axes.iter())
            .map(|&a| {
                let r = v.sites().ruler(a, &self.filter);
                let (range, log, ticks) = (r.span(), r.axis.log, r.ticks());
                (a, RulerInfo { range, log, ticks })
            })
            .collect();
        let shown = match v.pairs.is_empty() {
            true => auto.clone(),
            false => v.pairs.clone(),
        };
        ViewSource {
            name: t.modality.to_string(),
            sites: t.len(),
            n_genes: t.n_genes,
            criteria: v.criteria.clone(),
            columns: t.raw_columns().copied(),
            n_cells: v.n_cells.clone(),
            gene: t.gene_id.clone(),
            genes: v.genes.clone(),
            pins: v.genes.as_ref().map_or_else(Vec::new, |g| {
                (0..g.keys.len()).map(|g| self.list.pin_of(m, g)).collect()
            }),
            weight: self.weight,
            fails: v.fails.clone(),
            pairs: shown,
            auto,
            rulers,
        }
    }

    /// `Ctrl+E`: gather the report of the thresholds now and write it, on a
    /// thread of its own, to the first free `qc_report.html` in
    /// `report_dir`, or queue one more while one is under way: what to
    /// say. [`Self::end_report`] says where once it is written.
    pub(super) fn export(&mut self) -> String {
        if let Some(said) = self.report.busy() {
            return said;
        }
        let path = crate::tui::next_free_file(&self.report_dir, FILE.0, FILE.1);
        let (source, title, to) = (
            self.report_source(),
            format!("faba qc · {}", self.title),
            path.clone(),
        );
        self.report.spawn(path, move || {
            crate::report::write_new(&to, &title, FRAGMENT, &source.data())
        })
    }

    /// Say where the report went, or why not, once written, and write it
    /// again if that was asked meanwhile; true when it says so.
    pub(super) fn end_report(&mut self) -> bool {
        let Some(done) = self.report.poll() else {
            return false;
        };
        self.status = Some(done.then(|| self.export()));
        true
    }
}
