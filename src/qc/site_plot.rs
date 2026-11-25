//! `faba site-plot`, and the panels `faba qc --plots` writes beside its own
//! figures: one modality's site quantities as paper panels, a histogram
//! each and chosen pairs two at a time, with marginals on the pair's own
//! axes. Every figure of a quantity reads it on one [`Ruler`] (its
//! [`Scale`], range and bins), so a histogram, a 2D axis and its marginal
//! line up tick for tick. Thresholds are optional: given ones are drawn as
//! dashed lines, the sites they remove apart from those kept.
//!
//! It reads the site table and the `_site` matrices (cells per site over
//! every cell: no cell QC), and writes figures only.

use super::args::SiteFilterArgs;
use super::layout::{scan_input_dir, InputLayout, SITE_MODALITIES};
use super::matrix::{open_matrix, row_nnz_over_columns, select_columns};
use super::pairs::{self, Axis as Quantity, Cache, Sites};
use super::scale::{Ruler, Scale, Scales};
use super::site_tui::Figure;
use super::sites::{accumulate_site_cells, read_site_table, Criterion, SiteTable};
use crate::figure::style::{Sheet, StyleArgs};
use crate::figure::SizeArgs;
use clap::Args;
use data_beans::figure::panel::{DensityPanel, HistPanel, Look, Marginals};
use rustc_hash::FxHashMap;
use std::borrow::Cow;

/// A histogram's own size in points, when none is asked for.
pub const HIST_OWN: (f64, f64) = (252.0, 180.0);
/// A pair's own size in points, when none is asked for.
pub const PAIR_OWN: (f64, f64) = (252.0, 252.0);

/// Which panels, on which scales: `--hist`, `--pairs`, `--marginals`,
/// `--scale` and `--bins`, as `faba site-plot` and `faba qc --plots` take
/// them.
#[derive(Args, Debug, Clone, Default, PartialEq, serde::Serialize)]
pub struct PanelArgs {
    #[arg(
        long = "hist",
        value_name = "LIST",
        value_delimiter = ',',
        help = "Quantities to draw a histogram of: pv, fold, log_odds, coverage, converted, edit_ratio, cells (default: every one that applies; `none` for none)"
    )]
    pub hist: Vec<String>,

    #[arg(
        long = "pairs",
        value_name = "X:Y,...",
        value_delimiter = ',',
        help = "Pairs to draw two at a time, x:y, e.g. pv:cells,fold:pv (default: the pairs `faba qc` picks; `none` for none)"
    )]
    pub pairs: Vec<String>,

    #[arg(
        long = "marginals",
        value_name = "LIST",
        value_delimiter = ',',
        help = "Histograms beside each pair on its own axes: top (x), right (y), or both (default: none)"
    )]
    pub marginals: Vec<String>,

    #[arg(
        long = "scale",
        value_name = "Q=linear|log,...",
        value_delimiter = ',',
        help = "A quantity's scale, in every figure of it, e.g. pv=log,fold=linear (default: -log10 p, edit ratio and log odds linear; fold and counts log); log odds is linear only"
    )]
    pub scale: Vec<String>,

    #[arg(
        long = "bins",
        value_name = "N",
        value_parser = clap::builder::RangedU64ValueParser::<usize>::new().range(1..=1000),
        help = "Bins along each quantity's axis (default: about one per 6 points of the figure's width, 20 to 60)"
    )]
    pub bins: Option<usize>,
}

/// What was asked for, or Auto's.
#[derive(Debug, Clone, PartialEq)]
pub enum Choice<T> {
    Auto,
    Given(Vec<T>),
}

/// [`PanelArgs`] read and checked.
#[derive(Debug, Clone, PartialEq)]
pub struct Plan {
    pub hist: Choice<Quantity>,
    pub pairs: Choice<(Quantity, Quantity)>,
    pub marginals: Marginals,
    pub scales: Scales,
    pub bins: Option<usize>,
}

impl Plan {
    /// Bins along each quantity's axis on `sheet`: as asked, else about
    /// one per 6 points of a histogram's width.
    pub fn bins_on(&self, sheet: &Sheet) -> usize {
        (self.bins).unwrap_or_else(|| data_beans::figure::bins_for(sheet.fill(|| HIST_OWN).0))
    }
}

/// `items` read by `one`, Auto when there are none, none for `none`.
fn choice<T>(
    items: &[String],
    one: impl Fn(&str) -> anyhow::Result<T>,
) -> anyhow::Result<Choice<T>> {
    match items {
        [] => Ok(Choice::Auto),
        [none] if none.trim().eq_ignore_ascii_case("none") => Ok(Choice::Given(Vec::new())),
        _ => items
            .iter()
            .map(|s| one(s))
            .collect::<anyhow::Result<_>>()
            .map(Choice::Given),
    }
}

impl PanelArgs {
    /// The panels asked for, checked before anything is read.
    pub fn plan(&self) -> anyhow::Result<Plan> {
        let hist = choice(&self.hist, |s| {
            Quantity::parse(s).map_err(|e| anyhow::anyhow!("--hist: {e}"))
        })?;
        let pairs = choice(&self.pairs, |s| {
            let Some((x, y)) = s.split_once(':') else {
                anyhow::bail!("--pairs {s}: write a pair as x:y, e.g. pv:cells");
            };
            let q = |v: &str| Quantity::parse(v).map_err(|e| anyhow::anyhow!("--pairs {s}: {e}"));
            let (x, y) = (q(x)?, q(y)?);
            anyhow::ensure!(x != y, "--pairs {s}: two different quantities");
            Ok((x, y))
        })?;
        let mut marginals = Marginals::default();
        for m in &self.marginals {
            match m.trim().to_ascii_lowercase().as_str() {
                "top" => marginals.top = true,
                "right" => marginals.right = true,
                "none" => {}
                other => anyhow::bail!("--marginals {other}: top, right or both"),
            }
        }
        Ok(Plan {
            hist,
            pairs,
            marginals,
            scales: Scales::parse(&self.scale)?,
            bins: self.bins,
        })
    }
}

/// `faba qc`'s site thresholds, each optional: given ones are drawn and
/// cut the sites into kept and removed; with none, every site is kept.
#[derive(Args, Debug, Clone, Default, PartialEq, serde::Serialize)]
pub struct Thresholds {
    #[arg(
        long = "site-max-pv",
        help = "Draw and apply `faba qc`'s p-value threshold (marginal, no FDR)"
    )]
    pub site_max_pv: Option<f32>,
    #[arg(
        long = "site-min-log-odds",
        allow_hyphen_values = true,
        help = "m6A: draw and apply the raw log odds threshold"
    )]
    pub site_min_log_odds: Option<f32>,
    #[arg(
        long = "site-min-coverage",
        help = "Draw and apply the total coverage threshold"
    )]
    pub site_min_coverage: Option<u64>,
    #[arg(
        long = "site-min-converted",
        help = "Draw and apply the converted reads threshold"
    )]
    pub site_min_converted: Option<u64>,
    #[arg(
        long = "site-min-edit-ratio",
        help = "Draw and apply the lower edit ratio threshold"
    )]
    pub site_min_edit_ratio: Option<f32>,
    #[arg(
        long = "site-max-edit-ratio",
        help = "Draw and apply the upper edit ratio threshold"
    )]
    pub site_max_edit_ratio: Option<f32>,
    #[arg(
        long = "site-min-fold",
        help = "m6A: draw and apply the fold over control threshold"
    )]
    pub site_min_fold: Option<f32>,
    #[arg(
        long = "site-min-cells",
        help = "Draw and apply the cells threshold (cells with a converted read, of every cell)"
    )]
    pub site_min_cells: Option<usize>,
}

/// How sites are cut: the thresholds, those drawn and applied, and those
/// Auto ranks its pairs at.
#[derive(Debug, Clone)]
pub struct Cut {
    pub filter: SiteFilterArgs,
    pub on: Vec<Criterion>,
    pub auto: SiteFilterArgs,
}

impl Cut {
    /// Nothing cut; Auto at `faba qc`'s defaults.
    #[cfg(test)]
    pub fn none() -> Self {
        Thresholds::default().cut()
    }

    /// `faba qc`'s cut at `f`: every threshold that is not off.
    #[cfg(test)]
    pub fn of_qc(f: &SiteFilterArgs) -> Self {
        Self {
            filter: f.clone(),
            on: (Criterion::ALL.into_iter())
                .filter(|c| !c.is_off(f))
                .collect(),
            auto: f.clone(),
        }
    }
}

impl Thresholds {
    /// The thresholds given over ones that keep everything, applied; Auto
    /// ranks at them over `faba qc`'s defaults.
    pub fn cut(&self) -> Cut {
        let given = [
            (
                Criterion::MinCoverage,
                self.site_min_coverage.map(|v| v as f64),
            ),
            (
                Criterion::MinConverted,
                self.site_min_converted.map(|v| v as f64),
            ),
            (
                Criterion::MinEditRatio,
                self.site_min_edit_ratio.map(f64::from),
            ),
            (
                Criterion::MaxEditRatio,
                self.site_max_edit_ratio.map(f64::from),
            ),
            (Criterion::MinFold, self.site_min_fold.map(f64::from)),
            (Criterion::MinLogOdds, self.site_min_log_odds.map(f64::from)),
            (Criterion::MaxPv, self.site_max_pv.map(f64::from)),
            (Criterion::MinCells, self.site_min_cells.map(|v| v as f64)),
        ];
        let mut filter = SiteFilterArgs::permissive();
        let mut auto = SiteFilterArgs::default_values();
        let mut on = Vec::new();
        for (c, v) in given {
            if let Some(v) = v {
                c.set(&mut filter, v);
                c.set(&mut auto, v);
                on.push(c);
            }
        }
        Cut { filter, on, auto }
    }
}

/// One quantity as every panel of it draws it: its ruler, each site's
/// value on it (infinities at its ends), and its thresholds' edges.
struct Drawn<'a> {
    ruler: Ruler,
    values: Cow<'a, [f32]>,
    edges: Vec<f64>,
}

impl<'a> Drawn<'a> {
    /// Quantity `q` on `scale` over its shown values `shown`, spread as
    /// `spread` (see [`super::scale::spread`]), reaching the thresholds of `f` in `on`, in
    /// about `bins` bins; `shown` copied only when it has an infinity to
    /// put at an end.
    fn new(
        (q, scale): (Quantity, Scale),
        mut shown: Cow<'a, [f32]>,
        spread: Option<(f64, f64)>,
        (f, on): (&SiteFilterArgs, &[Criterion]),
        bins: usize,
    ) -> Self {
        let (ruler, edges) = Ruler::cut((q, scale), spread, f, |c| on.contains(&c), bins);
        if shown.iter().any(|v| v.is_infinite()) {
            ruler.clamp(shown.to_mut());
        }
        Self {
            ruler,
            values: shown,
            edges,
        }
    }
}

/// One modality's sites, laid out as panels.
pub struct Panels<'a> {
    modality: &'a str,
    hists: Vec<Quantity>,
    pairs: Vec<(Quantity, Quantity)>,
    marginals: Marginals,
    kept: Option<Vec<bool>>,
    drawn: [Option<Drawn<'a>>; 7],
    sheet: Sheet,
}

/// A pair of quantities drawn two at a time, `(x, y)`.
type Pair = (Quantity, Quantity);

/// The histograms and pairs `plan` asks for of `modality`'s sites, which
/// `criteria` apply to, Auto's pairs from `auto`. A quantity asked for that
/// the sites lack (fold without a control, cells without `_site` matrices)
/// is refused when `strict`, else left out with a warning.
fn chosen(
    modality: &str,
    criteria: &[Criterion],
    plan: &Plan,
    strict: bool,
    auto: impl FnOnce() -> Vec<Pair>,
) -> anyhow::Result<(Vec<Quantity>, Vec<Pair>)> {
    let applies = |q: Quantity| q.criteria().iter().any(|c| criteria.contains(c));
    let lacking = |q: Quantity| -> anyhow::Result<bool> {
        if applies(q) {
            return Ok(false);
        }
        let why = match q {
            Quantity::Cells => "no `_site` matrix to count cells from",
            _ => "no control arm",
        };
        let said = format!("{modality} sites have no {}: {why}", q.key());
        anyhow::ensure!(!strict, "{said}");
        log::warn!("{said}; not drawn");
        Ok(true)
    };
    let mut hists = Vec::new();
    for q in match &plan.hist {
        Choice::Auto => Quantity::ALL.into_iter().filter(|&q| applies(q)).collect(),
        Choice::Given(v) => v.clone(),
    } {
        if !lacking(q)? && !hists.contains(&q) {
            hists.push(q);
        }
    }
    let mut pairs = Vec::new();
    for (x, y) in match &plan.pairs {
        Choice::Auto => auto(),
        Choice::Given(v) => v.clone(),
    } {
        if !lacking(x)? && !lacking(y)? && !pairs.contains(&(x, y)) {
            pairs.push((x, y));
        }
    }
    Ok((hists, pairs))
}

/// The quantities `hists` and `pairs` draw, each once.
fn used(hists: &[Quantity], pairs: &[(Quantity, Quantity)]) -> Vec<Quantity> {
    let mut out: Vec<Quantity> = Vec::new();
    for q in (hists.iter().copied()).chain(pairs.iter().flat_map(|&(x, y)| [x, y])) {
        if !out.contains(&q) {
            out.push(q);
        }
    }
    out
}

impl<'a> Panels<'a> {
    /// The panels `plan` asks for over `t` (with its cells per site, when
    /// `_site` matrices gave them), cut by `cut`, laid out for `sheet`. A
    /// quantity asked for that the sites lack (fold without a control,
    /// cells without `_site` matrices) is refused when `strict`, else left
    /// out with a warning.
    pub fn new(
        t: &'a SiteTable,
        n_cells: Option<&[usize]>,
        cut: &Cut,
        plan: &Plan,
        sheet: Sheet,
        strict: bool,
    ) -> anyhow::Result<Self> {
        use rayon::prelude::*;
        let criteria: Vec<Criterion> = (Criterion::ALL.into_iter())
            .filter(|c| c.applies(t, n_cells.is_some()))
            .collect();
        let bins = plan.bins_on(&sheet);
        let cache = Cache::new(plan.scales, bins);
        // Auto ranks its pairs on the sites failing its thresholds.
        let fails = match plan.pairs {
            Choice::Auto => fails_at(t, n_cells, &criteria, &cut.auto),
            Choice::Given(_) => vec![0; t.len()],
        };
        let sites = Sites {
            table: t,
            n_cells,
            fails: &fails,
            criteria: &criteria,
            cache: &cache,
        };
        let (hists, pairs) = chosen(&t.modality, &criteria, plan, strict, || auto_ranked(&sites))?;
        let used = used(&hists, &pairs);
        used.par_iter().for_each(|&q| {
            sites.column(q);
        });
        let on: Vec<Criterion> = (cut.on.iter().copied())
            .filter(|c| criteria.contains(c))
            .collect();
        let kept = (!on.is_empty()).then(|| {
            (0..t.len())
                .into_par_iter()
                .map(|i| {
                    let cells = n_cells.map(|n| n[i]);
                    !on.iter().any(|c| c.fails(&cut.filter, t, i, cells))
                })
                .collect()
        });
        let mut columns = cache.into_columns();
        let mut drawn: [Option<Drawn>; 7] = Default::default();
        for q in used {
            let column = columns[q as usize]
                .take()
                .expect("a quantity drawn is read");
            drawn[q as usize] = Some(Drawn::new(
                (q, plan.scales.of(q)),
                Cow::Owned(column.shown),
                column.spread,
                (&cut.filter, &on),
                bins,
            ));
        }
        Ok(Self {
            modality: &t.modality,
            hists,
            pairs,
            marginals: plan.marginals,
            kept,
            drawn,
            sheet,
        })
    }

    /// The panels `plan` asks for over the sites a view holds, at its
    /// thresholds `f`, as [`Panels::new`] draws them at `faba qc`'s cut of
    /// `f` ([`Cut::of_qc`]), a quantity the sites lack left out: the
    /// values read off the view's columns, the sites kept off its failed
    /// checks, Auto's pairs ranked on them.
    pub fn of_view(
        s: &Sites<'a>,
        f: &SiteFilterArgs,
        plan: &Plan,
        sheet: Sheet,
    ) -> anyhow::Result<Self> {
        use rayon::prelude::*;
        let (hists, pairs) = chosen(&s.table.modality, s.criteria, plan, false, || {
            auto_ranked(s)
        })?;
        let on: Vec<Criterion> = (s.criteria.iter().copied())
            .filter(|c| !c.is_off(f))
            .collect();
        let bits = on.iter().fold(0, |b, c| b | c.bit());
        let kept = (!on.is_empty()).then(|| s.fails.par_iter().map(|&m| m & bits == 0).collect());
        let bins = plan.bins_on(&sheet);
        let mut drawn: [Option<Drawn>; 7] = Default::default();
        for q in used(&hists, &pairs) {
            let column = s.column(q);
            drawn[q as usize] = Some(Drawn::new(
                (q, plan.scales.of(q)),
                Cow::Borrowed(&column.shown),
                column.spread,
                (f, &on),
                bins,
            ));
        }
        Ok(Self {
            modality: &s.table.modality,
            hists,
            pairs,
            marginals: plan.marginals,
            kept,
            drawn,
            sheet,
        })
    }

    /// The histograms drawn, in order.
    #[cfg(test)]
    pub fn hists(&self) -> &[Quantity] {
        &self.hists
    }

    /// The pairs drawn, `(x, y)`, in order.
    #[cfg(test)]
    pub fn pairs(&self) -> &[(Quantity, Quantity)] {
        &self.pairs
    }

    fn drawn(&self, q: Quantity) -> &Drawn<'a> {
        self.drawn[q as usize]
            .as_ref()
            .expect("a quantity drawn is laid out")
    }

    /// How a figure `size` points looks: the font of the sheet, in its
    /// palette's kept and removed colours.
    fn look(&self, size: (f64, f64)) -> Look {
        let style = self.sheet.style(size);
        let roles = style.palette.roles();
        Look {
            font: style.font,
            panel_only: style.panel_only,
            kept: roles.kept.into(),
            removed: roles.removed.into(),
        }
    }

    /// The title of a figure, none with `--panel-only`.
    fn title(&self, what: String) -> String {
        match self.sheet.look.panel_only {
            true => String::new(),
            false => format!("{} {what}", self.modality),
        }
    }

    /// The histogram of `q`.
    pub fn hist(&self, q: Quantity) -> HistPanel<'_> {
        let d = self.drawn(q);
        HistPanel {
            axis: d.ruler.axis,
            bins: d.ruler.bins,
            values: &d.values,
            kept: self.kept.as_deref(),
            edges: d.edges.clone(),
            x_title: q.title().into(),
            y_title: "sites".into(),
            title: self.title(q.title().into()),
        }
    }

    /// The histogram of `q` as a figure `w` by `h` points.
    pub fn hist_figure(&self, q: Quantity, (w, h): (f64, f64)) -> String {
        self.hist(q).svg(&self.look((w, h)), w, h)
    }

    /// `y` against `x`, with the marginals asked for.
    pub fn pair(&self, x: Quantity, y: Quantity) -> DensityPanel<'_> {
        let (dx, dy) = (self.drawn(x), self.drawn(y));
        DensityPanel {
            x: dx.ruler.axis,
            x_bins: dx.ruler.bins,
            y: dy.ruler.axis,
            y_bins: dy.ruler.bins,
            xs: &dx.values,
            ys: &dy.values,
            kept: self.kept.as_deref(),
            x_edges: dx.edges.clone(),
            y_edges: dy.edges.clone(),
            x_title: x.title().into(),
            y_title: y.title().into(),
            title: self.title(format!("{} by {}", y.title(), x.title())),
            marginals: self.marginals,
        }
    }

    /// Every panel as a figure: `{modality}_hist_{q}`, then
    /// `{modality}_{y}_vs_{x}`, each the size asked (else its own).
    pub fn figures(&self) -> Vec<Figure> {
        use rayon::prelude::*;
        let (hw, hh) = self.sheet.fill(|| HIST_OWN);
        let (pw, ph) = self.sheet.fill(|| PAIR_OWN);
        let pl = self.look((pw, ph));
        let hists = (self.hists.par_iter()).map(|&q| Figure {
            stem: format!("{}_hist_{}", self.modality, q.key()),
            svg: self.hist_figure(q, (hw, hh)),
        });
        let pairs = (self.pairs.par_iter()).map(|&(x, y)| Figure {
            stem: format!("{}_{}_vs_{}", self.modality, y.key(), x.key()),
            svg: self.pair(x, y).svg(&pl, pw, ph),
        });
        hists.chain(pairs).collect()
    }
}

/// Per site of `t`, a [`Criterion::bit`] for each of `criteria` it fails
/// at `f`.
fn fails_at(
    t: &SiteTable,
    n_cells: Option<&[usize]>,
    criteria: &[Criterion],
    f: &SiteFilterArgs,
) -> Vec<u8> {
    use rayon::prelude::*;
    (0..t.len())
        .into_par_iter()
        .map(|i| {
            (criteria.iter())
                .filter(|c| c.fails(f, t, i, n_cells.map(|n| n[i])))
                .fold(0, |b, &c| b | c.bit())
        })
        .collect()
}

/// The pairs `faba qc` picks over `s`, at the thresholds its failed checks
/// were read at.
fn auto_ranked(s: &Sites) -> Vec<(Quantity, Quantity)> {
    (pairs::auto_pairs(s, pairs::PANELS).into_iter())
        .map(|p| (p.x, p.y))
        .collect()
}

/// Cells per site of `t`, aligned to it: every cell with a converted read
/// at the site, over every cell of every batch's `_site` matrix (no cell
/// QC); none without such a matrix.
pub fn site_cells(layout: &InputLayout, t: &SiteTable) -> anyhow::Result<Option<Vec<usize>>> {
    let mut acc = FxHashMap::default();
    let mut any = false;
    for m in (layout.matrices.iter()).filter(|m| m.site_modality() == Some(&*t.modality)) {
        any = true;
        let data = open_matrix(&m.path)?;
        let cols = select_columns(data.as_ref(), &data.column_names()?, None);
        let nnz = row_nnz_over_columns(data.as_ref(), &cols)?;
        accumulate_site_cells(&data.row_names()?, &nnz, t.converted_channel(), &mut acc);
    }
    Ok(any.then(|| t.cells_per_site(&acc)))
}

#[derive(Args, Debug, serde::Serialize)]
#[command(
    mut_arg("fig_width", |a| a.visible_alias("width")),
    mut_arg("fig_height", |a| a.visible_alias("height"))
)]
pub struct SitePlotArgs {
    #[arg(
        help = "A faba output directory: its {modality}_sites.parquet and {batch}_{modality}_site matrices"
    )]
    pub input_dir: Box<str>,

    #[arg(
        short = 'o',
        long = "output",
        help = "Directory the figures go to, new or empty: nothing is written over"
    )]
    pub output: Box<str>,

    #[arg(
        long = "modality",
        value_name = "m6a|atoi",
        help = "Which site table to draw (default: m6a, else atoi)"
    )]
    pub modality: Option<Box<str>>,

    #[command(flatten)]
    pub panels: PanelArgs,

    #[command(flatten)]
    pub site: Thresholds,

    #[command(flatten)]
    pub size: SizeArgs,

    #[command(flatten)]
    pub style: StyleArgs,
}

/// `--modality`, else the first site modality with a table in `layout`.
fn modality_of(
    layout: &InputLayout,
    asked: Option<&str>,
    dir: &str,
) -> anyhow::Result<&'static str> {
    let present: Vec<&'static str> = (SITE_MODALITIES.iter().copied())
        .filter(|m| layout.site_tables.contains_key(*m))
        .collect();
    match asked {
        Some(a) => {
            let m = (SITE_MODALITIES.iter().copied())
                .find(|m| m.eq_ignore_ascii_case(a))
                .ok_or_else(|| {
                    anyhow::anyhow!("--modality {a}: one of {}", SITE_MODALITIES.join(", "))
                })?;
            anyhow::ensure!(present.contains(&m), "{dir}: no {m}_sites.parquet");
            Ok(m)
        }
        None => present.first().copied().ok_or_else(|| {
            anyhow::anyhow!(
                "{dir}: no site table ({}_sites.parquet)",
                SITE_MODALITIES.join("_sites.parquet or ")
            )
        }),
    }
}

pub fn run_site_plot(args: &SitePlotArgs) -> anyhow::Result<()> {
    let plan = args.panels.plan()?;
    let sheet = Sheet::new(
        crate::figure::size::Figure::SitePanels,
        &args.size,
        &args.style,
    )?;
    let (input, out) = (&*args.input_dir, &*args.output);
    if let Some(problem) = crate::tui::output_problem(out) {
        anyhow::bail!("-o {problem}: faba site-plot never writes over files; name another");
    }
    let layout = scan_input_dir(input)?;
    let modality = modality_of(&layout, args.modality.as_deref(), input)?;
    let record = crate::run_record::RunRecord::start("site-plot", out)
        .input("fileset", Some(input))
        .options(args);
    let outcome = draw(&layout, modality, &args.site.cut(), &plan, sheet, out);
    record.finish(&outcome);
    for path in outcome? {
        println!("{path}");
    }
    Ok(())
}

/// The figures of `modality` in `layout`, saved under `out` in parallel;
/// the paths, in the figures' order.
fn draw(
    layout: &InputLayout,
    modality: &str,
    cut: &Cut,
    plan: &Plan,
    sheet: Sheet,
    out: &str,
) -> anyhow::Result<Vec<String>> {
    use rayon::prelude::*;
    let t = read_site_table(&layout.site_tables[modality], modality)?;
    anyhow::ensure!(
        !t.is_empty(),
        "{modality}_sites.parquet has no sites to draw"
    );
    let n_cells = site_cells(layout, &t)?;
    let figures = Panels::new(&t, n_cells.as_deref(), cut, plan, sheet, true)?.figures();
    anyhow::ensure!(
        !figures.is_empty(),
        "nothing to draw: --hist none and --pairs none"
    );
    std::fs::create_dir_all(out)?;
    let saved: Vec<Vec<String>> = (figures.par_iter())
        .map(|f| crate::figure::save_with_png(&f.svg, &format!("{out}/{}", f.stem)))
        .collect::<anyhow::Result<_>>()?;
    let paths: Vec<String> = saved.into_iter().flatten().collect();
    log::info!(
        "wrote {} figures of {modality} sites under {out}",
        figures.len()
    );
    Ok(paths)
}

#[cfg(test)]
#[path = "tests/site_plot.rs"]
mod tests;
