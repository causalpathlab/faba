//! `faba qc-view`: what a `faba qc` run kept and dropped, in its site view,
//! read only.
//!
//! A qc output carries the thresholds it applied (`options.site` of its
//! `qc.run.json`), the kept sites (`{modality}_sites.parquet`) and the
//! dropped ones with why (`{modality}_sites_dropped.parquet`, its `reason`
//! column), and each batch's cell QC (`{batch}_cell_qc_report.tsv`). The
//! view reads them back together, so every site the run saw is shown, kept
//! or dropped, with the applied thresholds drawn; it applies nothing and
//! writes nothing but the figures and the page asked for. Any other faba
//! output directory opens too, with no thresholds: every site is kept, as
//! `faba site-plot` draws it.

use std::path::{Path, PathBuf};

use anyhow::Context;
use arrow::array::{Array, StringArray};
use arrow::record_batch::RecordBatch;
use rustc_hash::{FxHashMap, FxHashSet};
use serde_json::Value;

use super::args::SiteFilterArgs;
use super::cells_plot::{Cells, Verdict};
use super::layout::{scan_input_dir, InputLayout, SITE_MODALITIES};
use super::matrix::{open_matrix, row_nnz_over_columns, select_columns};
use super::run::{cell_rule, counts_table, write_dry_plots, write_html};
use super::site_plot::PanelArgs;
use super::site_tui::{qc_flags, run_site_viewer, Batch, CellsShared};
use super::sites::{accumulate_site_cells, site_table_of, Criterion, SiteCells, SiteTable};
use crate::common::{read_lines, Args};
use crate::figure::style::{Sheet, StyleArgs};
use crate::figure::SizeArgs;

/// The run record a qc output carries.
const RECORD: &str = "qc.run.json";

#[derive(Args, Debug, serde::Serialize)]
pub struct QcViewArgs {
    #[arg(
        help = "A `faba qc` output directory (or its qc.run.json); any faba output directory opens too, with no thresholds"
    )]
    pub input_dir: Box<str>,

    #[arg(
        short = 'o',
        long = "output",
        requires = "plots",
        value_hint = clap::ValueHint::DirPath,
        help = "With --plots: the directory the figures go to, under qc_plots/ (new or empty)"
    )]
    pub output: Option<Box<str>>,

    #[arg(
        long = "batch-process",
        default_value_t = false,
        help = "Print the threshold table instead of opening the view",
        long_help = "Print, as TSV on stdout, what the cell QC kept per batch and what each\n\
                     threshold applied removes per modality (alone, and only by it), and the\n\
                     pairs of thresholds Auto picks to show two at a time, as `faba qc --dry-run`\n\
                     prints them, instead of opening the view; with --html, also the page, with\n\
                     --plots, also the figures. Without a terminal the table is printed as with it."
    )]
    pub batch_process: bool,

    /// `--no-detach` and `--session`: how the view runs as a session.
    #[command(flatten)]
    #[serde(skip)]
    pub session: legume_tui::cli::SessionArgs,

    #[arg(
        long = "plots",
        default_value_t = false,
        help = "Write the figures `faba qc --dry-run --plots` writes under -o DIR/qc_plots/, instead of opening the view",
        long_help = "Write, under -o DIR/qc_plots/ (DIR new or empty), the figures `faba qc\n\
                     --dry-run --plots` writes, at the thresholds applied: each modality's and\n\
                     knob's figure, each modality's pairs and each batch's cell QC, and the site\n\
                     panels under qc_plots/panels/ (--hist, --pairs, --marginals, --scale, --bins\n\
                     choose them), instead of opening the view."
    )]
    pub plots: bool,

    #[arg(
        long = "html",
        value_name = "FILE",
        value_hint = clap::ValueHint::FilePath,
        help = "Write the page Ctrl+E writes to FILE instead of opening the view (FILE must not exist)"
    )]
    pub html: Option<Box<str>>,

    /// The size of what `s` saves, and of the figures --plots writes.
    #[command(flatten)]
    pub size: SizeArgs,

    /// The font, panel and palette of what `s` saves and --plots writes.
    #[command(flatten)]
    pub style: StyleArgs,

    /// The site panels --plots writes, as `faba site-plot` draws them;
    /// `--scale` and `--bins` also set the view's axes, and `--marginals`
    /// what `s` saves of a pair.
    #[command(flatten)]
    pub panels: PanelArgs,

    #[arg(
        long = "gff",
        help = "Annotation for the view's metagene (default: the one recorded in the directory)"
    )]
    pub gff: Option<Box<str>>,

    #[arg(
        long = "genes",
        value_delimiter = ',',
        help = "Genes to list first in the view's gene list (symbols or gene keys)"
    )]
    pub genes: Vec<Box<str>>,
}

impl QcViewArgs {
    /// Whether it writes or prints instead of opening the view.
    fn batch(&self) -> bool {
        self.html.is_some() || self.plots || self.batch_process
    }

    /// The session the view runs in; none when it writes or prints
    /// instead. The directory is always named, so it starts a new one.
    pub fn view_session(&self) -> Option<legume_tui::Options> {
        (!self.batch()).then(|| self.session.options(true))
    }
}

/// A directory as the view opens it: every site of each modality, kept and
/// dropped, with its cells per site, and the thresholds applied.
pub struct Opened {
    pub dir: String,
    pub tables: FxHashMap<Box<str>, SiteTable>,
    /// Per modality, the kept cells with a converted read at each site.
    pub site_cells: SiteCells,
    /// As the qc run applied them; with no run, none.
    pub filter: SiteFilterArgs,
    /// Whether a qc run's record gave `filter`.
    pub applied: bool,
    /// The dropped sites read, per modality.
    pub dropped: FxHashMap<Box<str>, usize>,
    pub gff: Option<Box<str>>,
    /// Each batch's cell QC, from its report.
    pub cells: CellsShared,
    /// How the cell QC decided, as the table says it.
    pub cell_rule: String,
}

/// The thresholds a qc run recorded applying (`options.site`), each read
/// back as the `f32` or count it was; one written as null (an infinity)
/// is off.
pub fn recorded_thresholds(record: &Value) -> anyhow::Result<SiteFilterArgs> {
    let site = (record.get("options"))
        .and_then(|o| o.get("site"))
        .and_then(Value::as_object)
        .context("it records no site thresholds")?;
    let mut f = SiteFilterArgs::permissive();
    for c in Criterion::ALL {
        let key = c.flag().trim_start_matches("--").replace('-', "_");
        let v = match site.get(&key) {
            Some(Value::Null) => continue,
            Some(v) => v
                .as_f64()
                .with_context(|| format!("{key} is {v}, not a number"))?,
            None => anyhow::bail!("its site thresholds lack {key}"),
        };
        // The record holds each `f32` widened exactly, so it sets unrounded.
        c.set(&mut f, v);
    }
    Ok(f)
}

/// `kept`'s rows, then `dropped`'s, by `kept`'s columns (the `reason` of
/// a dropped site left out).
fn joined(kept: RecordBatch, dropped: &RecordBatch) -> anyhow::Result<RecordBatch> {
    let schema = kept.schema();
    let columns = (schema.fields().iter())
        .map(|f| {
            (dropped.column_by_name(f.name()).cloned())
                .with_context(|| format!("the dropped sites lack the `{}` column", f.name()))
        })
        .collect::<anyhow::Result<Vec<_>>>()?;
    let dropped = RecordBatch::try_new(schema.clone(), columns)?;
    Ok(arrow::compute::concat_batches(&schema, [&kept, &dropped])?)
}

/// The `reason` of each dropped site.
fn reasons(dropped: &RecordBatch) -> anyhow::Result<Vec<Box<str>>> {
    let r: &StringArray = crate::common::column(dropped, "reason")?;
    Ok((0..r.len()).map(|i| r.value(i).into()).collect())
}

/// Each batch's cells, by name.
type Names = FxHashMap<Box<str>, FxHashSet<Box<str>>>;

/// Per modality with a table, the cells of `layout`'s `_site` matrices
/// with a converted read at each site: of each batch, those `keep` names
/// (every cell with a read where it names none); and the same over the
/// `controls` batches alone, the control arm's.
fn count_site_cells(
    layout: &InputLayout,
    tables: &FxHashMap<Box<str>, SiteTable>,
    (keep, controls): (&Names, &FxHashSet<Box<str>>),
) -> anyhow::Result<(SiteCells, SiteCells)> {
    let mut out: SiteCells = FxHashMap::default();
    let mut control: SiteCells = FxHashMap::default();
    for m in &layout.matrices {
        let Some(t) = m.site_modality().and_then(|s| tables.get(s)) else {
            continue;
        };
        let data = open_matrix(&m.path)?;
        let cols = select_columns(data.as_ref(), &data.column_names()?, keep.get(&m.batch));
        let nnz = row_nnz_over_columns(data.as_ref(), &cols)?;
        let rows = data.row_names()?;
        let acc = out.entry(t.modality.clone()).or_default();
        accumulate_site_cells(&rows, &nnz, t.converted_channel(), acc);
        if controls.contains(&m.batch) {
            let acc = control.entry(t.modality.clone()).or_default();
            accumulate_site_cells(&rows, &nnz, t.converted_channel(), acc);
        }
    }
    Ok((out, control))
}

/// The batch a file in a qc output is the `suffix` of, when it is.
fn batch_of<'a>(path: &'a str, suffix: &str) -> Option<&'a str> {
    let name = Path::new(path).file_name()?.to_str()?;
    name.strip_suffix(suffix).filter(|b| !b.is_empty())
}

/// Each batch's kept cells, from the qc output's `{batch}_cells.tsv.gz`.
fn kept_cells(layout: &InputLayout) -> anyhow::Result<FxHashMap<Box<str>, FxHashSet<Box<str>>>> {
    let mut out = FxHashMap::default();
    for f in &layout.other_files {
        if let Some(b) = batch_of(f, "_cells.tsv.gz") {
            out.insert(b.into(), read_lines(f)?.into_iter().collect());
        }
    }
    Ok(out)
}

/// A batch's barcodes as its cell QC report (`write_qc_report`'s TSV)
/// gives them: counts, detected genes and the verdict.
fn read_cells(path: &str, batch: &str) -> anyhow::Result<Cells> {
    let lines = read_lines(path)?;
    let header: Vec<&str> = (lines.first().map_or("", |h| h.trim_start_matches('#')))
        .split('\t')
        .collect();
    let col = |name: &str| {
        (header.iter().position(|h| *h == name))
            .with_context(|| format!("{path}: no `{name}` column"))
    };
    let [genes, counts, empty, keep] = [
        col("n_genes")?,
        col("total_counts")?,
        col("near_empty")?,
        col("train_keep")?,
    ];
    let mut cells = Cells {
        batch: batch.into(),
        counts: Vec::new(),
        genes: Vec::new(),
        verdict: Vec::new(),
    };
    for line in lines.iter().skip(1).filter(|l| !l.is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        let at = |i: usize| f.get(i).copied().unwrap_or_default();
        let real = |i: usize| {
            (at(i).parse::<f32>()).with_context(|| format!("{path}: `{line}` is not a report row"))
        };
        cells.counts.push(real(counts)?);
        cells.genes.push(real(genes)?);
        cells.verdict.push(match (at(empty), at(keep)) {
            ("1", _) => Verdict::NearEmpty,
            (_, "0") => Verdict::Outlier,
            _ => Verdict::Kept,
        });
    }
    Ok(cells)
}

/// Warned when the verdicts at `filter` are not those recorded: the kept
/// sites (the first `n_kept` of `t`) pass, each dropped one fails first
/// the check its `reason` names.
fn check_reasons(
    t: &SiteTable,
    n_cells: Option<&[usize]>,
    filter: &SiteFilterArgs,
    (n_kept, recorded): (usize, &[Box<str>]),
) {
    let now = filter.reasons(t, n_cells);
    let differ = (now.iter().enumerate())
        .filter(|&(i, r)| {
            let was = (i >= n_kept).then(|| &*recorded[i - n_kept]);
            r.map(|r| r.label()) != was
        })
        .count();
    if differ > 0 {
        log::warn!(
            "{}: {differ} of {} sites are not cut at the thresholds recorded as qc cut them",
            t.modality,
            t.len()
        );
    }
}

/// `input` (a directory, or its qc.run.json) as the view opens it, the
/// metagene's annotation `gff` or the one recorded.
pub fn open(input: &str, gff: Option<&str>) -> anyhow::Result<Opened> {
    let dir = crate::run_record::output_dir_of(input)
        .with_context(|| format!("{input}: not a faba output directory or its {RECORD}"))?;
    let dir_s = dir.to_string_lossy().into_owned();
    let layout = scan_input_dir(&dir_s)?;
    let record_path = dir.join(RECORD);
    let record: Option<Value> = match record_path.is_file() {
        true => Some(
            serde_json::from_str(&std::fs::read_to_string(&record_path)?)
                .with_context(|| format!("reading {}", record_path.display()))?,
        ),
        false => None,
    };
    let filter = match &record {
        Some(r) => recorded_thresholds(r)
            .with_context(|| format!("{}: not a finished qc run's", record_path.display()))?,
        None => SiteFilterArgs::permissive(),
    };
    // Each modality's sites, the dropped ones after the kept with what
    // they were dropped for.
    let mut tables = FxHashMap::default();
    let mut dropped: FxHashMap<Box<str>, (usize, Vec<Box<str>>)> = FxHashMap::default();
    for m in SITE_MODALITIES {
        let Some(path) = layout.site_tables.get(*m) else {
            continue;
        };
        let mut batch = crate::common::read_record_batch(&**path)?;
        let gone = dir.join(format!("{m}_sites_dropped.parquet"));
        if record.is_some() && gone.is_file() {
            let d = crate::common::read_record_batch(&gone)?;
            dropped.insert((*m).into(), (batch.num_rows(), reasons(&d)?));
            batch = joined(batch, &d)?;
        }
        tables.insert(Box::from(*m), site_table_of(batch, m)?);
    }
    let controls = super::sites::control_batches(&dir);
    let (site_cells, control_cells) = match record
        .as_ref()
        .map(|_| site_cells_of_run(&dir, &layout, &tables, &controls))
    {
        Some(counted) => counted?,
        None => count_site_cells(&layout, &tables, (&FxHashMap::default(), &controls))?,
    };
    for (m, t) in tables.iter_mut() {
        super::sites::set_control_cells(t, control_cells.get(m));
    }
    for (m, (n_kept, why)) in &dropped {
        let t = &tables[m];
        let n = site_cells.get(m).map(|acc| t.cells_per_site(acc));
        check_reasons(t, n.as_deref(), &filter, (*n_kept, why));
    }
    let mut reports: Vec<(&str, &str)> = (layout.other_files.iter())
        .filter_map(|f| Some((batch_of(f, "_cell_qc_report.tsv")?, &**f)))
        .collect();
    reports.sort_unstable();
    let cells = (reports.iter())
        .map(|(b, f)| read_cells(f, b))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let option = |key: &str| (record.as_ref()).and_then(|r| r.get("options")?.get(key).cloned());
    let rule = cell_rule(
        option("no_cell_qc")
            .and_then(|v| v.as_bool())
            .unwrap_or(false),
        option("qc_mads").and_then(|v| v.as_f64()).unwrap_or(5.0) as f32,
    );
    Ok(Opened {
        gff: crate::run_record::explicit_or_recorded(gff, &dir_s, "gff", "annotation"),
        dir: dir_s,
        tables,
        site_cells,
        filter,
        applied: record.is_some(),
        dropped: dropped
            .into_iter()
            .map(|(m, (_, r))| (m, r.len()))
            .collect(),
        cells: cells.into(),
        cell_rule: rule,
    })
}

/// Cells per site as the qc run in `dir` counted them: over its input's
/// `_site` matrices, each batch's cells those it kept. With the input not
/// on this machine, over `dir`'s own, which hold the kept sites only:
/// every dropped site then reads 0 cells, as it is warned.
fn site_cells_of_run(
    dir: &Path,
    layout: &InputLayout,
    tables: &FxHashMap<Box<str>, SiteTable>,
    controls: &FxHashSet<Box<str>>,
) -> anyhow::Result<(SiteCells, SiteCells)> {
    let input = crate::run_record::input_of(&dir.join(RECORD), "fileset")
        .map(|p| PathBuf::from(&*p))
        .filter(|p| p.is_dir());
    let read = input
        .as_deref()
        .map(|p| scan_input_dir(&p.to_string_lossy()))
        .transpose()?;
    match read {
        Some(input) => count_site_cells(&input, tables, (&kept_cells(layout)?, controls)),
        None => {
            let own = count_site_cells(layout, tables, (&FxHashMap::default(), controls))?;
            if !own.0.is_empty() {
                log::warn!(
                    "{}: the directory qc read is not on this machine; cells per site are \
                     counted off this one's `_site` matrices, so every dropped site reads 0",
                    dir.display()
                );
            }
            Ok(own)
        }
    }
}

/// Refused, before anything is read, when `--html` names a file that
/// exists (never overwritten) or `--plots` a folder that is not new or
/// empty, or the figures' page is too small for them.
fn problems(args: &QcViewArgs, sheet: &Sheet) -> anyhow::Result<()> {
    if let Some(html) = args.html.as_deref() {
        anyhow::ensure!(
            !Path::new(html).exists(),
            "--html {html} already exists: faba qc-view never overwrites a file; name another path"
        );
    }
    if args.plots {
        let Some(out) = args.output.as_deref() else {
            anyhow::bail!(
                "`faba qc-view --plots` writes its figures under -o DIR: name one, new or empty"
            );
        };
        if let Some(problem) = crate::tui::output_problem(out) {
            anyhow::bail!("-o {problem}: faba qc-view never writes over files; name another");
        }
        use super::site_tui::{knob_least, knob_room};
        knob_room(sheet, knob_least(sheet, None), "a qc figure")?;
    }
    Ok(())
}

pub fn run_qc_view(args: &QcViewArgs) -> anyhow::Result<()> {
    let plan = args.panels.plan()?;
    let sheet = Sheet::new(crate::figure::size::Figure::Qc, &args.size, &args.style)?;
    problems(args, &sheet)?;
    let o = open(&args.input_dir, args.gff.as_deref())?;
    if !o.applied {
        log::info!("{}: no qc run here; no thresholds, every site kept", o.dir);
    }
    for (m, n) in &o.dropped {
        let (all, at) = (o.tables[m].len(), qc_flags(&o.filter));
        log::info!("{m}: {all} sites, {n} of them dropped by qc at {at}");
    }
    let genes = (o.gff.as_deref(), &args.genes[..]);
    let tables = (&o.tables, &o.site_cells);
    // Without a terminal, the table is printed as with --batch-process.
    let writes = args.html.is_some() || args.plots;
    if !writes && !crate::tui::batch(args.batch_process, "faba qc-view") {
        let viewed = (plan, sheet);
        return run_site_viewer(&o.dir, tables, &o.filter, genes, viewed, o.cells.clone());
    }
    let scales = (plan.scales, plan.bins_on(&sheet));
    let mut batch = Batch::new(
        &o.dir,
        &o.tables,
        &o.site_cells,
        &o.filter,
        genes,
        (o.cells.clone(), Some(scales)),
    )?;
    if args.batch_process || !writes {
        let kept: Vec<(&str, usize, usize)> = (o.cells.iter())
            .map(|c| {
                (
                    &*c.batch,
                    c.tally()[Verdict::Kept as usize],
                    c.verdict.len(),
                )
            })
            .collect();
        print!("{}", counts_table((&o.cell_rule, &kept), &batch));
    }
    if let Some(out) = args.output.as_deref().filter(|_| args.plots) {
        write_dry_plots(out, &mut batch, (&plan, sheet))?;
    }
    write_html(args.html.as_deref(), &mut batch)
}

#[cfg(test)]
#[path = "tests/view.rs"]
mod tests;
