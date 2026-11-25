//! `faba qc`: read a faba output directory, decide cells, features and sites,
//! and write a new, filtered fileset. Never in place.
//!
//! Order of operations, per batch:
//! 1. cells are decided once, on `{batch}_count`, and the same keep set is
//!    applied to every matrix of the batch so modalities stay column-aligned;
//! 2. sites are decided on the parquet columns plus the per-site kept-cell
//!    count read off the `_site` matrices (pooled over batches);
//! 3. `_site` matrices keep both channels of every kept site;
//! 4. gene-level `{batch}_m6a` / `{batch}_atoi` are RE-POOLED from the filtered
//!    site matrix, so they cannot disagree with the site cut;
//! 5. every other matrix takes the cell keep set and `--row-nnz-cutoff`;
//! 6. every other file is copied through;
//! 7. with the view, its figure for every modality and knob goes to
//!    `qc_plots/` as PDF and PNG.

use std::sync::Arc;

use crate::common::*;
use crate::run_record::{explicit_or_recorded, find_input, RunRecord};
use data_beans::qc_lib::{compute_qc, write_qc_report, QcConfig, QcReport};
use data_beans::sparse_io_vector::SparseIoVec;
use rustc_hash::{FxHashMap, FxHashSet};

use super::args::{QcArgs, SiteFilterArgs};
use super::layout::{file_name, scan_input_dir, InputLayout, MatrixFile, SITE_MODALITIES};
use super::matrix::{
    open_matrix, row_nnz_over_columns, select_columns, shape, write_subset, Backend, OutSpec,
    Written,
};
use super::path_tui::ask_paths;
use super::progress::{FinishOnDrop, Progress};
use super::repool::repool_gene_level;
use super::site_plot::Plan;
use super::site_tui::{qc_flags, run_site_picker, Batch, CellsShared, Figure, Picked, Writer};
use super::sites::{
    accumulate_site_cells, read_site_table, site_matrix_rows, write_site_tables, SiteTable,
};
use legume_tui::view::terminal_available;

struct SummaryRow {
    file: Box<str>,
    before: (usize, usize, usize),
    after: (usize, usize, usize),
}

/// One batch's cell verdicts, held until nothing can cancel the run.
struct CellDecision {
    batch: Box<str>,
    names: Vec<Box<str>>,
    report: QcReport,
    kept: Vec<Box<str>>,
}

impl CellDecision {
    fn write(&self, out_dir: &str) -> anyhow::Result<()> {
        let report_path = format!("{out_dir}/{}_cell_qc_report.tsv", self.batch);
        write_qc_report(&report_path, &self.names, &self.report)?;
        write_lines(
            &self.kept,
            &format!("{out_dir}/{}_cells.tsv.gz", self.batch),
        )
    }

    fn keep_set(&self) -> FxHashSet<Box<str>> {
        self.kept.iter().cloned().collect()
    }
}

/// Decide the cells of one batch on its `_count` matrix, through the
/// data-beans cell QC in one pass: a near-empty floor (the larger of
/// `--column-nnz-cutoff` and `--qc-min-cell-nnz`), the 2-means suggestion
/// under `--auto-cutoff`, and the MAD-outlier drops unless `--no-cell-qc`.
/// The per-cell verdicts and kept barcodes are returned for [`CellDecision::write`].
fn decide_cells(m: &MatrixFile, args: &QcArgs) -> anyhow::Result<CellDecision> {
    let data = open_matrix(&m.path)?;
    let names = data.column_names()?;
    let mut vec = SparseIoVec::new();
    vec.push(Arc::from(data), None)?;
    let floor = if args.no_cell_qc {
        1
    } else {
        args.qc_min_cell_nnz
    };
    let cfg = QcConfig {
        n_mads: args.qc_mads,
        min_cell_nnz: args.column_nnz_cutoff.max(floor),
        mad_on_n_genes: !args.no_cell_qc,
        mad_on_counts: !args.no_cell_qc,
        auto_cell_cutoff: args.auto_cutoff,
        qc_histogram: args.show_histogram,
        ..QcConfig::default()
    };
    let report = compute_qc(&vec, &cfg, args.block_size)?;

    let idx = report.emit_idx_unmasked();
    anyhow::ensure!(
        !idx.is_empty(),
        "{}: no cell survives the cell QC; relax --column-nnz-cutoff / --qc-min-cell-nnz or pass --no-cell-qc",
        m.batch
    );
    let kept: Vec<Box<str>> = idx.iter().map(|&c| names[c].clone()).collect();
    info!(
        "{}: kept {} of {} cells (near-empty floor {}, MAD QC {})",
        m.batch,
        kept.len(),
        names.len(),
        cfg.min_cell_nnz,
        if args.no_cell_qc { "off" } else { "on" }
    );
    Ok(CellDecision {
        batch: m.batch.clone(),
        names,
        report,
        kept,
    })
}

/// `data-beans squeeze`'s rule for a feature-axis nnz cutoff: an explicit
/// non-zero value wins, else the 2-means suggestion under `--auto-cutoff`,
/// else nothing. Empty rows always drop.
fn row_cutoff(args: &QcArgs, nnz: &[usize], label: &str) -> usize {
    let nnz_f: Vec<f32> = nnz.iter().map(|&x| x as f32).collect();
    let suggested = (args.auto_cutoff || args.show_histogram)
        .then(|| suggest_nnz_cutoff(&nnz_f))
        .flatten();
    let cutoff = if args.row_nnz_cutoff > 0 {
        args.row_nnz_cutoff
    } else if args.auto_cutoff {
        suggested.unwrap_or(0)
    } else {
        0
    };
    if args.show_histogram {
        print_nnz_summary(label, "nnz", &nnz_f, cutoff, suggested);
    } else if args.auto_cutoff {
        info!("{label}: row nnz cutoff {cutoff} (2-means suggestion {suggested:?})");
    }
    cutoff.max(1)
}

fn out_spec<'a>(m: &'a MatrixFile, out_dir: &'a str, stem: &'a str, no_zip: bool) -> OutSpec<'a> {
    OutSpec {
        out_dir,
        stem,
        backend: &m.backend,
        zip: m.zipped && !no_zip,
    }
}

fn record(
    summary: &mut Vec<SummaryRow>,
    stem: &str,
    before: (usize, usize, usize),
    w: Option<Written>,
) {
    match w {
        Some(w) => {
            info!(
                "{stem}: {}x{} ({} nnz) -> {}x{} ({} nnz)",
                before.0, before.1, before.2, w.nrow, w.ncol, w.nnz
            );
            summary.push(SummaryRow {
                file: file_name(&w.target),
                before,
                after: (w.nrow, w.ncol, w.nnz),
            });
        }
        None => {
            log::warn!("{stem}: nothing survives the cut; not written");
            summary.push(SummaryRow {
                file: format!("{stem} (not written)").into(),
                before,
                after: (0, 0, 0),
            });
        }
    }
}

/// A matrix opened once with its axis names and shape.
struct Opened {
    data: Backend,
    row_names: Vec<Box<str>>,
    col_names: Vec<Box<str>>,
    before: (usize, usize, usize),
}

fn open_with_names(m: &MatrixFile) -> anyhow::Result<Opened> {
    let data = open_matrix(&m.path)?;
    let row_names = data.row_names()?;
    let col_names = data.column_names()?;
    let before = shape(data.as_ref());
    Ok(Opened {
        data,
        row_names,
        col_names,
        before,
    })
}

/// A `_site` matrix opened in step 2 and written in steps 3-4.
struct SiteMatrix {
    m: MatrixFile,
    opened: Opened,
    cols: Vec<usize>,
}

/// The input and output directories: as given, else asked for in a pop-up.
fn resolve_paths(args: &QcArgs) -> anyhow::Result<(String, String)> {
    if let (Some(input), Some(output)) = (&args.input_dir, &args.output) {
        return Ok((input.to_string(), output.to_string()));
    }
    anyhow::ensure!(
        !args.batch_process,
        "`faba qc --batch-process` needs INPUT_DIR and -o/--output"
    );
    let Some((input, output)) = ask_paths(args.input_dir.as_deref(), args.output.as_deref())?
    else {
        anyhow::bail!("cancelled at the directories; nothing written");
    };
    info!("faba qc {input} -o {output}");
    Ok((input, output))
}

/// Refused, before anything is read, when `--html` cannot be written as
/// asked: without --batch-process or --dry-run (the view writes its page on
/// Ctrl+E), or to a file that exists (never overwritten; the write itself
/// refuses one that appears meanwhile).
fn html_problem(args: &QcArgs) -> anyhow::Result<()> {
    let Some(html) = args.html.as_deref() else {
        return Ok(());
    };
    anyhow::ensure!(
        args.batch_process || args.dry_run,
        "--html writes the page with --batch-process or --dry-run; in the view, Ctrl+E writes it"
    );
    anyhow::ensure!(
        !std::path::Path::new(html).exists(),
        "--html {html} already exists: faba qc never overwrites a file; name another path"
    );
    Ok(())
}

/// Refused, before anything is read, when a dry run's figures cannot go
/// where asked: `--dry-run --plots` writes them to `-o DIR`'s `qc_plots/`,
/// DIR new or empty (never a fileset, nothing else in it).
fn dry_plots_problem(args: &QcArgs) -> anyhow::Result<()> {
    if !(args.dry_run && args.plots) {
        return Ok(());
    }
    let Some(out) = args.output.as_deref() else {
        anyhow::bail!(
            "`faba qc --dry-run --plots` writes its figures under -o DIR: name one, new or empty"
        );
    };
    if let Some(problem) = crate::tui::output_problem(out) {
        anyhow::bail!("-o {problem}: faba qc never writes over files; name another");
    }
    Ok(())
}

/// A dry run's figures: each modality's knobs and pairs, each batch's
/// cells, and the site panels, under `-o DIR/qc_plots` (checked new or empty before anything
/// was read, see [`dry_plots_problem`]).
fn write_dry_plots(args: &QcArgs, batch: &mut Batch, (plan, sheet): Figures) -> anyhow::Result<()> {
    let out = args.output.as_deref().unwrap_or_default();
    let figures = batch.figures(Some(plan), sheet)?;
    let failed = save_figures(out, &figures, &Progress::default());
    if let Some((stem, e)) = failed.first() {
        anyhow::bail!("qc_plots/{stem}: not saved: {e}");
    }
    info!("wrote {} figures under {out}/qc_plots", figures.len());
    Ok(())
}

/// The figures' page and look, checked before anything is read: a knob
/// figure's table and plots need a width (their height, which depends on
/// the plots shown, is checked before the view opens and as each is
/// drawn), unless --panel-only draws the plot alone.
fn figure_sheet(args: &QcArgs) -> anyhow::Result<crate::figure::style::Sheet> {
    use crate::qc::site_tui::{knob_least, knob_room};
    let sheet = crate::figure::style::Sheet::new(&args.size, &args.style)?;
    knob_room(&sheet, knob_least(&sheet, None), "a qc figure")?;
    Ok(sheet)
}

/// The site panels asked for and the page and look of every figure, read
/// and checked once per run.
type Figures<'a> = (&'a Plan, crate::figure::style::Sheet);

/// [`Plan::bins_on`] its sheet: the scales and bins every figure reads.
fn scales_of((plan, sheet): Figures) -> (crate::qc::scale::Scales, usize) {
    (plan.scales, plan.bins_on(&sheet))
}

/// `--html`'s page, at the thresholds given, its folder made first; never
/// over a file (see [`html_problem`]).
fn write_html(args: &QcArgs, batch: &mut Batch) -> anyhow::Result<()> {
    let Some(html) = args.html.as_deref() else {
        return Ok(());
    };
    let path = std::path::Path::new(html);
    if let Some(parent) = path.parent().filter(|p| !p.as_os_str().is_empty()) {
        std::fs::create_dir_all(parent)?;
    }
    if batch.write_page(path)? {
        info!("wrote {html}");
    }
    Ok(())
}

pub fn run_qc(args: &QcArgs) -> anyhow::Result<()> {
    html_problem(args)?;
    let plan = args.panels.plan()?;
    let sheet = figure_sheet(args)?;
    dry_plots_problem(args)?;
    if args.dry_run {
        return dry_run(args, (&plan, sheet));
    }
    anyhow::ensure!(
        args.batch_process || terminal_available(),
        "`faba qc` picks the site thresholds in a full-screen view, which needs stdin and \
         stdout on a terminal; pass --batch-process to cut with the --site-* values as given"
    );
    let (input, output) = resolve_paths(args)?;
    let (input, out_dir) = (input.as_str(), output.as_str());
    if std::path::Path::new(out_dir).exists() && std::fs::read_dir(out_dir)?.next().is_some() {
        anyhow::bail!("output directory {out_dir} already contains files; choose an empty one");
    }
    // Carry the annotation and genome forward, so tools reading the new
    // fileset find them as they would in the original.
    let mut gff = explicit_or_recorded(args.gff.as_deref(), input, "gff", "annotation");
    let genome = find_input(input, "genome", "genome");
    let record = RunRecord::start("qc", out_dir)
        .input("fileset", Some(input))
        .input("genome", genome.as_deref())
        .options(args);
    let given = gff.clone();
    let outcome = filter_fileset(args, input, out_dir, &mut gff, (&plan, sheet));
    // Recorded after the view, which may have read another annotation.
    let mut record = record.input("gff", gff.as_deref());
    if gff != given {
        record.set_option("gff", &gff);
    }
    // A run cancelled before writing leaves the directory empty, for a rerun.
    let wrote = std::fs::read_dir(out_dir).is_ok_and(|mut d| d.next().is_some());
    if wrote {
        // The directories, which may have been asked for, and the thresholds
        // applied, which the view may have changed.
        record.set_option("input_dir", &input);
        record.set_option("output", &out_dir);
        if let Ok(Some(site)) = &outcome {
            record.set_option("site", site);
        }
        record.finish(&outcome);
    }
    outcome.map(|_| ())
}

/// `faba qc --dry-run`: the cell verdicts per batch and each site
/// threshold's counts per modality, as TSV on stdout; then, after a blank
/// line, the pairs Auto picks to show two at a time, a TSV of their own.
/// No fileset; the figures and the page when asked.
fn dry_run(args: &QcArgs, figures: Figures) -> anyhow::Result<()> {
    let input = args
        .input_dir
        .as_deref()
        .ok_or_else(|| anyhow::anyhow!("`faba qc --dry-run` needs INPUT_DIR"))?;
    let prep = prepare(args, input)?;
    let mut out = String::from("level\tname\trule\tthreshold\tremoves\tonly_this\tkept\tof\n");
    for d in &prep.decisions {
        let (n, kept) = (d.names.len(), d.kept.len());
        out.push_str(&format!(
            "cells\t{}\tcell QC\t{}\t{}\t\t{kept}\t{n}\n",
            d.batch,
            if args.no_cell_qc {
                "floor only".to_string()
            } else {
                format!("{} MADs", args.qc_mads)
            },
            n - kept,
        ));
    }
    let gff = explicit_or_recorded(args.gff.as_deref(), input, "gff", "annotation");
    let mut batch = Batch::new(
        input,
        &prep.tables,
        &prep.site_cells,
        &args.site,
        (gff.as_deref(), &args.genes),
        (cells_of(&prep.decisions), Some(scales_of(figures))),
    )?;
    let (counts, pairs) = batch.counts();
    for t in counts {
        out.push_str(&format!(
            "sites\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n",
            t.modality, t.flag, t.threshold, t.removes, t.only, t.kept, t.sites
        ));
    }
    // The pairs Auto picks for the 2D panels: `split`, the sites one of
    // the pair's thresholds keeps and the other drops.
    out.push_str("\nmodality\tx\ty\trank_correlation\tsplit\n");
    for (m, p) in pairs {
        out.push_str(&format!(
            "{m}\t{}\t{}\t{:.2}\t{}\n",
            p.x.short(),
            p.y.short(),
            p.rho,
            p.split
        ));
    }
    print!("{out}");
    if args.plots {
        write_dry_plots(args, &mut batch, figures)?;
    }
    write_html(args, &mut batch)
}

/// `figures` saved under `out/qc_plots/` as PDF and PNG (a stem's
/// folder made first), in parallel, each a step on `progress`; those not
/// saved, and why.
fn save_figures(out: &str, figures: &[Figure], progress: &Progress) -> Vec<(String, String)> {
    use rayon::prelude::*;
    if figures.is_empty() {
        return Vec::new();
    }
    let dir = format!("{out}/qc_plots");
    let made = std::fs::create_dir_all(&dir).map_err(|e| e.to_string());
    (figures.par_iter())
        .filter_map(|f| {
            progress.next(format!("qc_plots/{}", f.stem));
            let prefix = format!("{dir}/{}", f.stem);
            let saved = made.clone().and_then(|()| {
                let folder = std::path::Path::new(&prefix).parent();
                folder
                    .map_or(Ok(()), std::fs::create_dir_all)
                    .map_err(|e| e.to_string())?;
                crate::figure::save_with_png(&f.svg, &prefix).map_err(|e| e.to_string())
            });
            saved.err().map(|e| (f.stem.clone(), e))
        })
        .collect()
}

/// Each batch's barcodes and the cell QC's verdicts, for the figures and
/// the page.
fn cells_of(decisions: &[CellDecision]) -> CellsShared {
    (decisions.iter())
        .map(|d| crate::qc::cells_plot::Cells::of(&d.batch, &d.report))
        .collect()
}

/// Everything read before the site thresholds are decided: the input's
/// layout, the cell verdicts, the site tables and the `_site` matrices with
/// their kept cells per site.
struct Prepared {
    layout: InputLayout,
    cells: FxHashMap<Box<str>, FxHashSet<Box<str>>>,
    decisions: Vec<CellDecision>,
    tables: FxHashMap<Box<str>, SiteTable>,
    site_matrices: Vec<SiteMatrix>,
    site_cells: FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>,
}

/// [`run_qc`] once the output directory is known to be empty, with `gff` for
/// the view's metagene, which the view may replace. Returns the site
/// thresholds applied, `None` when nothing was written.
fn filter_fileset(
    args: &QcArgs,
    input: &str,
    out_dir: &str,
    gff: &mut Option<Box<str>>,
    (plan, sheet): Figures,
) -> anyhow::Result<Option<SiteFilterArgs>> {
    let scales = scales_of((plan, sheet));
    let prep = prepare(args, input)?;
    let has_sites = prep.tables.values().any(|t| !t.is_empty());
    if args.batch_process || !has_sites {
        if !has_sites && !args.batch_process {
            log::warn!(
                "no editing site table to pick thresholds on; cutting cells and features only"
            );
        }
        let cells = cells_of(&prep.decisions);
        let mut batch = match args.plots || args.html.is_some() {
            true => Some(Batch::new(
                input,
                &prep.tables,
                &prep.site_cells,
                &args.site,
                (gff.as_deref(), &args.genes),
                (cells, Some(scales)),
            )?),
            false => None,
        };
        let figures = match (&mut batch, args.plots) {
            (Some(b), true) => b.figures(Some(plan), sheet)?,
            _ => Vec::new(),
        };
        write_fileset(
            &prep,
            args,
            out_dir,
            &args.site,
            &figures,
            &Progress::default(),
        )?;
        if let Some(b) = &mut batch {
            write_html(args, b)?;
        }
        return Ok(Some(args.site.clone()));
    }
    let progress = Progress::default();
    let picked = std::thread::scope(|s| -> anyhow::Result<_> {
        let mut writer = None;
        let (picked, picked_gff) = run_site_picker(
            input,
            &prep.tables,
            &prep.site_cells,
            &args.site,
            (gff.as_deref(), &args.genes),
            scales,
            Writer {
                output: out_dir,
                progress: &progress,
                sheet,
                panels: Some(plan.clone()),
                cells: cells_of(&prep.decisions),
                start: Some(Box::new(|site: SiteFilterArgs, figures: Vec<Figure>| {
                    let (prep, progress) = (&prep, &progress);
                    writer = Some(s.spawn(move || {
                        let _finish = FinishOnDrop(progress);
                        write_fileset(prep, args, out_dir, &site, &figures, progress)
                    }));
                })),
            },
        )?;
        *gff = picked_gff;
        if let Some(w) = writer {
            w.join()
                .map_err(|_| anyhow::anyhow!("writing the filtered fileset panicked"))??;
        }
        Ok(picked)
    })?;
    match picked {
        Picked::Apply(f) => {
            info!("site thresholds: {}", qc_flags(&f));
            Ok(Some(f))
        }
        Picked::PrintOnly(f) => {
            let command = batch_command(args, input, out_dir, &f, gff.as_deref());
            crate::tui::work(|| println!("{command}"))?;
            Ok(None)
        }
        Picked::Cancelled => anyhow::bail!("cancelled at the site thresholds; nothing written"),
    }
}

/// The `faba qc --batch-process` command that cuts `input` into `out_dir` as
/// this run would, with the site thresholds `site` and the annotation `gff`
/// the view ended on.
fn batch_command(
    args: &QcArgs,
    input: &str,
    out_dir: &str,
    site: &SiteFilterArgs,
    gff: Option<&str>,
) -> String {
    let mut cmd = format!(
        "faba qc {} -o {} --batch-process -r {} -c {} --qc-mads {} --qc-min-cell-nnz {}",
        input,
        out_dir,
        args.row_nnz_cutoff,
        args.column_nnz_cutoff,
        args.qc_mads,
        args.qc_min_cell_nnz
    );
    for (on, flag) in [
        (args.auto_cutoff, "--auto-cutoff"),
        (args.no_cell_qc, "--no-cell-qc"),
        (args.show_histogram, "--show-histogram"),
        (args.no_zip, "--no-zip"),
    ] {
        if on {
            cmd += &format!(" {flag}");
        }
    }
    if let Some(b) = args.block_size {
        cmd += &format!(" --block-size {b}");
    }
    if let Some(gff) = gff {
        cmd += &format!(" --gff {gff}");
    }
    format!("{cmd} {}", qc_flags(site))
}

/// Steps 1-2 up to the site thresholds: read, decide cells, and count kept
/// cells per site. Nothing is written.
fn prepare(args: &QcArgs, input: &str) -> anyhow::Result<Prepared> {
    let layout: InputLayout = scan_input_dir(input)?;
    let batches = layout.batches();
    info!(
        "{}: {} matrices over {} batches, {} site tables, {} other files",
        input,
        layout.matrices.len(),
        batches.len(),
        layout.site_tables.len(),
        layout.other_files.len()
    );

    // 1. cells, per batch, on `_count`.
    let mut cells: FxHashMap<Box<str>, FxHashSet<Box<str>>> = FxHashMap::default();
    let mut decisions: Vec<CellDecision> = Vec::new();
    for b in &batches {
        match layout.matrix(b, "count") {
            Some(m) => {
                let d = decide_cells(m, args)?;
                cells.insert(b.clone(), d.keep_set());
                decisions.push(d);
            }
            None => log::warn!(
                "{b}: no `_count` matrix; keeping every non-empty cell of its other matrices"
            ),
        }
    }

    // 2. sites: parquet columns first, then cells per site from the `_site`
    //    matrices (pooled over batches).
    let mut tables: FxHashMap<Box<str>, SiteTable> = FxHashMap::default();
    for (modality, path) in &layout.site_tables {
        tables.insert(modality.clone(), read_site_table(path, modality)?);
    }
    let mut site_matrices: Vec<SiteMatrix> = Vec::new();
    let mut site_cells: FxHashMap<Box<str>, FxHashMap<Box<str>, usize>> = FxHashMap::default();
    for m in &layout.matrices {
        let Some(modality) = m.site_modality() else {
            continue;
        };
        let Some(t) = tables.get(modality) else {
            log::warn!(
                "{}: no {modality}_sites.parquet; the matrix is not written",
                m.stem()
            );
            continue;
        };
        let opened = open_with_names(m)?;
        let cols = select_columns(opened.data.as_ref(), &opened.col_names, cells.get(&m.batch));
        let row_nnz = row_nnz_over_columns(opened.data.as_ref(), &cols)?;
        accumulate_site_cells(
            &opened.row_names,
            &row_nnz,
            t.converted_channel(),
            site_cells.entry(modality.into()).or_default(),
        );
        site_matrices.push(SiteMatrix {
            m: m.clone(),
            opened,
            cols,
        });
    }
    Ok(Prepared {
        layout,
        cells,
        decisions,
        tables,
        site_matrices,
        site_cells,
    })
}

/// Steps 3-6 under the thresholds `site_args`, reporting each file to
/// `progress`.
fn write_fileset(
    prep: &Prepared,
    args: &QcArgs,
    out_dir: &str,
    site_args: &SiteFilterArgs,
    figures: &[Figure],
    progress: &Progress,
) -> anyhow::Result<()> {
    let Prepared {
        layout,
        cells,
        decisions,
        tables,
        site_matrices,
        site_cells,
    } = prep;
    let other_matrices: Vec<&MatrixFile> = {
        let repooled: FxHashSet<(&str, &str)> = site_matrices
            .iter()
            .map(|sm| (&*sm.m.batch, sm.m.site_modality().unwrap_or_default()))
            .collect();
        layout
            .matrices
            .iter()
            .filter(|m| m.site_modality().is_none() && !repooled.contains(&(&*m.batch, &*m.kind)))
            .collect()
    };
    let n_site_tables = SITE_MODALITIES
        .iter()
        .filter(|m| tables.contains_key(**m))
        .count();
    // One step per file group below, in order.
    let n_steps = decisions.len()
        + n_site_tables
        + 2 * site_matrices.len()
        + other_matrices.len()
        + figures.len()
        + 2;
    progress.plan(n_steps);
    let mut summary: Vec<SummaryRow> = Vec::new();

    std::fs::create_dir_all(out_dir)?;
    for d in decisions {
        progress.next(format!("{}: cell verdicts", d.batch));
        d.write(out_dir)?;
    }
    let mut kept_sites: FxHashMap<Box<str>, FxHashSet<Box<str>>> = FxHashMap::default();
    for modality in SITE_MODALITIES {
        let Some(t) = tables.get(*modality) else {
            continue;
        };
        progress.next(format!("{modality}_sites.parquet"));
        let n_cells = site_cells.get(*modality).map(|acc| t.cells_per_site(acc));
        let reasons = site_args.reasons(t, n_cells.as_deref());
        let kept: FxHashSet<Box<str>> = t
            .key
            .iter()
            .zip(&reasons)
            .filter(|(_, r)| r.is_none())
            .map(|(k, _)| k.clone())
            .collect();
        let (n_kept, n_dropped) = write_site_tables(
            t,
            &reasons,
            &format!("{out_dir}/{modality}_sites.parquet"),
            &format!("{out_dir}/{modality}_sites_dropped.parquet"),
        )?;
        let mut by_reason: FxHashMap<&str, usize> = FxHashMap::default();
        for r in reasons.iter().flatten() {
            *by_reason.entry(r.label()).or_insert(0) += 1;
        }
        let mut by_reason: Vec<_> = by_reason.into_iter().collect();
        by_reason.sort();
        info!(
            "{modality}: kept {n_kept} sites, dropped {n_dropped} {by_reason:?}; cells per site {}",
            if n_cells.is_some() {
                "from the _site matrices"
            } else {
                "unavailable (no _site matrix), cells check skipped"
            }
        );
        kept_sites.insert((*modality).into(), kept);
    }

    // 3 + 4. site matrices, then the re-pooled gene-level matrices.
    for sm in site_matrices {
        let modality = sm.m.site_modality().unwrap_or_default();
        let rows = site_matrix_rows(&sm.opened.row_names, kept_sites.get(modality));
        let stem = sm.m.stem();
        progress.next(&*stem);
        let w = write_subset(
            sm.opened.data.as_ref(),
            &sm.cols,
            &rows,
            &sm.opened.row_names,
            &sm.opened.col_names,
            &out_spec(&sm.m, out_dir, &stem, args.no_zip),
        )?;
        record(&mut summary, &stem, sm.opened.before, w);

        let gene_stem = format!("{}_{}", sm.m.batch, modality);
        progress.next(format!("{gene_stem} (re-pooled)"));
        let w = repool_gene_level(
            sm.opened.data.as_ref(),
            &sm.cols,
            &rows,
            &sm.opened.row_names,
            &sm.opened.col_names,
            args.row_nnz_cutoff,
            &out_spec(&sm.m, out_dir, &gene_stem, args.no_zip),
        )?;
        let before = layout
            .matrix(&sm.m.batch, modality)
            .and_then(|gm| open_matrix(&gm.path).ok())
            .map(|d| shape(d.as_ref()))
            .unwrap_or((0, 0, 0));
        record(&mut summary, &gene_stem, before, w);
    }

    // 5. everything else: cells + row nnz.
    for m in other_matrices {
        let stem = m.stem();
        progress.next(&*stem);
        let opened = open_with_names(m)?;
        let cols = select_columns(opened.data.as_ref(), &opened.col_names, cells.get(&m.batch));
        let row_nnz = row_nnz_over_columns(opened.data.as_ref(), &cols)?;
        let cutoff = row_cutoff(args, &row_nnz, &format!("{stem} rows"));
        let rows: Vec<usize> = (0..row_nnz.len())
            .filter(|&r| row_nnz[r] >= cutoff)
            .collect();
        let w = write_subset(
            opened.data.as_ref(),
            &cols,
            &rows,
            &opened.row_names,
            &opened.col_names,
            &out_spec(m, out_dir, &stem, args.no_zip),
        )?;
        record(&mut summary, &stem, opened.before, w);
    }

    // 6. copy-through, except the per-batch cell lists `qc` rewrote.
    progress.next("other files");
    for f in &layout.other_files {
        let name = file_name(f);
        if name
            .strip_suffix("_cells.tsv.gz")
            .is_some_and(|b| cells.contains_key(b))
        {
            continue;
        }
        std::fs::copy(f.as_ref(), format!("{out_dir}/{name}"))?;
    }

    progress.next("qc_summary.tsv");
    let mut lines: Vec<Box<str>> = vec![
        "#file\trows_before\tcols_before\tnnz_before\trows_after\tcols_after\tnnz_after".into(),
    ];
    for s in &summary {
        lines.push(
            format!(
                "{}\t{}\t{}\t{}\t{}\t{}\t{}",
                s.file, s.before.0, s.before.1, s.before.2, s.after.0, s.after.1, s.after.2
            )
            .into(),
        );
    }
    write_lines(&lines, &format!("{out_dir}/qc_summary.tsv"))?;
    // 7. the view's figures, as it showed them when the cut was applied.
    //    The fileset is complete without them, so a failed one only warns.
    for (stem, e) in save_figures(out_dir, figures, progress) {
        log::warn!("qc_plots/{stem}: not saved: {e}");
    }
    debug_assert_eq!(
        progress.started(),
        n_steps,
        "the progress plan missed a step"
    );
    info!("done: {out_dir}");
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(cmd: &str) -> QcArgs {
        #[derive(clap::Parser)]
        struct Cli {
            #[command(flatten)]
            qc: QcArgs,
        }
        let words = cmd.split(' ').skip(2); // `faba qc`
        <Cli as clap::Parser>::parse_from(std::iter::once("faba").chain(words)).qc
    }

    /// Each step of the 2D view timed on a real faba directory, named by
    /// `FABA_QC_TIMING_DIR` (read only; what it writes goes to a temporary
    /// folder): `cargo test --release qc_timing -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn qc_timing() {
        let Ok(input) = std::env::var("FABA_QC_TIMING_DIR") else {
            eprintln!("FABA_QC_TIMING_DIR is not set");
            return;
        };
        let args = parse(&format!("faba qc {input} --batch-process -o unused"));
        let t = std::time::Instant::now();
        let prep = prepare(&args, &input).unwrap();
        let read = t.elapsed().as_secs_f64();
        let out = tempfile::tempdir().unwrap();
        let cells = cells_of(&prep.decisions);
        let steps = crate::qc::site_tui::time_steps(
            &prep.tables,
            &prep.site_cells,
            &args.site,
            cells,
            out.path(),
        );
        eprintln!("{read:>9.3} s  reading the site tables, the _site matrices and the cell QC");
        for (what, s) in steps {
            eprintln!("{s:>9.3} s  {what}");
        }
    }

    /// One modality's sampled sites, as the snapshot keeps them.
    #[derive(serde::Serialize, serde::Deserialize)]
    /// (Reals as their bits: JSON has no NaN or infinity.)
    struct Snap {
        modality: String,
        pv: Vec<u32>,
        coverage: Vec<u64>,
        converted: Vec<u64>,
        control_coverage: Vec<u64>,
        edit_ratio: Vec<u32>,
        fold: Vec<u32>,
        raw_log_odds: Vec<u64>,
        gene_id: Vec<u32>,
        n_genes: usize,
        n_cells: Option<Vec<usize>>,
    }

    /// The figures and the page redrawn from a snapshot in
    /// `FABA_QC_SNAPSHOT` (a JSON file), into `FABA_QC_REDRAW_OUT`; the
    /// snapshot made first, once, from the faba directory in
    /// `FABA_QC_TIMING_DIR` (read only), sampling at most 60,000 sites per
    /// modality: `cargo test qc_redraw -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn qc_redraw() {
        let (Ok(snap), Ok(out)) = (
            std::env::var("FABA_QC_SNAPSHOT"),
            std::env::var("FABA_QC_REDRAW_OUT"),
        ) else {
            eprintln!("FABA_QC_SNAPSHOT and FABA_QC_REDRAW_OUT are not set");
            return;
        };
        type Saved = (Vec<Snap>, Vec<crate::qc::cells_plot::Cells>);
        if !std::path::Path::new(&snap).exists() {
            let input = std::env::var("FABA_QC_TIMING_DIR").expect("FABA_QC_TIMING_DIR");
            let args = parse(&format!("faba qc {input} --batch-process -o unused"));
            let prep = prepare(&args, &input).unwrap();
            let mut snaps = Vec::new();
            for (m, t) in &prep.tables {
                let n_cells = prep.site_cells.get(m).map(|acc| t.cells_per_site(acc));
                let step = t.len().div_ceil(60_000).max(1);
                let pick = |i: &usize| i.is_multiple_of(step);
                let idx: Vec<usize> = (0..t.len()).filter(pick).collect();
                let take = |v: &[f32]| idx.iter().map(|&i| v[i].to_bits()).collect::<Vec<_>>();
                let takeu = |v: &[u64]| idx.iter().map(|&i| v[i]).collect::<Vec<_>>();
                snaps.push(Snap {
                    modality: m.to_string(),
                    pv: take(&t.pv),
                    coverage: takeu(&t.coverage),
                    converted: takeu(&t.converted),
                    control_coverage: takeu(&t.control_coverage),
                    edit_ratio: take(&t.edit_ratio),
                    fold: take(&t.fold),
                    raw_log_odds: idx.iter().map(|&i| t.raw_log_odds[i].to_bits()).collect(),
                    gene_id: idx.iter().map(|&i| t.gene_id[i]).collect(),
                    n_genes: t.n_genes,
                    n_cells: n_cells.map(|n| idx.iter().map(|&i| n[i]).collect()),
                });
            }
            let saved: Saved = (snaps, cells_of(&prep.decisions).to_vec());
            std::fs::write(&snap, serde_json::to_string(&saved).unwrap()).unwrap();
        }
        let (snaps, cells): Saved =
            serde_json::from_str(&std::fs::read_to_string(&snap).unwrap()).unwrap();
        let tables: Vec<(SiteTable, Option<Vec<usize>>)> = (snaps.into_iter())
            .map(|s| {
                let t = SiteTable {
                    modality: s.modality.into(),
                    batch: arrow::record_batch::RecordBatch::new_empty(std::sync::Arc::new(
                        arrow::datatypes::Schema::empty(),
                    )),
                    key: (0..s.pv.len()).map(|i| format!("s{i}").into()).collect(),
                    gene_id: s.gene_id,
                    n_genes: s.n_genes,
                    pv: s.pv.iter().map(|&b| f32::from_bits(b)).collect(),
                    coverage: s.coverage,
                    converted: s.converted,
                    control_coverage: s.control_coverage,
                    edit_ratio: s.edit_ratio.iter().map(|&b| f32::from_bits(b)).collect(),
                    fold: s.fold.iter().map(|&b| f32::from_bits(b)).collect(),
                    raw_log_odds: s.raw_log_odds.iter().map(|&b| f64::from_bits(b)).collect(),
                };
                (t, s.n_cells)
            })
            .collect();
        let mut tables = tables;
        tables.sort_by_key(|(t, _)| SITE_MODALITIES.iter().position(|m| **m == *t.modality));
        std::fs::create_dir_all(&out).unwrap();
        let t = std::time::Instant::now();
        let f = SiteFilterArgs::default_values();
        let written = crate::qc::site_tui::redraw(&tables, &f, cells, std::path::Path::new(&out));
        eprintln!("redrawn in {:.3} s:", t.elapsed().as_secs_f64());
        written.iter().for_each(|p| eprintln!("  {}", p.display()));
    }

    /// `--html` goes only to a file that is not there, with --batch-process
    /// or --dry-run: an existing one is refused before anything is read
    /// (the input here does not exist), naming it and asking for another
    /// path, and left as it was.
    #[test]
    fn html_refuses_a_file_that_exists_before_any_work() {
        let dir = tempfile::tempdir().unwrap();
        let taken = dir.path().join("qc.html");
        std::fs::write(&taken, "mine").unwrap();
        let t = taken.display();
        for mode in ["--batch-process -o out", "--dry-run"] {
            let args = parse(&format!("faba qc /no/such/input {mode} --html {t}"));
            let err = run_qc(&args).unwrap_err().to_string();
            assert!(
                err.contains(&t.to_string()) && err.contains("another path"),
                "{err}"
            );
        }
        assert_eq!(std::fs::read_to_string(&taken).unwrap(), "mine");
        let free = dir.path().join("new.html");
        let view = parse(&format!("faba qc /no/such/input --html {}", free.display()));
        assert!(html_problem(&view)
            .unwrap_err()
            .to_string()
            .contains("Ctrl+E"));
        let batch = parse(&format!("faba qc in --dry-run --html {}", free.display()));
        assert!(html_problem(&batch).is_ok());
        // A folder not there yet, inside a cut's output folder or not, is
        // made when the page is written.
        let out = dir.path().join("cut");
        let o = out.display();
        for page in [format!("{o}/qc.html"), format!("{o}/pages/qc.html")] {
            let inside = parse(&format!("faba qc in --batch-process -o {o} --html {page}"));
            assert!(html_problem(&inside).is_ok(), "{page}");
        }
    }

    /// A dry run's figures go under -o DIR, new or empty: a DIR with files
    /// is refused before anything is read, and one is needed.
    #[test]
    fn dry_run_plots_refuse_a_folder_with_files() {
        let dir = tempfile::tempdir().unwrap();
        std::fs::write(dir.path().join("mine.txt"), "mine").unwrap();
        let d = dir.path().display();
        let args = parse(&format!("faba qc /no/such/input --dry-run --plots -o {d}"));
        let err = run_qc(&args).unwrap_err().to_string();
        assert!(err.contains("already contains files"), "{err}");
        let none = parse("faba qc /no/such/input --dry-run --plots");
        assert!(dry_plots_problem(&none)
            .unwrap_err()
            .to_string()
            .contains("-o DIR"));
        let empty = tempfile::tempdir().unwrap();
        let fresh = parse(&format!(
            "faba qc in --dry-run --plots -o {}",
            empty.path().display()
        ));
        assert!(dry_plots_problem(&fresh).is_ok());
    }

    /// The figures' page is checked before anything is read: narrower than
    /// a knob figure's table is refused, naming --panel-only, which draws
    /// the plot alone at any size.
    #[test]
    fn a_narrow_page_is_refused_unless_the_panel_is_drawn_alone() {
        let narrow = parse("faba qc in --batch-process -o out --plots --fig-width 2.5");
        let e = figure_sheet(&narrow).unwrap_err().to_string();
        assert!(e.contains("--panel-only"), "{e}");
        let alone = parse(
            "faba qc in --batch-process -o out --plots --fig-width 2.5 --panel-only --font-size 6",
        );
        let sheet = figure_sheet(&alone).unwrap();
        assert_eq!(sheet.page, (Some(180.0), None));
        assert!(sheet.look.panel_only && sheet.look.font_size == Some(6.0));
        let bad = parse("faba qc in --batch-process -o out --plots --font-size 0");
        assert!(figure_sheet(&bad).is_err());
    }

    #[test]
    fn the_printed_command_reruns_the_same_cut_in_batch() {
        let args = parse("faba qc in -o out -r 3 --no-cell-qc --auto-cutoff --block-size 64");
        let mut site = args.site.clone();
        site.site_max_pv = 0.01;
        site.site_min_cells = 4;
        let again = parse(&batch_command(
            &args,
            "in",
            "out",
            &site,
            args.gff.as_deref(),
        ));
        assert!(again.batch_process);
        assert_eq!(
            (again.input_dir.as_deref(), again.output.as_deref()),
            (Some("in"), Some("out"))
        );
        assert_eq!(qc_flags(&again.site), qc_flags(&site));
        assert_eq!(
            serde_json::to_value(QcArgs {
                batch_process: false,
                site: site.clone(),
                ..again
            })
            .unwrap(),
            serde_json::to_value(QcArgs { site, ..args }).unwrap()
        );
    }
}
