//! `faba cluster`: cells clustered by their own site reads.
//!
//! A blockmodel over (cell, gene) links of site reads (see [`site_vi`]), or
//! over a count matrix's cells and genes (see [`count`]), each fitted by
//! variational inference from several k-means starts, the best by ELBO
//! kept. Cells from inputs outside `--fit` are placed against the fit
//! afterwards, as a control.

use crate::run_record::RunRecord;
use log::{info, warn};
use std::borrow::Cow;
use std::io::IsTerminal;
use std::path::Path;

mod args;
pub(crate) mod batches;
mod browse;
mod command;
mod count;
mod links;
mod math;
mod means;
mod progress;
mod read;
/// Each model's fit, from the read and filtered input to its tables.
mod site_vi;
mod run {
    pub(super) mod count;
    pub(super) mod sites;
    pub(super) mod vi;
}
mod settings;
mod vi;
mod weights;
mod within;
mod write;

pub use args::ClusterArgs;
use args::{check, inputs, output_dir, with_inputs};
pub(crate) use args::{refused, REMOVED_SAID};
#[cfg(feature = "mcp")]
pub(crate) use browse::opens;
pub use browse::{fitting_hint, run_modules, ModulesArgs};
use command::command;

#[cfg(test)]
mod tests;

/// The run record of a site fit of modality `m`, `{m}_blockmodel.run.json`
/// (the command saved beside it, `.cmd.sh` for `.run.json`).
fn record_name(m: &str) -> String {
    format!("{m}_blockmodel.run.json")
}

/// Where faba before 0.19 kept the record of a mixture fit (`--model crp`),
/// and before that of every site fit: `{m}_cluster.run.json`.
fn older_record_name(m: &str) -> String {
    format!("{m}_cluster.run.json")
}

/// Bytes the fit (not the reading) will hold at its peak, roughly: the
/// reads (raw, a copy of the fitted ones when some are left out, and
/// weighted); the links' copy of the reads and each link's gene; per
/// restart each cell's role weights and level, each gene's role weights and
/// the blocks'; the best restart's kept beside the one running; and the
/// tables: each cell's and gene's links and each site's reads per block,
/// with the statistics' fixed parts over the genes.
fn memory_estimate(
    entries: usize,
    copied: bool,
    (sites, cells, links): (usize, usize, usize),
    args: &ClusterArgs,
) -> usize {
    let reads = entries * (12 + 16 + if copied { 12 } else { 0 });
    let s = args::blockmodel(args);
    let k = s.cell_roles * s.gene_roles;
    let genes = sites.min(links);
    let restart = cells * (s.cell_roles + 2) * 8 + genes * s.gene_roles * 8;
    // The start's profiles (each link's residual and its f32 row entry)
    // shared by the restarts' starts; the statistics' 64 fixed parts.
    reads
        + entries * 20
        + links * (8 + 4)
        + links * (8 + 4)
        + 2 * restart
        + (cells + genes) * k * 8
        + sites * k * 16
        + 64 * (genes * s.gene_roles + 2 * k + 1) * 8
}

/// Whether the restarts are drawn in the live view as they run: on a
/// terminal, unless --batch-process says a bar.
fn live(args: &ClusterArgs) -> bool {
    !args.batch_process
        && std::io::stderr().is_terminal()
        && data_beans::interactive::tui_available()
}

/// `work` run counting on a progress of its own, which comes back with its
/// result: on a terminal in the live view under `title` (see
/// [`progress::view`]), where `q` cancels and Ctrl+C stops the restarts
/// early, the view staying up after a fit until a key when `stay`; else
/// with a bar on stderr.
fn watched<T: Send>(
    args: &ClusterArgs,
    title: &str,
    stay: bool,
    work: impl FnOnce(&progress::Progress) -> anyhow::Result<T> + Send,
) -> (progress::Progress, anyhow::Result<T>) {
    if live(args) {
        let progress = progress::Progress::interruptible(args.chains);
        let out = progress::view::watch(title, &progress, stay, work);
        (progress, out)
    } else {
        let progress = progress::Progress::on_stderr(args.chains);
        let out = work(&progress);
        (progress, out)
    }
}

/// What a run that ended in `error` comes to: a cancel asked for in the
/// live view is not a failure (nothing written, said so, no modality);
/// any other error as it is.
fn unless_cancelled(args: &ClusterArgs, error: anyhow::Error) -> anyhow::Result<Option<String>> {
    if live(args) && error.is::<progress::Cancelled>() {
        eprintln!("faba cluster: {error}");
        return Ok(None);
    }
    Err(error)
}

/// Cluster as `args` say; the modality fitted, none when no inputs were
/// chosen. On a terminal the restarts are drawn in a view as they run, else
/// on a bar on stderr (see [`watched`]). A fit ends with its restarts'
/// ELBO traces and advice on stderr (see [`progress::spark`]). Flags of the
/// MCMC fits faba no longer has are refused before anything is read or
/// written (see [`args::refuse_removed`]).
/// With `--within`, one fit per cell module of a count fit (see [`within`]).
pub fn run_cluster(args: &ClusterArgs) -> anyhow::Result<Option<String>> {
    args::refuse_removed(args)?;
    if args.within.is_some() {
        return within::run(args);
    }
    let offered = if args.is_count() {
        (
            "count matrices",
            crate::qc::layout::is_count_matrix_name as fn(&str) -> bool,
        )
    } else {
        (
            "site matrices",
            crate::qc::layout::is_site_matrix_name as fn(&str) -> bool,
        )
    };
    // The inputs are asked for first, in a browser of their own.
    let Some(data_files) = with_inputs(&args.data_files, "cluster", false, offered)? else {
        return Ok(None);
    };
    let args = &ClusterArgs {
        data_files,
        ..args.clone()
    };
    let model = if args.is_count() {
        "count blockmodel"
    } else {
        "site blockmodel"
    };
    let dir = output_dir(args.output.as_deref(), &args.data_files);
    let title = format!("fitting the {model} in {}", dir.display());
    let (progress, fitted) = watched(args, &title, true, |p| run_cluster_watched(args, p));
    match fitted {
        Ok(modality) => {
            progress::spark::tell(&progress);
            Ok(Some(modality))
        }
        Err(e) => unless_cancelled(args, e),
    }
}

/// [`run_cluster`] of `args`, their inputs named, the restarts counting
/// their steps on `progress`; the modality fitted. Cancelled, it writes nothing,
/// not even its record.
fn run_cluster_watched(
    args: &ClusterArgs,
    progress: &progress::Progress,
) -> anyhow::Result<String> {
    let dir = output_dir(args.output.as_deref(), &args.data_files);
    let record = RunRecord::start("cluster", &dir.to_string_lossy())
        .inputs("matrices", &args.data_files)
        .options(args);
    recorded(
        record,
        |record| cluster(args, &dir, record, progress),
        |written| save_command(args, &dir, &written.record),
    )
    .map(|written| written.modality)
}

/// `fit` run under `record`, which is then written under the name the fit
/// gives it, `then` first (so what it writes is among the record's
/// outputs); a fit that fails is recorded as failed. Cancelled, nothing is
/// written, not even the record.
fn recorded(
    mut record: RunRecord,
    fit: impl FnOnce(&mut RunRecord) -> anyhow::Result<Written>,
    then: impl FnOnce(&Written),
) -> anyhow::Result<Written> {
    let outcome = match fit(&mut record) {
        Err(e) if e.is::<progress::Cancelled>() => return Err(e),
        outcome => outcome,
    };
    if let Ok(written) = &outcome {
        record = record.file_name(&written.record);
        then(written);
    }
    record.finish(&outcome);
    outcome
}

/// The command of `args`, to read, edit and run again, saved in `dir`
/// beside the run record named `record` (`.cmd.sh` for `.run.json`).
fn save_command(args: &ClusterArgs, dir: &Path, record: &str) {
    let name = record.replace(".run.json", ".cmd.sh");
    let script = crate::script::Script {
        name: &name,
        made_by: "faba cluster",
        guard: None,
        cd: false,
        head: 1,
        replace: true,
    };
    // The fit is done and recorded whether or not the script is written.
    if let Err(e) = crate::script::write(dir, &script, &command(args)) {
        warn!("the command was not saved: {e:#}");
    }
}

/// What a fit wrote: the modality fitted (lowercased, as its tables are
/// named) and its run record's file name.
struct Written {
    modality: String,
    record: String,
}

/// Fit and write, the restarts' steps counted on `progress`; what was
/// written.
fn cluster(
    args: &ClusterArgs,
    dir: &Path,
    record: &mut RunRecord,
    progress: &progress::Progress,
) -> anyhow::Result<Written> {
    check(args)?;
    if args.is_count() {
        return run::count::count_fit(args, dir, record, progress);
    }
    let all = read::read_inputs(&inputs(args)?, progress)?;
    site_fit(args, dir, record, progress, all)
}

/// The site fit of `args` over the cells and reads of `all`: filtered,
/// fitted, summarised and written to `dir`; what was written.
fn site_fit(
    args: &ClusterArgs,
    dir: &Path,
    record: &mut RunRecord,
    progress: &progress::Progress,
    all: read::SiteData,
) -> anyhow::Result<Written> {
    let (data, kept) = read::filter(&all, args.min_cell_sites, args.min_site_cells);
    progress.check()?;
    // Only the input cells are needed from here: let their reads go.
    let read_sites = all.sites.len();
    let read::SiteData {
        cells: all_cells, ..
    } = all;
    let (fitted, controls): (Vec<usize>, Vec<usize>) =
        (0..data.cells.len()).partition(|&c| data.cells[c].fit);
    anyhow::ensure!(
        fitted.len() >= 2 && !data.sites.is_empty(),
        "{} cells and {} sites pass the filters: too few to cluster",
        fitted.len(),
        data.sites.len()
    );
    info!(
        "{}: {} of {} cells kept ({} fitted), {} of {} sites",
        data.modality,
        data.cells.len(),
        all_cells.len(),
        fitted.len(),
        data.sites.len(),
        read_sites
    );
    // The fitted cells' reads: all of them, or a copy of those fitted.
    let fit_reads: Cow<read::Reads> = if fitted.len() < data.cells.len() {
        Cow::Owned(data.reads.subset(&fitted))
    } else {
        Cow::Borrowed(&data.reads)
    };
    let links = links::Links::count(&fit_reads, &run::sites::number_genes(&data.sites).0);
    let memory = memory_estimate(
        fit_reads.site.len(),
        matches!(fit_reads, Cow::Owned(_)),
        (data.sites.len(), fitted.len(), links),
        args,
    );
    info!("about {} MB at the peak of the fit", memory >> 20);
    progress.check()?;

    let (weighted, batch_stats) = weigh_batches(args, &data, &fitted, &fit_reads);
    progress.check()?;

    record.set_result(
        "cells",
        &serde_json::json!({
            "read": all_cells.len(),
            "kept": data.cells.len(),
            "fitted": fitted.len(),
            "assigned": data.cells.len() - fitted.len(),
        }),
    );
    record.set_result(
        "sites",
        &serde_json::json!({"read": read_sites, "kept": data.sites.len()}),
    );
    record.set_result("batch_weights", &batch_stats);
    record.set_result("fit_memory_estimate_bytes", &memory);
    let fit = FitData {
        data: &data,
        fitted: &fitted,
        weighted: &weighted,
        all_cells: &all_cells,
        kept: &kept,
        controls: &controls,
        progress,
    };
    run::sites::blocks_fit(args, dir, record, fit)
}

/// Batches taken out of the fitted reads by importance weights, when there
/// are several and the correction is on: the reads weighted, and what the
/// weights did, for the record.
fn weigh_batches<'r>(
    args: &ClusterArgs,
    data: &read::SiteData,
    fitted: &[usize],
    fit_reads: &'r read::Reads,
) -> (read::Weighted<'r>, Option<serde_json::Value>) {
    let batches: Vec<&str> = {
        let mut b: Vec<&str> = fitted
            .iter()
            .map(|&c| data.cells[c].batch.as_ref())
            .collect();
        b.sort_unstable();
        b.dedup();
        b
    };
    if !args.no_batch_correction && batches.len() > 1 {
        let batch_of: Vec<u16> = fitted
            .iter()
            .map(|&c| {
                batches
                    .binary_search(&data.cells[c].batch.as_ref())
                    .unwrap() as u16
            })
            .collect();
        let w = weights::batch_weights(fit_reads, &batch_of, batches.len(), data.sites.len());
        info!(
            "importance weights take out {} batches: effective share of reads {:?}, {:.2}% capped",
            batches.len(),
            w.effective_share
                .iter()
                .map(|e| format!("{e:.2}"))
                .collect::<Vec<_>>(),
            100.0 * w.capped_share
        );
        let stats = serde_json::json!({
            "batches": batches,
            "effective_share": w.effective_share,
            "capped_share": w.capped_share,
            "max_weight": w.max_weight,
        });
        (
            read::Weighted::new(fit_reads, Some((&w.on, &w.off))),
            Some(stats),
        )
    } else {
        (read::Weighted::new(fit_reads, None), None)
    }
}

/// What the site fit needs from the read and filtered input: the kept
/// cells' data, those fitted and those left out of the fit (the controls,
/// both in order), the fitted reads as weighted, every input
/// cell, and each kept cell's index among them.
struct FitData<'a> {
    data: &'a read::SiteData,
    fitted: &'a [usize],
    weighted: &'a read::Weighted<'a>,
    all_cells: &'a [read::Cell],
    kept: &'a [usize],
    controls: &'a [usize],
    /// Where the restarts count their steps.
    progress: &'a progress::Progress,
}
