//! `faba cluster`: cells clustered by their own site reads.
//!
//! A blockmodel over (cell, gene) links (see [`unified`]), its link and
//! blocks presets, or a Chinese-restaurant-process mixture over cells (see
//! [`crp`]), run as several chains and summarised without a cells x cells
//! matrix (see [`summary`]). Cells from inputs outside `--fit` are assigned
//! to the fitted roles or clusters afterwards, as a control.

use crate::run_record::RunRecord;
use legume_numeric::mcmc::engine::McmcConfig;
use log::{info, warn};
use std::borrow::Cow;
use std::io::IsTerminal;
use std::path::Path;

mod args;
mod browse;
mod command;
mod count;
mod crp;
mod embed;
mod fast;
mod fit;
mod progress;
mod read;
/// Each model's fit, from the read and filtered input to its tables.
mod run {
    pub(super) mod blockmodel;
    pub(super) mod count;
    pub(super) mod crp;
}
mod summary;
mod unified;
mod within;
mod write;

use args::{check, inputs, output_dir, with_inputs};
pub use args::{ClusterArgs, ModelKind};
#[cfg(feature = "mcp")]
pub(crate) use browse::opens;
pub use browse::{fitting_hint, run_modules, ModulesArgs};
use command::command;

#[cfg(test)]
mod tests;

/// The run record of `model`'s fit of modality `m`: the mixture's
/// `{m}_cluster.run.json`, the blockmodel's `{m}_blockmodel.run.json`, each
/// its own so that neither replaces the other's (the command saved beside
/// each, `.cmd.sh` for `.run.json`).
fn record_name(m: &str, model: ModelKind) -> String {
    match model {
        ModelKind::Crp => format!("{m}_cluster.run.json"),
        ModelKind::Link | ModelKind::Blocks => format!("{m}_blockmodel.run.json"),
    }
}

/// Split R-hat above this, 1.05, counts as not converged.
const RHAT_OK: f64 = 1.05;

/// Bytes the fit (not the reading) will hold at its peak, roughly: the
/// reads (raw, a copy of the fitted ones when some are left out, and
/// weighted), then the model's own per-chain state. `links`: the
/// blockmodel's (cell, gene) links (unused by the mixture).
fn memory_estimate(
    entries: usize,
    copied: bool,
    (sites, cells, links): (usize, usize, usize),
    kept: usize,
    args: &ClusterArgs,
) -> usize {
    let reads = entries * (12 + 16 + if copied { 12 } else { 0 });
    let (shared, chain) = match args.model {
        // The picture chains start from; per chain two f64 per site for
        // each live cluster (a guess of 32), split-merge scratch, labels
        // and kept draws.
        ModelKind::Crp => (
            sites * 32 * 4 + cells * 32 * 4,
            sites * 16 * 32 + sites * 33 + cells * 4 + kept * cells * 2,
        ),
        // The links: their own copy of the weighted reads, each link's
        // cell, entries and gene, its unconverted reads and their baselines
        // summed, and the genes' links; the sites' baselines with their
        // exponentials and whether each is a Gamma site. Per chain: role
        // counts (u32) and their sums (f64) for cells and genes (at most
        // one gene per site), two labels per link, at finite strength
        // reads, profile, its `ln B`, the `ln Gamma` cache and holders per
        // site and block (with the shifts' list of each block's sites), at
        // infinite strength the shifts' list of each block's links; and the
        // tally's links per cell and gene and reads per site in each block.
        ModelKind::Link | ModelKind::Blocks => {
            let s = args::blockmodel(args);
            let k = s.cell_roles * s.gene_roles;
            let blocks = sites * k;
            let per_block = if s.strength.is_finite() {
                blocks * (16 + 16 + 8 + 24 + 4 + 8)
            } else {
                links * 8
            };
            (
                entries * 20 + links * (12 + 16 + 12) + sites * (8 + 8 + 1),
                (cells * s.cell_roles + sites * s.gene_roles) * 12
                    + links * 4
                    + per_block
                    + (cells + sites) * k * 8
                    + blocks * 16,
            )
        }
    };
    reads + shared + args.chains * chain
}

/// Whether the chains are drawn in the live view as they run: on a
/// terminal, unless --batch-process says a bar.
fn live(args: &ClusterArgs) -> bool {
    !args.batch_process
        && std::io::stderr().is_terminal()
        && data_beans::interactive::tui_available()
}

/// `work` run counting on a progress of its own, which comes back with its
/// result: on a terminal in the live view under `title` (see
/// [`progress::view`]), where `q` cancels and Ctrl+C stops the chains
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
/// chosen. On a terminal the chains are drawn in a view as they run, else
/// on a bar on stderr (see [`watched`]). A blockmodel fit ends with its
/// chains' traces, R-hat and advice on stderr (see [`progress::spark`]).
/// With `--within`, one fit per cell module of a count fit (see [`within`]).
pub fn run_cluster(args: &ClusterArgs) -> anyhow::Result<Option<String>> {
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
        "count blockmodel".to_string()
    } else {
        format!("{} model", args.model.name())
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

/// [`run_cluster`] of `args`, their inputs named, the chains counting their
/// sweeps on `progress`; the modality fitted. Cancelled, it writes nothing,
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

/// The chains' schedule for `args`: `--burn-in` sweeps, then the
/// `(sweeps - burn-in) / thin` kept draws every `thin`-th sweep, the
/// engine stopping on the last of them. A `thin` that does not divide the
/// sweeps after the burn-in leaves the remainder unrun; the progress shows
/// the sweeps run ([`progress::sweeps`]), not `--sweeps`.
fn mcmc_config(args: &ClusterArgs) -> McmcConfig {
    McmcConfig {
        n_samples: kept_draws(args).max(1),
        warmup: args.burn_in,
        thin: args.thin,
        seed: args.seed,
    }
}

/// The draws `args` keep: `(sweeps - burn-in) / thin`, a `thin` of 0 read
/// as 1; none when the sweeps end with the burn-in.
fn kept_draws(args: &ClusterArgs) -> usize {
    (args.sweeps - args.burn_in) / args.thin.max(1)
}

/// What a fit wrote: the modality fitted (lowercased, as its tables are
/// named) and its run record's file name.
struct Written {
    modality: String,
    record: String,
}

/// Fit, summarise and write, the sweeps counted on `progress`; what was
/// written.
fn cluster(
    args: &ClusterArgs,
    dir: &Path,
    record: &mut RunRecord,
    progress: &progress::Progress,
) -> anyhow::Result<Written> {
    check(args)?;
    if args.is_count() {
        return run::count::count_fit(args, &mcmc_config(args), dir, record, progress);
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
    // Entries held out of the fit, the same for either model.
    let (fit_reads, held): (Cow<read::Reads>, Option<read::Reads>) = if args.holdout > 0.0 {
        let (fit, held) = read::hold_out(&fit_reads, args.holdout, args.seed);
        (Cow::Owned(fit), Some(held))
    } else {
        (fit_reads, None)
    };
    let links = match args.model {
        ModelKind::Crp => 0,
        ModelKind::Link | ModelKind::Blocks => {
            unified::Links::count(&fit_reads, &run::blockmodel::number_genes(&data.sites).0)
        }
    };
    let memory = memory_estimate(
        fit_reads.site.len(),
        matches!(fit_reads, Cow::Owned(_)),
        (data.sites.len(), fitted.len(), links),
        kept_draws(args),
        args,
    );
    info!("about {} MB at the peak of the fit", memory >> 20);
    progress.check()?;

    let (weighted, batch_stats) = weigh_batches(args, &data, &fitted, &fit_reads);
    progress.check()?;

    let config = mcmc_config(args);
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
        held: held.as_ref(),
        all_cells: &all_cells,
        kept: &kept,
        controls: &controls,
        progress,
    };
    match args.model {
        ModelKind::Crp => run::crp::crp_fit(args, &config, dir, record, fit, &fit_reads),
        ModelKind::Link | ModelKind::Blocks => {
            run::blockmodel::blockmodel_fit(args, &config, dir, record, fit)
        }
    }
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
        let w = fit::batch_weights(fit_reads, &batch_of, batches.len(), data.sites.len());
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

/// What each model needs from the read and filtered input: the kept
/// cells' data, those fitted and those left out of the fit (the controls,
/// both in order), the fitted reads as weighted and held out, every input
/// cell, and each kept cell's index among them.
struct FitData<'a> {
    data: &'a read::SiteData,
    fitted: &'a [usize],
    weighted: &'a read::Weighted<'a>,
    held: Option<&'a read::Reads>,
    all_cells: &'a [read::Cell],
    kept: &'a [usize],
    controls: &'a [usize],
    /// Where the chains count their sweeps.
    progress: &'a progress::Progress,
}
