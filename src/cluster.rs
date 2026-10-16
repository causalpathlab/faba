//! `faba cluster`: cells clustered by their own site reads.
//!
//! A Chinese-restaurant-process mixture over cells (see [`crp`]), run as
//! several chains and summarised without a cells x cells matrix (see
//! [`summary`]). Cells from inputs outside `--fit` are assigned to the fitted
//! clusters afterwards, as a control.

use crate::run_record::RunRecord;
use clap::Args;
use legume_numeric::mcmc::engine::{run_mcmc_parallel, McmcConfig};
use log::{info, warn};
use rustc_hash::{FxHashMap, FxHashSet};
use std::borrow::Cow;
use std::path::{Path, PathBuf};

mod crp;
mod embed;
mod fit;
mod read;
mod summary;
mod write;

#[cfg(test)]
mod tests;

/// Split R-hat above this is reported as not converged.
const RHAT_OK: f64 = 1.05;

#[derive(Args, Debug, Clone, serde::Serialize)]
pub struct ClusterArgs {
    #[arg(
        required = true,
        num_args = 1..,
        value_name = "MATRIX|DIR",
        help = "Site matrices (zarr or h5), or a faba output directory (or a run record in it)",
        long_help = "Site matrices (zarr or h5): `{batch}_m6a_site` or `{batch}_atoi_site`,\n\
                     or a faba output directory (or a run record in it), whose site matrices\n\
                     of --modality are taken, one batch each."
    )]
    pub data_files: Vec<Box<str>>,

    #[arg(
        long = "modality",
        value_name = "m6a|atoi",
        help = "Which sites to take from an output directory (default: m6a, else atoi)"
    )]
    modality: Option<Box<str>>,

    #[arg(
        short = 'o',
        long = "output",
        value_name = "DIR",
        help = "Where to write the tables and the run record (default: the input's directory)"
    )]
    output: Option<Box<str>>,

    #[arg(
        long = "fit",
        value_name = "PATTERN",
        help = "Fit only the cells of matrices whose path matches PATTERN (`*`, `?`); assign the rest",
        long_help = "Fit only the cells of matrices whose path matches PATTERN (`*` any run,\n\
                     `?` one character), e.g. `--fit '*_wt_*'`. The other matrices' cells are\n\
                     assigned to the fitted clusters afterwards, as a control, from their reads\n\
                     as read: batch weights are not applied to them, so a control keeps its own\n\
                     level, and its `offset` measures that level against the fitted clusters.\n\
                     Default: fit all."
    )]
    fit: Option<Box<str>>,

    #[arg(
        long,
        default_value_t = false,
        help = "Do not weight reads to take the fitted batches out",
        long_help = "Do not weight reads to take the fitted batches out. By default, with more\n\
                     than one fitted batch, each read is weighted by the ratio of its site's\n\
                     pooled fraction to its batch's there (importance weights, capped and\n\
                     self-normalised per batch and site), so clusters cannot be made of\n\
                     batches. That also removes a site's difference between batches that\n\
                     differ in their cells."
    )]
    no_batch_correction: bool,

    #[arg(long, default_value_t = 50, help = "Drop cells covering fewer sites")]
    min_cell_sites: usize,

    #[arg(
        long,
        default_value_t = 20,
        help = "Drop sites covered by fewer fitted cells"
    )]
    min_site_cells: usize,

    #[arg(
        long,
        default_value_t = 5.0,
        help = "Prior strength: pseudo-reads pulling a cluster's site fraction toward the site's pooled fraction"
    )]
    prior_strength: f64,

    #[arg(
        long,
        default_value_t = 10,
        help = "Smallest cluster reported: smaller ones' cells go to their most probable larger cluster"
    )]
    min_cluster_cells: usize,

    #[arg(
        long,
        default_value_t = 4,
        help = "Independent chains, run in parallel"
    )]
    chains: usize,

    #[arg(
        long,
        default_value_t = 150,
        help = "Sweeps per chain, burn-in included"
    )]
    sweeps: usize,

    #[arg(
        long,
        default_value_t = 50,
        help = "Sweeps discarded at the start of each chain"
    )]
    burn_in: usize,

    #[arg(
        long,
        default_value_t = 2,
        help = "Keep every this-many-th sweep after the burn-in"
    )]
    thin: usize,

    #[arg(long, default_value_t = 10, help = "Split-merge proposals per sweep")]
    split_merge: usize,

    #[arg(
        long,
        default_value_t = 2,
        help = "Restricted Gibbs scans before each split-merge proposal"
    )]
    restricted_scans: usize,

    #[arg(
        long,
        default_value_t = 20,
        help = "Chains start from k-means into a number of clusters drawn from 1 up to this"
    )]
    start_clusters: usize,

    #[arg(
        long,
        default_value_t = 1,
        help = "Random seed (chain i uses seed + i)"
    )]
    seed: u64,
}

/// The matrices to read: each with its batch and whether it is fitted; and
/// the modality asked for.
fn inputs(args: &ClusterArgs) -> anyhow::Result<Vec<read::Input>> {
    let mut out = Vec::new();
    let fits = |path: &str| {
        args.fit
            .as_deref()
            .is_none_or(|p| crate::common::wildcard(p, path))
    };
    for f in &args.data_files {
        let folder: Box<str> = crate::run_record::default_output_dir(f)
            .to_string_lossy()
            .into();
        match crate::run_record::output_dir_of(f) {
            Some(dir) => {
                let sites = crate::site_analysis::output_dir::output_site_matrices(
                    &dir,
                    args.modality.as_deref(),
                )?;
                for (batch, path) in sites.matrices {
                    out.push(read::Input {
                        fit: fits(&path),
                        path,
                        folder: folder.clone(),
                        batch,
                    });
                }
            }
            None => {
                let name = crate::qc::layout::file_name(f);
                out.push(read::Input {
                    path: f.clone(),
                    folder,
                    batch: crate::qc::layout::matrix_batch(&name).unwrap_or(name),
                    fit: fits(f),
                });
            }
        }
    }
    anyhow::ensure!(
        out.iter().any(|i| i.fit),
        "no input matrix matches --fit {}",
        args.fit.as_deref().unwrap_or("")
    );
    qualify_batches(&mut out);
    Ok(out)
}

/// A batch name met in several folders is a different batch in each: those
/// are named `{folder}/{batch}`.
fn qualify_batches(inputs: &mut [read::Input]) {
    let mut folders: FxHashMap<Box<str>, FxHashSet<Box<str>>> = FxHashMap::default();
    for i in inputs.iter() {
        folders
            .entry(i.batch.clone())
            .or_default()
            .insert(i.folder.clone());
    }
    for i in inputs.iter_mut() {
        if folders[&i.batch].len() > 1 {
            i.batch = format!("{}/{}", i.folder, i.batch).into();
        }
    }
}

/// `-o`, else the first input directory, else the first matrix's directory.
fn output_dir(args: &ClusterArgs) -> PathBuf {
    if let Some(o) = &args.output {
        return PathBuf::from(o.as_ref());
    }
    crate::run_record::default_output_dir(&args.data_files[0])
}

/// Settings that would make the fit meaningless, refused up front.
fn check(args: &ClusterArgs) -> anyhow::Result<()> {
    anyhow::ensure!(
        args.prior_strength > 0.0 && args.prior_strength.is_finite(),
        "--prior-strength must be positive"
    );
    anyhow::ensure!(args.chains >= 1, "--chains must be at least 1");
    // Split R-hat needs at least two draws in each half of every chain.
    anyhow::ensure!(
        args.sweeps.saturating_sub(args.burn_in) / args.thin.max(1) >= 4,
        "--sweeps less --burn-in must leave at least 4 draws per chain (after --thin)"
    );
    Ok(())
}

/// Bytes the fit (not the reading) will hold at its peak, roughly: the
/// reads (raw, a copy of the fitted ones when some are left out, and
/// weighted), the picture chains start from, and per chain two f64 per site
/// for each live cluster (a guess of 32 clusters), split-merge scratch,
/// labels and kept draws.
fn memory_estimate(
    entries: usize,
    copied: bool,
    sites: usize,
    cells: usize,
    kept: usize,
    args: &ClusterArgs,
) -> usize {
    let reads = entries * (12 + 16 + if copied { 12 } else { 0 });
    let picture = sites * 32 * 4 + cells * 32 * 4;
    let chain = sites * 16 * 32 + sites * 33 + cells * 4 + kept * cells * 2;
    reads + picture + args.chains * chain
}

pub fn run_cluster(args: &ClusterArgs) -> anyhow::Result<()> {
    let dir = output_dir(args);
    let mut record = RunRecord::start("cluster", &dir.to_string_lossy())
        .inputs("matrices", &args.data_files)
        .options(args);
    let outcome = cluster(args, &dir, &mut record);
    if let Ok(name) = &outcome {
        record = record.file_name(name);
    }
    record.finish(&outcome);
    outcome.map(|_| ())
}

/// Fit, summarise and write; the run record's file name.
fn cluster(args: &ClusterArgs, dir: &Path, record: &mut RunRecord) -> anyhow::Result<String> {
    check(args)?;
    let inputs = inputs(args)?;
    let all = read::read_inputs(&inputs)?;
    let (data, kept) = read::filter(&all, args.min_cell_sites, args.min_site_cells);
    // Only the input cells are needed from here: let their reads go.
    let read_sites = all.sites.len();
    let read::SiteData {
        cells: all_cells, ..
    } = all;
    let fitted: Vec<usize> = (0..data.cells.len())
        .filter(|&c| data.cells[c].fit)
        .collect();
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
    let kept_draws = (args.sweeps - args.burn_in) / args.thin.max(1);
    let memory = memory_estimate(
        fit_reads.site.len(),
        matches!(fit_reads, Cow::Owned(_)),
        data.sites.len(),
        fitted.len(),
        kept_draws,
        args,
    );
    info!("about {} MB at the peak of the fit", memory >> 20);

    // Batches taken out of the fitted reads by importance weights.
    let batches: Vec<&str> = {
        let mut b: Vec<&str> = fitted
            .iter()
            .map(|&c| data.cells[c].batch.as_ref())
            .collect();
        b.sort_unstable();
        b.dedup();
        b
    };
    let (weighted, batch_stats) =
        if !args.no_batch_correction && batches.len() > 1 {
            let batch_of: Vec<u16> = fitted
                .iter()
                .map(|&c| {
                    batches
                        .binary_search(&data.cells[c].batch.as_ref())
                        .unwrap() as u16
                })
                .collect();
            let w = fit::batch_weights(&fit_reads, &batch_of, batches.len(), data.sites.len());
            info!(
            "importance weights take out {} batches: effective share of reads {:?}, {:.2}% capped",
            batches.len(),
            w.effective_share.iter().map(|e| format!("{e:.2}")).collect::<Vec<_>>(),
            100.0 * w.capped_share
        );
            let stats = serde_json::json!({
                "batches": batches,
                "effective_share": w.effective_share,
                "capped_share": w.capped_share,
            "max_weight": w.max_weight,
            });
            (
                read::Weighted::new(&fit_reads, Some((&w.on, &w.off))),
                Some(stats),
            )
        } else {
            (read::Weighted::new(&fit_reads, None), None)
        };

    let settings = crp::Settings {
        strength: args.prior_strength,
        split_merge: args.split_merge,
        restricted: args.restricted_scans,
        start_clusters: args.start_clusters,
        gibbs: true,
        fixed_alpha: None,
    };
    let model = crp::Model::new(&weighted, data.sites.len(), settings);
    let config = McmcConfig {
        n_samples: kept_draws.max(1),
        warmup: args.burn_in,
        thin: args.thin,
        seed: args.seed,
    };
    let chains = run_mcmc_parallel(&model, &config, args.chains);
    let posterior = summary::posterior(&chains, args.min_cluster_cells);
    if posterior.rhat_k > RHAT_OK || posterior.rhat_log_posterior > RHAT_OK {
        warn!(
            "chains disagree (split R-hat: K {:.2}, log posterior {:.2}): run more --sweeps or --chains",
            posterior.rhat_k, posterior.rhat_log_posterior
        );
    }
    let (part, tiles) = fit::settle(
        &fit_reads,
        &weighted,
        &summary::partition(&chains),
        args.min_cluster_cells,
        data.sites.len(),
        args.prior_strength,
    );
    let sizes = summary::sizes(&part.z, part.k);

    // Every kept cell's cluster, confidence and offset: a fitted cell's from
    // the fit (its offset from its reads as the fit saw them), any other
    // placed by its own reads.
    let mut fit_index = vec![None; data.cells.len()];
    for (i, &c) in fitted.iter().enumerate() {
        fit_index[c] = Some(i);
    }
    let placed: Vec<fit::Placement> = {
        use rayon::prelude::*;
        (0..data.cells.len())
            .into_par_iter()
            .map(|c| match fit_index[c] {
                Some(i) => {
                    let (s, m, u) = weighted.of(i);
                    let k = part.z[i];
                    fit::Placement {
                        cluster: k,
                        probability: part.confidence[i],
                        offset: tiles.offset(s, m, u, k as usize),
                    }
                }
                None => {
                    // A control cell, from its reads as read.
                    let (s, m, u) = data.reads.of(c);
                    tiles.place(s, m, u, &sizes)
                }
            })
            .collect()
    };
    let mut rows: Vec<write::CellRow> = all_cells
        .iter()
        .map(|cell| write::CellRow { cell, placed: None })
        .collect();
    for (p, &orig) in placed.iter().zip(&kept) {
        rows[orig].placed = Some(p);
    }

    std::fs::create_dir_all(dir)?;
    let m = &data.modality.to_ascii_lowercase();
    write::cells(&rows, &dir.join(format!("{m}_cell_clusters.parquet")))?;
    write::sites(
        &data.sites,
        &tiles,
        &dir.join(format!("{m}_site_clusters.parquet")),
    )?;

    record.set_result("posterior", &posterior);
    record.set_result("cluster_sizes", &sizes);
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
    record.set_result(
        "converged",
        &(posterior.rhat_k <= RHAT_OK && posterior.rhat_log_posterior <= RHAT_OK),
    );
    record.set_result("fit_memory_estimate_bytes", &memory);
    record.set_result("batch_weights", &batch_stats);
    info!(
        "{} clusters (sizes {:?}); K posterior {:?}; split R-hat K {:.2}, log posterior {:.2}",
        part.k, sizes, posterior.k, posterior.rhat_k, posterior.rhat_log_posterior
    );
    Ok(format!("{m}_cluster.run.json"))
}
