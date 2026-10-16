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
                     assigned to the fitted clusters afterwards, as a control. Default: fit all."
    )]
    fit: Option<Box<str>>,

    #[arg(
        long,
        default_value_t = false,
        help = "Do not cancel each fitted batch's shift at each site",
        long_help = "Do not cancel each fitted batch's shift at each site. By default, with\n\
                     more than one fitted batch, each batch's logit offset from the site's\n\
                     pooled fraction is cancelled by reweighting its reads before clustering,\n\
                     so clusters cannot be made of batches. That also removes a site's\n\
                     difference between batches that differ in their cells."
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
        default_value_t = 15,
        help = "Nearest cells offered as a split-merge anchor's partner (0: random anchors only)"
    )]
    anchor_neighbours: usize,

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
            .is_none_or(|p| crate::site_analysis::pileup::wildcard(p, path))
    };
    for f in &args.data_files {
        match crate::run_record::output_dir_of(f) {
            Some(dir) => {
                let sites =
                    crate::site_analysis::output_dir::output_sites(&dir, args.modality.as_deref())?;
                anyhow::ensure!(
                    !sites.matrices.is_empty(),
                    "{} has no {}_site matrices",
                    dir.display(),
                    sites.modality
                );
                for (batch, path) in sites.matrices {
                    out.push(read::Input {
                        fit: fits(&path),
                        path,
                        batch,
                    });
                }
            }
            None => {
                let name = crate::qc::layout::file_name(f);
                out.push(read::Input {
                    path: f.clone(),
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
    Ok(out)
}

/// `-o`, else the first input directory, else the first matrix's directory.
fn output_dir(args: &ClusterArgs) -> PathBuf {
    if let Some(o) = &args.output {
        return PathBuf::from(o.as_ref());
    }
    let first = &args.data_files[0];
    crate::run_record::output_dir_of(first).unwrap_or_else(|| {
        Path::new(first.as_ref())
            .parent()
            .map(Path::to_path_buf)
            .unwrap_or_default()
    })
}

/// Bytes the fit (not the reading) will hold at its peak, roughly: the reads, a pair of counts
/// per site for each live cluster in each chain (a guess of 32 clusters),
/// split-merge scratch, and the kept draws.
fn memory_estimate(entries: usize, sites: usize, cells: usize, args: &ClusterArgs) -> usize {
    let kept = args.sweeps.saturating_sub(args.burn_in) / args.thin.max(1) + 1;
    let per_chain = sites * 8 * 32 + sites * 17 + cells * 4 + kept * cells * 2;
    entries * 12 + args.chains * per_chain
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
    let inputs = inputs(args)?;
    let all = read::read_inputs(&inputs)?;
    let (data, kept) = read::filter(&all, args.min_cell_sites, args.min_site_cells);
    // Only the input cells are needed from here: let their reads go.
    let read_sites = all.sites.len();
    let all_cells = all.cells;
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
    let fit_reads = read::Reads::from_cells(
        fitted
            .iter()
            .map(|&c| {
                let (s, m, u) = data.reads.of(c);
                (0..s.len()).map(|i| (s[i], m[i], u[i])).collect()
            })
            .collect(),
    );
    let memory = memory_estimate(fit_reads.site.len(), data.sites.len(), fitted.len(), args);
    info!("about {} MB at the peak of the fit", memory >> 20);

    let settings = crp::Settings {
        strength: args.prior_strength,
        split_merge: args.split_merge,
        restricted: args.restricted_scans,
        start_clusters: args.start_clusters,
        neighbours: args.anchor_neighbours,
        gibbs: true,
        fixed_alpha: None,
    };
    let batches: Vec<&str> = {
        let mut b: Vec<&str> = fitted
            .iter()
            .map(|&c| data.cells[c].batch.as_ref())
            .collect();
        b.sort_unstable();
        b.dedup();
        b
    };
    let shift = (!args.no_batch_correction && batches.len() > 1).then(|| {
        let batch_of: Vec<u16> = fitted
            .iter()
            .map(|&c| {
                batches
                    .binary_search(&data.cells[c].batch.as_ref())
                    .unwrap() as u16
            })
            .collect();
        info!(
            "cancelling the shifts of {} batches at each site",
            batches.len()
        );
        fit::batch_shifts(&fit_reads, &batch_of, batches.len(), data.sites.len())
    });
    let model = crp::Model::new(&fit_reads, data.sites.len(), settings, shift.as_deref());
    let kept_draws = args.sweeps.saturating_sub(args.burn_in) / args.thin.max(1);
    let config = McmcConfig {
        n_samples: kept_draws.max(1),
        warmup: args.burn_in,
        thin: args.thin,
        seed: args.seed,
    };
    let chains = run_mcmc_parallel(&model, &config, args.chains.max(1));
    let posterior = summary::posterior(&chains, args.min_cluster_cells);
    let draws_part = summary::partition(&chains);
    if posterior.rhat_k > RHAT_OK || posterior.rhat_log_posterior > RHAT_OK {
        warn!(
            "chains disagree (split R-hat: K {:.2}, log posterior {:.2}): run more --sweeps or --chains",
            posterior.rhat_k, posterior.rhat_log_posterior
        );
    }

    // Clusters' reads from the fitted cells; every kept cell's cluster and offset.
    let (part, tiles) = fit::settle(
        &fit_reads,
        &draws_part,
        args.min_cluster_cells,
        data.sites.len(),
        args.prior_strength,
    );
    let mut sizes = vec![0usize; part.k];
    part.z.iter().for_each(|&k| sizes[k as usize] += 1);
    let mut cluster_of: Vec<Option<(u16, f64)>> = vec![None; data.cells.len()];
    for (i, &c) in fitted.iter().enumerate() {
        cluster_of[c] = Some((part.z[i], part.confidence[i]));
    }
    for c in (0..data.cells.len()).filter(|&c| !data.cells[c].fit) {
        let (s, m, u) = data.reads.of(c);
        let p = tiles.assign(s, m, u, &sizes);
        let (k, conf) = p
            .iter()
            .enumerate()
            .max_by(|a, b| a.1.total_cmp(b.1))
            .map(|(k, &v)| (k as u16, v))
            .unwrap();
        cluster_of[c] = Some((k, conf));
    }
    let mut rows: Vec<write::CellRow> = all_cells
        .iter()
        .map(|cell| write::CellRow {
            cell,
            kept: false,
            cluster: None,
            confidence: None,
            offset: None,
        })
        .collect();
    for (c, &orig) in kept.cells.iter().enumerate() {
        let (k, conf) = cluster_of[c].unwrap();
        let (s, m, u) = data.reads.of(c);
        rows[orig].kept = true;
        rows[orig].cluster = Some(k);
        rows[orig].confidence = Some(conf);
        rows[orig].offset = Some(tiles.offset(s, m, u, k as usize));
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
    record.set_result(
        "batches_corrected",
        &shift.as_ref().map_or(0, |_| batches.len()),
    );
    info!(
        "{} clusters (sizes {:?}); K posterior {:?}; split R-hat K {:.2}, log posterior {:.2}",
        part.k, sizes, posterior.k, posterior.rhat_k, posterior.rhat_log_posterior
    );
    Ok(format!("{m}_cluster.run.json"))
}
