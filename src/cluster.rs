//! `faba cluster`: cells clustered by their own site reads.
//!
//! A Chinese-restaurant-process mixture over cells (see [`crp`]), run as
//! several chains and summarised without a cells x cells matrix (see
//! [`summary`]). Cells from inputs outside `--fit` are assigned to the fitted
//! clusters afterwards, as a control.

use crate::run_record::RunRecord;
use legume_numeric::mcmc::engine::{run_mcmc_parallel, McmcConfig};
use log::{info, warn};
use rustc_hash::FxHashMap;
use std::borrow::Cow;
use std::path::Path;

mod args;
mod blocks;
mod browse;
mod command;
mod crp;
mod embed;
mod fit;
mod link;
mod read;
mod summary;
mod write;

use args::{check, inputs, output_dir, with_inputs};
pub use args::{ClusterArgs, ModelKind};
pub use browse::{run_modules, ModulesArgs};
use command::command;

#[cfg(test)]
mod tests;

/// Split R-hat above this is reported as not converged.
const RHAT_OK: f64 = 1.05;

/// Bytes the fit (not the reading) will hold at its peak, roughly: the
/// reads (raw, a copy of the fitted ones when some are left out, and
/// weighted), then the model's own per-chain state.
fn memory_estimate(
    entries: usize,
    copied: bool,
    sites: usize,
    cells: usize,
    kept: usize,
    args: &ClusterArgs,
) -> usize {
    let reads = entries * (12 + 16 + if copied { 12 } else { 0 });
    // Links hold their own copy of the weighted reads.
    let links = entries * 20;
    let (shared, chain) = match args.model {
        // The picture chains start from; per chain two f64 per site for
        // each live cluster (a guess of 32), split-merge scratch, labels
        // and kept draws.
        ModelKind::Crp => (
            sites * 32 * 4 + cells * 32 * 4,
            sites * 16 * 32 + sites * 33 + cells * 4 + kept * cells * 2,
        ),
        // Per chain: two f64 per site and community, a u32 and an f64 sum
        // per cell and community, and each entry's link label twice (the
        // chain's and its tally's).
        ModelKind::Link => {
            let k = args.communities;
            (links, sites * k * 16 + cells * k * 12 + entries * 4)
        }
        // Per chain: role counts (u32) and their sums (f64) for cells and
        // genes (at most one gene per site), two labels per link, and the
        // tally's links per cell and gene and reads per site in each block.
        ModelKind::Blocks => {
            let (a, b) = (args.cell_roles, args.gene_roles);
            (
                links,
                (cells * a + sites * b) * 12
                    + entries * 4
                    + (cells + sites) * a * b * 8
                    + sites * a * b * 16,
            )
        }
    };
    reads + shared + args.chains * chain
}

/// Cluster as `args` say; the modality fitted, none when no inputs were
/// chosen.
pub fn run_cluster(args: &ClusterArgs) -> anyhow::Result<Option<String>> {
    let Some(args) = &with_inputs(args, "cluster", false)? else {
        return Ok(None);
    };
    let dir = output_dir(args);
    let mut record = RunRecord::start("cluster", &dir.to_string_lossy())
        .inputs("matrices", &args.data_files)
        .options(args);
    let outcome = cluster(args, &dir, &mut record);
    if let Ok(name) = &outcome {
        record = record.file_name(name);
        // The command, to read, edit and run again, beside its record.
        let name = name.replace(".run.json", ".cmd.sh");
        let script = crate::script::Script {
            name: &name,
            made_by: "faba cluster",
            guard: None,
            cd: false,
            head: 1,
            replace: true,
        };
        // The fit is done and recorded whether or not the script is written.
        if let Err(e) = crate::script::write(&dir, &script, &command(args)) {
            warn!("the command was not saved: {e:#}");
        }
    }
    record.finish(&outcome);
    outcome.map(|name| name.strip_suffix("_cluster.run.json").map(String::from))
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
    // Entries held out of the fit, the same for either model.
    let (fit_reads, held): (Cow<read::Reads>, Option<read::Reads>) = if args.holdout > 0.0 {
        let (fit, held) = read::hold_out(&fit_reads, args.holdout, args.seed);
        (Cow::Owned(fit), Some(held))
    } else {
        (fit_reads, None)
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

    let (weighted, batch_stats) = weigh_batches(args, &data, &fitted, &fit_reads);

    let config = McmcConfig {
        n_samples: kept_draws.max(1),
        warmup: args.burn_in,
        thin: args.thin,
        seed: args.seed,
    };
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
    let fit = LinkFit {
        data: &data,
        fitted: &fitted,
        weighted: &weighted,
        held: held.as_ref(),
        all_cells: &all_cells,
        kept: &kept,
    };
    match args.model {
        ModelKind::Crp => crp_fit(args, &config, dir, record, fit, &fit_reads),
        ModelKind::Link => link_fit(args, &config, dir, record, fit),
        ModelKind::Blocks => blocks_fit(args, &config, dir, record, fit),
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

/// Fit the mixture, write its tables, and fill the record; the run record's
/// file name. `fit_reads` are the fitted cells' reads before any weights.
fn crp_fit(
    args: &ClusterArgs,
    config: &McmcConfig,
    dir: &Path,
    record: &mut RunRecord,
    fit: LinkFit,
    fit_reads: &read::Reads,
) -> anyhow::Result<String> {
    let LinkFit {
        data,
        fitted,
        weighted,
        held,
        all_cells,
        kept,
    } = fit;
    let settings = crp::Settings {
        strength: args.prior_strength,
        split_merge: args.split_merge,
        restricted: args.restricted_scans,
        start_clusters: args.start_clusters,
        gibbs: true,
        fixed_alpha: None,
    };
    let mut model = crp::Model::new(weighted, data.sites.len(), settings);
    if let Some(h) = held {
        model = model.with_heldout(h);
    }
    let chains = run_mcmc_parallel(&model, config, args.chains);
    let posterior = summary::posterior(&chains, args.min_cluster_cells);
    if posterior.rhat_k > RHAT_OK || posterior.rhat_log_posterior > RHAT_OK {
        warn!(
            "chains disagree (split R-hat: K {:.2}, log posterior {:.2}): run more --sweeps or --chains",
            posterior.rhat_k, posterior.rhat_log_posterior
        );
    }
    let (part, tiles) = fit::settle(
        fit_reads,
        weighted,
        &summary::partition(&chains),
        args.min_cluster_cells,
        data.sites.len(),
        args.prior_strength,
    );
    let sizes = summary::sizes(&part.z, part.k);
    let placed = place_cells(data, fitted, weighted, &part, &tiles, &sizes);
    let mut rows: Vec<write::CellRow> = all_cells
        .iter()
        .map(|cell| write::CellRow { cell, placed: None })
        .collect();
    for (p, &orig) in placed.iter().zip(kept) {
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
        "converged",
        &(posterior.rhat_k <= RHAT_OK && posterior.rhat_log_posterior <= RHAT_OK),
    );
    info!(
        "{} clusters (sizes {:?}); K posterior {:?}; split R-hat K {:.2}, log posterior {:.2}",
        part.k, sizes, posterior.k, posterior.rhat_k, posterior.rhat_log_posterior
    );
    Ok(format!("{m}_cluster.run.json"))
}

/// Every kept cell's cluster, confidence and offset: a fitted cell's from
/// the fit (its offset from its reads as the fit saw them), any other
/// placed by its own reads.
fn place_cells(
    data: &read::SiteData,
    fitted: &[usize],
    weighted: &read::Weighted,
    part: &summary::Partition,
    tiles: &fit::Tiles,
    sizes: &[usize],
) -> Vec<fit::Placement> {
    let mut fit_index = vec![None; data.cells.len()];
    for (i, &c) in fitted.iter().enumerate() {
        fit_index[c] = Some(i);
    }
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
                tiles.place(s, m, u, sizes)
            }
        })
        .collect()
}

/// What each model needs from the read and filtered input.
struct LinkFit<'a> {
    data: &'a read::SiteData,
    fitted: &'a [usize],
    weighted: &'a read::Weighted<'a>,
    held: Option<&'a read::Reads>,
    all_cells: &'a [read::Cell],
    kept: &'a [usize],
}

impl<'a> LinkFit<'a> {
    /// Every input cell's row: the fitted cells' memberships (and levels)
    /// from `fitted`, the placed ones' from `placed` (`k` per cell, in the
    /// order of [`controls`]), and none for cells not kept.
    fn rows<'b>(
        &self,
        k: usize,
        fitted: (&'b [f32], Option<&[f64]>),
        placed: (&'b [f32], Option<&[f64]>),
    ) -> Vec<write::CommunityRow<'a>>
    where
        'b: 'a,
    {
        let mut rows: Vec<write::CommunityRow> = self
            .all_cells
            .iter()
            .map(|cell| write::CommunityRow {
                cell,
                membership: None,
                level: None,
            })
            .collect();
        let controls: Vec<usize> = (0..self.data.cells.len())
            .filter(|&c| !self.data.cells[c].fit)
            .collect();
        for (cells, (weights, levels)) in [(self.fitted, fitted), (&controls[..], placed)] {
            for (i, &c) in cells.iter().enumerate() {
                let row = &mut rows[self.kept[c]];
                row.membership = Some(&weights[i * k..(i + 1) * k]);
                row.level = levels.map(|v| v[i] as f32);
            }
        }
        rows
    }
}

/// A warning when the chains disagree.
fn warn_unconverged(rhat: f64) {
    if rhat > RHAT_OK {
        warn!("chains disagree (split R-hat of the log joint {rhat:.2}): run more --sweeps or --chains");
    }
}

/// Fit the link-community model, write its tables, and fill the record;
/// the run record's file name.
fn link_fit(
    args: &ClusterArgs,
    config: &McmcConfig,
    dir: &Path,
    record: &mut RunRecord,
    fit: LinkFit,
) -> anyhow::Result<String> {
    let data = fit.data;
    let (gene_of, names) = number_genes(&data.sites);
    let links = link::Links::new(fit.weighted, &gene_of);
    info!(
        "{} (cell, gene) links over {} genes",
        links.len(),
        names.len()
    );
    let settings = link::Settings {
        k: args.communities,
        alpha: args.membership_alpha,
        gamma: args.community_gamma,
        strength: args.prior_strength,
    };
    let mut model = link::Model::new(&links, data.sites.len(), settings);
    if let Some(h) = fit.held {
        model = model.with_heldout(h);
    }
    let chains = run_mcmc_parallel(&model, config, args.chains);
    let summary = link::summarize(&model, &links, &chains);
    warn_unconverged(summary.rhat_log_joint);
    let best = chains[summary.best_chain]
        .last()
        .expect("a chain has draws")
        .means();
    let k = args.communities;
    // Cells kept but not fitted, placed against the best chain's communities.
    let profile = model.profile(&best.z);
    let control_reads = controls(data);
    let control_weighted = read::Weighted::new(&control_reads, None);
    let control_links = link::Links::new(&control_weighted, &gene_of);
    let placed = model.place(&control_links, &profile, &summary.global_weights, args.seed);
    let rows = fit.rows(k, (&best.membership, None), (&placed, None));
    std::fs::create_dir_all(dir)?;
    let m = &data.modality.to_ascii_lowercase();
    write::communities(
        &rows,
        k,
        "community",
        &dir.join(format!("{m}_cell_communities.parquet")),
    )?;
    write::gene_communities(
        &links.genes_of(&gene_of),
        &names,
        &best.z,
        k,
        &dir.join(format!("{m}_gene_communities.parquet")),
    )?;
    write::community_sites(
        &data.sites,
        &profile,
        &dir.join(format!("{m}_community_sites.parquet")),
    )?;
    record.set_result("link", &summary);
    record.set_result("converged", &(summary.rhat_log_joint <= RHAT_OK));
    info!(
        "link communities: read shares {:?}; split R-hat of the log joint {:.2}; held-out {:?} per read",
        summary.read_shares, summary.rhat_log_joint, summary.heldout_per_read
    );
    Ok(format!("{m}_cluster.run.json"))
}
/// The reads, as read, of the kept cells left out of the fit, in order.
fn controls(data: &read::SiteData) -> read::Reads {
    let cells: Vec<usize> = (0..data.cells.len())
        .filter(|&c| !data.cells[c].fit)
        .collect();
    data.reads.subset(&cells)
}

/// Each site's gene, numbered in order of first appearance, and the genes'
/// names in that order.
fn number_genes(sites: &[read::Site]) -> (Vec<u32>, Vec<&str>) {
    let mut index: FxHashMap<&str, u32> = FxHashMap::default();
    let mut names = Vec::new();
    let gene_of = sites
        .iter()
        .map(|s| {
            *index.entry(s.gene.as_ref()).or_insert_with(|| {
                names.push(s.gene.as_ref());
                names.len() as u32 - 1
            })
        })
        .collect();
    (gene_of, names)
}

/// Fit the two-role blockmodel, write its tables, and fill the record; the
/// run record's file name.
fn blocks_fit(
    args: &ClusterArgs,
    config: &McmcConfig,
    dir: &Path,
    record: &mut RunRecord,
    fit: LinkFit,
) -> anyhow::Result<String> {
    let data = fit.data;
    let (gene_of, names) = number_genes(&data.sites);
    let links = link::Links::new(fit.weighted, &gene_of);
    info!(
        "{} (cell, gene) links over {} genes",
        links.len(),
        names.len()
    );
    let settings = blocks::Settings {
        cell_roles: args.cell_roles,
        gene_roles: args.gene_roles,
        alpha: args.membership_alpha,
        gene_alpha: args.gene_alpha,
        shift_sd: args.shift_sd,
        cell_sd: args.cell_sd,
    };
    let mut model = blocks::Model::new(&links, &gene_of, data.sites.len(), settings);
    if let Some(h) = fit.held {
        model = model.with_heldout(h, &gene_of);
    }
    let chains = run_mcmc_parallel(&model, config, args.chains);
    let summary = blocks::summarize(&model, &chains);
    warn_unconverged(summary.rhat_log_joint);
    let best = chains[summary.best_chain]
        .last()
        .expect("a chain has draws")
        .means(true);
    let a = args.cell_roles;
    // Cells kept but not fitted, placed against the best chain.
    let control_reads = controls(data);
    let control_weighted = read::Weighted::new(&control_reads, None);
    let control_links = link::Links::new(&control_weighted, &gene_of);
    let (placed, levels) = model.place(&control_links, &gene_of, &best, args.seed);
    let rows = fit.rows(
        a,
        (&best.cell, Some(&best.cell_base)),
        (&placed, Some(&levels)),
    );
    std::fs::create_dir_all(dir)?;
    let m = &data.modality.to_ascii_lowercase();
    write::communities(
        &rows,
        a,
        "role",
        &dir.join(format!("{m}_cell_roles.parquet")),
    )?;
    write::gene_roles(
        &names,
        &best.gene,
        args.gene_roles,
        &dir.join(format!("{m}_gene_roles.parquet")),
    )?;
    let (k, b) = (best.blocks(), args.gene_roles);
    let block_links: Vec<f64> = (0..k)
        .map(|j| best.gene_block.iter().skip(j).step_by(k).sum())
        .collect();
    write::blocks(
        &summary.shifts,
        &summary.shift_sd,
        &block_links,
        b,
        &dir.join(format!("{m}_blocks.parquet")),
    )?;
    write::block_members(
        &[("gene", names.clone())],
        &best.gene_block,
        &best.gene_enrichment(),
        k,
        b,
        &dir.join(format!("{m}_block_genes.parquet")),
    )?;
    let fitted_cells: Vec<&read::Cell> = fit.fitted.iter().map(|&c| &data.cells[c]).collect();
    write::block_members(
        &[
            (
                "barcode",
                fitted_cells.iter().map(|c| c.barcode.as_ref()).collect(),
            ),
            (
                "batch",
                fitted_cells.iter().map(|c| c.batch.as_ref()).collect(),
            ),
        ],
        &best.cell_block,
        &best.cell_enrichment(),
        k,
        b,
        &dir.join(format!("{m}_block_cells.parquet")),
    )?;
    write::block_sites(
        &data.sites,
        &best.site_on,
        &best.site_off,
        &best.site_lift(model.baseline()),
        k,
        b,
        &dir.join(format!("{m}_block_sites.parquet")),
    )?;
    write::site_stats(
        &data.sites,
        &read::SiteStats::of(&data.reads, data.sites.len()),
        &dir.join(format!("{m}_site_stats.parquet")),
    )?;
    record.set_result("blocks", &summary);
    record.set_result("converged", &(summary.rhat_log_joint <= RHAT_OK));
    info!(
        "blocks: split R-hat of the log joint {:.2}; chain agreement cells {:?}, genes {:?}; strong gene roles {:?}, top-gene overlap {:?}; held-out {:?} per read",
        summary.rhat_log_joint, summary.cell_agreement, summary.gene_agreement, summary.strong_gene_roles, summary.top_gene_overlap, summary.heldout_per_read
    );
    Ok(format!("{m}_cluster.run.json"))
}
