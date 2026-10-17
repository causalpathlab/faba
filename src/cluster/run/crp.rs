//! The mixture's fit: chains, the settled partition, every kept cell
//! placed, and its tables.

use super::super::{args, crp, fit, read, record_name, summary, write, FitData, Written, RHAT_OK};
use crate::run_record::RunRecord;
use legume_numeric::mcmc::engine::{run_mcmc_parallel, McmcConfig};
use log::{info, warn};
use std::path::Path;

/// Fit the mixture, write its tables, and fill the record; what was
/// written. `fit_reads` are the fitted cells' reads before any weights.
pub(in crate::cluster) fn crp_fit(
    args: &args::ClusterArgs,
    config: &McmcConfig,
    dir: &Path,
    record: &mut RunRecord,
    fit: FitData,
    fit_reads: &read::Reads,
) -> anyhow::Result<Written> {
    let FitData {
        data,
        fitted,
        weighted,
        held,
        all_cells,
        kept,
        ..
    } = fit;
    let settings = crp::Settings {
        strength: args::crp_strength(args),
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
        args::crp_strength(args),
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
    let m = data.modality.to_ascii_lowercase();
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
    Ok(Written {
        record: record_name(&m, args.model),
        modality: m,
    })
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
