//! The site blockmodel (`--model blocks`) fitted by variational inference
//! (see [`super::super::site_vi`]): restarts from k-means starts, the best
//! by ELBO kept, the cells kept but not fitted placed by a local step, and
//! the blockmodel's tables.

use super::super::crp::{logit, pooled_fractions};
use super::super::site_vi::{self, Data, Fit, Prior};
use super::super::{args, read, record_name, unified, write, FitData, Written};
use super::blockmodel::{number_genes, tables_entry, write_tables};
use super::vi::End;
use crate::cluster::count::vi::KMEANS_ITERS;
use crate::run_record::RunRecord;
use log::info;
use std::path::Path;

/// Fit the blockmodel of `args` (the blocks preset: infinite strength, the
/// shift, gene roles) by variational inference, write its tables and fill
/// the record; what was written.
pub(in crate::cluster) fn blocks_fit(
    args: &args::ClusterArgs,
    dir: &Path,
    record: &mut RunRecord,
    fit: FitData,
) -> anyhow::Result<Written> {
    let data = fit.data;
    let progress = fit.progress;
    let (gene_of, names) = number_genes(&data.sites);
    progress.check()?;
    let links = unified::Links::new(fit.weighted, &gene_of);
    info!(
        "{} (cell, gene) links over {} genes",
        links.len(),
        names.len()
    );
    // The role weights' global weights are fixed (uniform) in this fit.
    let mut recorded = args::recorded(args);
    recorded.settings.weights = unified::Weights::Symmetric;
    record.set_result("model", &recorded);
    let s = &recorded.settings;
    let prior = Prior {
        cell_roles: s.cell_roles,
        gene_roles: s.gene_roles,
        alpha: s.alpha,
        gene_alpha: s.gene_alpha,
        shift_sd: s.shift_sd.expect("the blocks preset has the shift"),
        cell_sd: s.cell_sd,
    };
    let options = super::vi::options(args);
    record.set_result("vi", &options);
    let mu: Vec<f64> = pooled_fractions(&links.site, &links.on, &links.off, data.sites.len())
        .into_iter()
        .map(logit)
        .collect();
    let vi_data = Data::new(&links, &gene_of, names.len(), mu);
    let model = Fit {
        data: &vi_data,
        prior,
    };
    progress.check()?;
    let warmup = KMEANS_ITERS;
    progress.begin_vi(warmup + options.epochs, warmup);
    let (best, ends) = super::super::progress::barred(progress, || {
        // The profiles serve the starts alone: moved into theirs, dropped
        // with it once they are made.
        let profiles = site_vi::Profiles::new(&model);
        let model = &model;
        super::vi::restarts(
            progress,
            (args.seed, args.chains),
            (warmup, options.epochs),
            move |rng, tick, stop| site_vi::kmeans_start(model, &profiles, rng, tick, stop),
            |start, rng, hooks| site_vi::fit(model, &options, start, rng, hooks),
            |f| End {
                elbo: f.elbo,
                per_unit: f.elbo_per_read(&vi_data),
                epochs: f.epochs,
                svi_epochs: f.svi_epochs,
                converged: f.converged,
            },
        )
    })?;
    if progress.cancelled() {
        return Err(super::super::progress::Cancelled.into());
    }
    let Some((best, fitted)) = best else {
        log::warn!("stopped before any restart finished an epoch: nothing to write");
        return Err(super::super::progress::Cancelled.into());
    };
    let kept = ends[best].expect("the best restart ended");
    // As the run stood when the fits ended: a stop during the tables'
    // passes after does not make the fit one stopped early.
    let stopped = progress.stopped_early().is_some();
    let means = model.means(&fitted);
    // Cells kept but not fitted, each a local step against the globals.
    let control_reads = data.reads.subset(fit.controls);
    let control_weighted = read::Weighted::new(&control_reads, None);
    let control_links = unified::Links::new(&control_weighted, &gene_of);
    let control_data = Data::new(&control_links, &gene_of, names.len(), vi_data.mu.clone());
    let (placed, levels) = model.place(&control_data, &fitted);
    let (a, b) = (prior.cell_roles, prior.gene_roles);
    let rows = fit.rows(a, (&means.cell, &means.level), (&placed, &levels));
    super::vi::conclude(
        progress,
        record,
        ("blocks", "read"),
        &ends,
        (best, kept),
        (warmup, options.epochs, options.tol),
        stopped,
        serde_json::json!({"empty_roles": super::vi::empty_roles(&means.block_links(), b)}),
    )?;
    std::fs::create_dir_all(dir)?;
    let m = data.modality.to_ascii_lowercase();
    let path = |name: &str| dir.join(format!("{m}_{name}.parquet"));
    write::cell_roles(&rows, a, &path("cell_roles"))?;
    write_tables(&path, &fit, (&names, b), &means, &vi_data.mu)?;
    write::trace(&progress.traces(), &path("trace"))?;
    record.set_result("tables", &tables_entry(&m));
    Ok(Written {
        record: record_name(&m, args.model),
        modality: m,
    })
}
