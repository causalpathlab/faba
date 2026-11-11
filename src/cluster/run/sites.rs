//! The site blockmodel (`--model blocks`) fitted by variational inference
//! (see [`super::super::site_vi`]): restarts from k-means starts, the best
//! by ELBO kept, the cells kept but not fitted placed by a local step, and
//! the blockmodel's tables.

use super::super::links::Links;
use super::super::math::{logit, pooled_fractions};
use super::super::means::Means;
use super::super::settings::Weights;
use super::super::site_vi::{self, Data, Fit, Prior};
use super::super::{args, read, record_name, write, FitData, Written};
use super::vi::End;
use crate::cluster::count::vi::KMEANS_ITERS;
use crate::run_record::RunRecord;
use log::info;
use rustc_hash::FxHashMap;
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
    let links = Links::new(fit.weighted, &gene_of);
    info!(
        "{} (cell, gene) links over {} genes",
        links.len(),
        names.len()
    );
    // The role weights' global weights are fixed (uniform) in this fit.
    let mut recorded = args::recorded(args);
    recorded.settings.weights = Weights::Symmetric;
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
    let control_links = Links::new(&control_weighted, &gene_of);
    let control_data = Data::new(&control_links, &gene_of, names.len(), vi_data.mu.clone());
    let (placed, levels) = model.place(&control_data, &fitted);
    let (a, b) = (prior.cell_roles, prior.gene_roles);
    let rows = fit.rows(a, (&means.cell, &means.level), (&placed, &levels));
    super::vi::conclude(
        progress,
        record,
        "read",
        &ends,
        (best, kept),
        (warmup, options.epochs, options.tol),
        stopped,
        super::vi::empty_roles(&means.block_links(), b),
    )?;
    std::fs::create_dir_all(dir)?;
    let m = data.modality.to_ascii_lowercase();
    let path = |name: &str| dir.join(format!("{m}_{name}.parquet"));
    write::cell_roles(&rows, a, &path("cell_roles"))?;
    write_tables(&path, &fit, (&names, b), &means, &vi_data.mu)?;
    write::trace(&progress.traces(), &path("trace"))?;
    record.set_result("tables", &tables_entry(&m));
    Ok(Written {
        record: record_name(&m),
        modality: m,
    })
}

impl<'a> FitData<'a> {
    /// Every input cell's row: the fitted cells' role weights and levels
    /// from `fitted`, the placed ones' from `placed` (`k` weights per cell,
    /// in the order of the controls; levels empty when there are none), and
    /// none for cells not kept.
    fn rows<'b>(
        &self,
        k: usize,
        fitted: (&'b [f32], &[f64]),
        placed: (&'b [f32], &[f64]),
    ) -> Vec<write::RoleRow<'a>>
    where
        'b: 'a,
    {
        let mut rows: Vec<write::RoleRow> = self
            .all_cells
            .iter()
            .map(|cell| write::RoleRow {
                cell,
                membership: None,
                level: None,
            })
            .collect();
        for (cells, (weights, levels)) in [(self.fitted, fitted), (self.controls, placed)] {
            for (i, &c) in cells.iter().enumerate() {
                let row = &mut rows[self.kept[c]];
                row.membership = Some(&weights[i * k..(i + 1) * k]);
                row.level = levels.get(i).map(|&v| v as f32);
            }
        }
        rows
    }
}

/// Each site's gene, numbered in order of first appearance, and the genes'
/// names in that order.
pub(in crate::cluster) fn number_genes(sites: &[read::Site]) -> (Vec<u32>, Vec<&str>) {
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

/// The fit's tables beside the cells' roles, each at `path(name)`: the
/// genes' (`names`, in `b` gene roles) roles, the blocks, their genes,
/// cells and sites, the sites' reads per cell role, and the sites' own
/// statistics, from the fit's means `best` (`baseline`: each site's
/// pooled logit).
fn write_tables(
    path: &dyn Fn(&str) -> std::path::PathBuf,
    fit: &FitData,
    (names, b): (&[&str], usize),
    best: &Means,
    baseline: &[f64],
) -> anyhow::Result<()> {
    let data = fit.data;
    let k = best.blocks();
    write::gene_roles(names, &best.gene, b, &path("gene_roles"))?;
    write::blocks(
        &best.shift,
        &best.shift_sd,
        &best.block_links(),
        b,
        &path("blocks"),
    )?;
    write::block_members(
        &[("gene", names.to_vec())],
        &best.gene_block,
        &best.gene_enrichment(),
        k,
        b,
        &path("block_genes"),
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
        &path("block_cells"),
    )?;
    write::block_sites(
        &data.sites,
        &best.site_on,
        &best.site_off,
        &best.site_lift(baseline),
        k,
        b,
        &path("block_sites"),
    )?;
    write::site_cells(
        &data.sites,
        &best.site_on,
        &best.site_off,
        k,
        b,
        &path("site_cells"),
    )?;
    write::site_stats(
        &data.sites,
        &read::SiteStats::of(&data.reads, data.sites.len()),
        &path("site_stats"),
    )
}

/// The tables a blockmodel fit of modality `m` writes, as its run record
/// names them (under `tables` in its results): the cells' side (each cell's
/// roles, each block's cells), the genes' side (each gene's roles, the
/// blocks' shifts, genes and sites, the sites' pileups and counts) and the
/// restarts' trace, file names beside the record.
pub(crate) fn tables_entry(m: &str) -> serde_json::Value {
    let names = |names: &[&str]| -> Vec<String> {
        names.iter().map(|n| format!("{m}_{n}.parquet")).collect()
    };
    serde_json::json!({
        "modality": m,
        "cells": names(&["cell_roles", "block_cells"]),
        "genes": names(&["gene_roles", "blocks", "block_genes", "block_sites", "site_cells", "site_stats"]),
        "trace": format!("{m}_trace.parquet"),
    })
}
