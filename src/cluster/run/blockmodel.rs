//! The blockmodel's fit: chains, their tallies pooled where they agree
//! (else the best chain's alone), the cells left out placed against those
//! means, and its tables.

use super::super::{args, read, record_name, unified, write, FitData, Written};
use crate::run_record::RunRecord;
use legume_numeric::mcmc::engine::McmcConfig;
use log::{info, warn};
use rustc_hash::FxHashMap;
use std::path::Path;

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

/// Which chains the tables pool, those left out, and the quantities short
/// of convergence.
fn report_pooling(pooling: &unified::Pooling) {
    for e in &pooling.excluded {
        warn!("chain {} left out of the tables: {}", e.chain, e.reason);
    }
    match &pooling.fallback {
        Some(why) => warn!(
            "tables from the best chain ({}) alone: {why}",
            pooling.reference
        ),
        None => info!("tables from chains {:?} pooled", pooling.pooled),
    }
    if !pooling.not_converged.is_empty() {
        warn!(
            "not converged across the chains pooled: {}",
            pooling.not_converged.join(", ")
        );
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

/// Fit the blockmodel `args` set (the link or blocks preset, each switch
/// given overriding it), write its tables from the chains' pooled
/// posterior means (the best chain's when they do not pool; never a chain
/// the convergence verdict finds lagging) and each chain's trace, and fill
/// the record; what was written.
pub(in crate::cluster) fn blockmodel_fit(
    args: &args::ClusterArgs,
    config: &McmcConfig,
    dir: &Path,
    record: &mut RunRecord,
    fit: FitData,
) -> anyhow::Result<Written> {
    let data = fit.data;
    let (gene_of, names) = number_genes(&data.sites);
    // The links and the model take a while to build: not once cancelled.
    fit.progress.check()?;
    let links = unified::Links::new(fit.weighted, &gene_of);
    info!(
        "{} (cell, gene) links over {} genes",
        links.len(),
        names.len()
    );
    // The settings fitted, as the reports that draw them read them.
    let recorded = args::recorded(args);
    record.set_result("model", &recorded);
    let (a, b) = (recorded.settings.cell_roles, recorded.settings.gene_roles);
    fit.progress.check()?;
    let mut model = unified::Model::new(&links, &gene_of, data.sites.len(), recorded.settings)?;
    if let Some(h) = fit.held {
        model = model.with_heldout(h);
    }
    let chains = super::super::progress::run(&model, config, args.chains, fit.progress)?;
    let summary = unified::summarize(&model, &chains);
    // Whether it converged is the verdict on the traces and the R-hat, a
    // chain it finds lagging not pooled, and on what pooling found short of
    // convergence; the fitting view, or `progress::spark` after the fit,
    // tells it.
    fit.progress.summarised(summary.rhat_log_joint)?;
    let lagging = fit
        .progress
        .verdict()
        .expect("summarised")
        .lagging()
        .to_vec();
    let (summary, best, pooling) = unified::pooled(&model, &chains, summary, &lagging);
    report_pooling(&pooling);
    fit.progress.pooled(pooling.unsettled.clone());
    let verdict = fit.progress.verdict().expect("summarised");
    // Cells kept but not fitted, placed against the tables' means.
    let control_reads = data.reads.subset(fit.controls);
    let control_weighted = read::Weighted::new(&control_reads, None);
    let control_links = unified::Links::new(&control_weighted, &gene_of);
    let (placed, levels) = model.place(&control_links, &best, args.seed);
    let rows = fit.rows(a, (&best.cell, &best.level), (&placed, &levels));
    std::fs::create_dir_all(dir)?;
    let m = data.modality.to_ascii_lowercase();
    let path = |name: &str| dir.join(format!("{m}_{name}.parquet"));
    write::cell_roles(&rows, a, &path("cell_roles"))?;
    write_tables(&path, &fit, (&names, b), &best, &model)?;
    // Each chain's log joint as it ran, and the schedule it ran in.
    write::trace(&fit.progress.traces(), &path("trace"))?;
    record.set_result(
        "trace",
        &serde_json::json!({
            "sweeps": super::super::progress::sweeps(config),
            "burn_in": config.warmup,
        }),
    );
    record.set_result("blockmodel", &summary);
    record.set_result("pooling", &pooling);
    record.set_result("converged", &verdict.converged());
    info!(
        "blockmodel: split R-hat of the log joint {:.2}; chain agreement cells {:?}, genes {:?}; strong gene roles {:?}, top-gene overlap {:?}; held-out {:?} per read",
        summary.rhat_log_joint, summary.cell_agreement, summary.gene_agreement, summary.strong_gene_roles, summary.top_gene_overlap, summary.heldout_per_read
    );
    Ok(Written {
        record: record_name(&m, args.model),
        modality: m,
    })
}

/// The fit's tables beside the cells' roles, each at `path(name)`: the
/// genes' (`names`, in `b` gene roles) roles, the blocks, their genes,
/// cells and sites, the sites' reads per cell role, and the sites' own
/// statistics, from the tables' means `best`.
fn write_tables(
    path: &dyn Fn(&str) -> std::path::PathBuf,
    fit: &FitData,
    (names, b): (&[&str], usize),
    best: &unified::Means,
    model: &unified::Model,
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
        &best.site_lift(model.baseline()),
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
