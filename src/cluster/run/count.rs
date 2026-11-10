//! `faba cluster --modality count`: the count blockmodel's fit (see
//! [`crate::cluster::count`]) of gene count matrices, from reading them to
//! its tables.

use super::super::command::COUNT_MODALITY;
use super::super::count::{self, vi, Counts};
use super::super::read::{walk_matrix, Visit};
use super::super::summary::main_labels;
use super::super::{args, progress, read, write, Written};
use super::vi::End;
use crate::run_record::RunRecord;
use data_beans::aux::feature_rows::split_count_row;
use log::info;
use rustc_hash::FxHashMap;
use std::path::Path;

/// The count fit's run record, which `faba cluster --within` and the
/// viewer look for by this name.
pub(in crate::cluster) const RECORD: &str = "count_blockmodel.run.json";

/// The genes' names, the cells, and each cell's `(gene, count)`s.
type CountsRead = (Vec<Box<str>>, Vec<read::Cell>, Vec<Vec<(u32, u32)>>);

/// A count matrix's reader: its spliced rows matched to the genes of the
/// inputs before it, its cells and their counts appended.
struct CountVisit<'v> {
    input: &'v read::Input,
    rows: Vec<Option<u32>>,
    first: usize,
    read: &'v mut CountsRead,
    index: &'v mut FxHashMap<Box<str>, u32>,
}

impl Visit for CountVisit<'_> {
    fn rows(&mut self, names: &[Box<str>]) -> anyhow::Result<()> {
        let (genes, ..) = &mut *self.read;
        for name in names {
            let gene = match split_count_row(name) {
                Some((gene, false)) => gene,
                _ => {
                    self.rows.push(None);
                    continue;
                }
            };
            let next = genes.len() as u32;
            self.rows
                .push(Some(*self.index.entry(gene.into()).or_insert_with(|| {
                    genes.push(gene.into());
                    next
                })));
        }
        anyhow::ensure!(
            self.rows.iter().any(Option::is_some),
            "{} holds no spliced count rows",
            self.input.path
        );
        Ok(())
    }

    fn columns(&mut self, barcodes: &[Box<str>]) {
        let (_, cells, per_cell) = &mut *self.read;
        cells.extend(barcodes.iter().map(|b| read::Cell {
            barcode: b.clone(),
            batch: self.input.batch.clone(),
            fit: self.input.fit,
        }));
        per_cell.resize_with(cells.len(), Vec::new);
    }

    fn entry(&mut self, row: usize, col: usize, n: u32) {
        if let Some(Some(g)) = self.rows.get(row) {
            self.read.2[self.first + col].push((*g, n));
        }
    }
}

/// Every cell's spliced counts over the genes of every input.
fn read_counts(
    inputs: &[read::Input],
    progress: &progress::Progress,
) -> anyhow::Result<CountsRead> {
    let mut read: CountsRead = Default::default();
    let mut index = FxHashMap::default();
    for input in inputs {
        let mut visit = CountVisit {
            input,
            rows: Vec::new(),
            first: read.1.len(),
            read: &mut read,
            index: &mut index,
        };
        walk_matrix(&input.path, progress, &mut visit)?;
    }
    Ok(read)
}

/// The `n` genes of highest dispersion (variance over mean of the counts,
/// each cell's scaled to the median total) among those with any count,
/// in their order; all of them when there are fewer.
fn variable_genes(per_cell: &[&[(u32, u32)]], genes: usize, n: usize) -> Vec<usize> {
    let totals: Vec<f64> = per_cell
        .iter()
        .map(|c| c.iter().map(|e| e.1 as f64).sum())
        .collect();
    let mut sorted = totals.clone();
    sorted.sort_by(f64::total_cmp);
    let median = sorted
        .get(sorted.len() / 2)
        .copied()
        .unwrap_or(1.0)
        .max(1.0);
    let (mut sum, mut sq) = (vec![0.0f64; genes], vec![0.0f64; genes]);
    for (cell, &t) in per_cell.iter().zip(&totals) {
        let scale = median / t.max(1.0);
        for &(g, x) in cell.iter() {
            let y = x as f64 * scale;
            sum[g as usize] += y;
            sq[g as usize] += y * y;
        }
    }
    let cells = per_cell.len().max(1) as f64;
    let mut ranked: Vec<(usize, f64)> = (0..genes)
        .filter(|&g| sum[g] > 0.0)
        .map(|g| {
            let mean = sum[g] / cells;
            (g, (sq[g] / cells - mean * mean) / mean)
        })
        .collect();
    ranked.sort_by(|x, y| y.1.total_cmp(&x.1));
    let mut chosen: Vec<usize> = ranked.into_iter().take(n).map(|(g, _)| g).collect();
    chosen.sort_unstable();
    chosen
}

/// Fit `args`' count matrices by variational inference (see
/// [`count::vi`]), write the tables in `dir` from the restart of highest
/// ELBO, and fill the record.
pub(in crate::cluster) fn count_fit(
    args: &args::ClusterArgs,
    dir: &Path,
    record: &mut RunRecord,
    progress: &progress::Progress,
) -> anyhow::Result<Written> {
    let inputs = args::inputs(args)?;
    let (all_genes, cells, per_cell) = read_counts(&inputs, progress)?;
    // Genes chosen over the fitted cells alone (`--fit`).
    let fitted_rows: Vec<&[(u32, u32)]> = (per_cell.iter().zip(&cells))
        .filter(|(_, c)| c.fit)
        .map(|(r, _)| r.as_slice())
        .collect();
    let chosen = variable_genes(&fitted_rows, all_genes.len(), args.count_genes);
    drop(fitted_rows);
    let mut batches: Vec<&str> = (cells.iter().filter(|c| c.fit))
        .map(|c| c.batch.as_ref())
        .collect();
    batches.sort_unstable();
    batches.dedup();
    if batches.len() > 1 {
        log::warn!(
            "count: {} batches fitted together, with no batch correction: a batch effect \
             can come out as cell roles",
            batches.len()
        );
    }
    let mut renumber = vec![u32::MAX; all_genes.len()];
    for (i, &g) in chosen.iter().enumerate() {
        renumber[g] = i as u32;
    }
    let names: Vec<&str> = chosen.iter().map(|&g| all_genes[g].as_ref()).collect();
    // Each cell over the chosen genes; a cell with too few is not fitted.
    let restricted: Vec<Vec<(u32, u32)>> = per_cell
        .iter()
        .map(|c| {
            c.iter()
                .filter(|e| renumber[e.0 as usize] != u32::MAX)
                .map(|&(g, x)| (renumber[g as usize], x))
                .collect()
        })
        .collect();
    // The fitted cells (`--fit`) with enough of the chosen genes.
    let kept: Vec<usize> = (0..cells.len())
        .filter(|&c| cells[c].fit && restricted[c].len() >= args.min_cell_sites)
        .collect();
    anyhow::ensure!(
        kept.len() >= 2 && !names.is_empty(),
        "{} cells and {} genes pass the filters: too few to cluster",
        kept.len(),
        names.len()
    );
    drop(per_cell);
    let counts = Counts::from_cells(names.len(), kept.iter().map(|&c| restricted[c].as_slice()));
    info!(
        "count: {} of {} cells kept, {} of {} genes (most variable), {} nonzero counts",
        kept.len(),
        cells.len(),
        names.len(),
        all_genes.len(),
        counts.x.len()
    );
    let settings = count::Settings {
        cell_roles: args.cell_roles,
        gene_roles: args.gene_roles.unwrap_or(8),
        ..count::Settings::default()
    };
    let options = super::vi::options(args);
    record.set_result("model", &settings);
    record.set_result("vi", &options);
    record.set_result(
        "cells",
        &serde_json::json!({"read": cells.len(), "kept": kept.len(), "fitted": kept.len()}),
    );
    record.set_result(
        "genes",
        &serde_json::json!({"read": all_genes.len(), "kept": names.len()}),
    );
    let (aa, bb) = (settings.cell_roles, settings.gene_roles);
    progress.check()?;
    // Each restart's steps: its k-means start's iterations, then its epochs.
    let warmup = vi::KMEANS_ITERS;
    progress.begin_vi(warmup + options.epochs, warmup);
    let (best, ends) = progress::barred(progress, || {
        // The profiles serve the starts alone: moved into theirs, dropped
        // with it once they are made.
        let profiles = vi::Profiles::new(&counts);
        let (counts, settings) = (&counts, &settings);
        super::vi::restarts(
            progress,
            (args.seed, args.chains),
            (warmup, options.epochs),
            move |rng, tick, stop| vi::kmeans_start(counts, &profiles, settings, rng, tick, stop),
            |start, rng, hooks| vi::fit(counts, settings, &options, start, rng, hooks),
            |f| End {
                elbo: f.elbo,
                per_unit: f.elbo_per_count(counts),
                epochs: f.epochs,
                svi_epochs: f.svi_epochs,
                converged: f.converged,
            },
        )
    })?;
    if progress.cancelled() {
        return Err(progress::Cancelled.into());
    }
    // Stopped before any restart finished an epoch: nothing to write, as
    // a cancel.
    let Some((best, fitted)) = best else {
        log::warn!("stopped before any restart finished an epoch: nothing to write");
        return Err(progress::Cancelled.into());
    };
    let ended = ends[best].expect("the best restart ended");
    let stopped = progress.stopped_early().is_some();
    let means = vi::summarize(&fitted);
    // The verdict, and the commit to writing: nothing is written before.
    super::vi::conclude(
        progress,
        record,
        ("count", "count"),
        &ends,
        (best, ended),
        (warmup, options.epochs, options.tol),
        stopped,
        serde_json::json!({"empty_roles": super::vi::empty_roles(&means.units, bb)}),
    )?;
    std::fs::create_dir_all(dir)?;
    let path = |name: &str| dir.join(format!("count_{name}.parquet"));
    let mut rows: Vec<write::RoleRow> = cells
        .iter()
        .map(|cell| write::RoleRow {
            cell,
            membership: None,
            level: None,
        })
        .collect();
    for (i, &c) in kept.iter().enumerate() {
        rows[c].membership = Some(&means.cell[i * aa..(i + 1) * aa]);
    }
    write::cell_roles(&rows, aa, &path("cell_roles"))?;
    write::gene_roles(&names, &means.gene, bb, &path("gene_roles"))?;
    write::blocks(
        &means.log_rate,
        &means.log_rate_sd,
        &means.units,
        bb,
        &path("blocks"),
    )?;
    // Each gene's UMIs in each cell module (each cell's counts by its
    // membership) and the modules' cells, for the module view.
    let mut umis = vec![0.0f64; names.len() * aa];
    let mut module_cells = vec![0.0f64; aa];
    for (i, &c) in kept.iter().enumerate() {
        let row = &restricted[c];
        let theta = &means.cell[i * aa..(i + 1) * aa];
        for (m, &t) in module_cells.iter_mut().zip(theta) {
            *m += t as f64;
        }
        for &(g, x) in row {
            for (a, &t) in theta.iter().enumerate() {
                umis[g as usize * aa + a] += x as f64 * t as f64;
            }
        }
    }
    let main: Vec<usize> = main_labels(&means.gene, bb)
        .into_iter()
        .map(usize::from)
        .collect();
    write::count_genes(&names, &main, &umis, &module_cells, &path("block_sites"))?;
    write::trace(&progress.traces(), &path("trace"))?;
    record.set_result(
        "tables",
        &serde_json::json!({
            "modality": COUNT_MODALITY,
            "cells": ["count_cell_roles.parquet"],
            "genes": ["count_gene_roles.parquet", "count_blocks.parquet", "count_block_sites.parquet"],
            "trace": "count_trace.parquet",
        }),
    );
    Ok(Written {
        modality: COUNT_MODALITY.into(),
        record: RECORD.into(),
    })
}
