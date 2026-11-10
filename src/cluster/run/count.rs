//! `faba cluster --modality count`: the count blockmodel's fit (see
//! [`crate::cluster::count`]) of gene count matrices, from reading them to
//! its tables.

use super::super::command::COUNT_MODALITY;
use super::super::count::{self, vi, Counts};
use super::super::progress::advice;
use super::super::read::{walk_matrix, Visit};
use super::super::summary::main_labels;
use super::super::{args, progress, read, write, Written};
use crate::run_record::RunRecord;
use data_beans::aux::feature_rows::split_count_row;
use log::info;
use rand::rngs::SmallRng;
use rand::SeedableRng;
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
    let options = vi::Options {
        svi_epochs: args.svi_epochs,
        batch_size: args.batch_size,
        epochs: args.epochs,
        tol: args.tol,
        ..vi::Options::default()
    };
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
    let restarts = args.chains;
    // Each restart's steps: its k-means start's iterations, then its epochs.
    let warmup = vi::KMEANS_ITERS;
    progress.begin_vi(warmup + options.epochs, warmup);
    let (best, ends) = progress::barred(progress, || {
        restarts_fit(&counts, &settings, &options, args.seed, restarts, progress)
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
    let per: Vec<f64> = ends.iter().flatten().map(|e| e.per_count).collect();
    let spread = per.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        - per.iter().copied().fold(f64::INFINITY, f64::min);
    let stopped = progress.stopped_early().is_some();
    if stopped {
        let epochs: Vec<Option<usize>> =
            ends.iter().map(|e| e.as_ref().map(|e| e.epochs)).collect();
        log::warn!(
            "stopped early: the tables are from the epochs run so far ({epochs:?} by restart)"
        );
        record.set_result("stopped_early", &serde_json::json!({"epochs": epochs}));
    }
    progress.optimised(advice::Verdict::vi(
        fitted.converged,
        stopped,
        (fitted.epochs, options.epochs),
        per.len(),
        spread,
    ))?;
    let means = vi::summarize(&fitted);
    info!(
        "count: each restart's ELBO per count {:?}; tables from restart {best} ({} epochs: {})",
        per.iter()
            .map(|x| (x * 1e4).round() / 1e4)
            .collect::<Vec<_>>(),
        fitted.epochs,
        progress.verdict().expect("optimised").said(),
    );
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
    // Roles a fit by VI keeps but gives no units to.
    let empty = |units: Vec<f64>| -> Vec<usize> {
        (units.iter().enumerate())
            .filter(|(_, &u)| u < 1.0)
            .map(|(k, _)| k + 1)
            .collect()
    };
    let cell_units = (means.units.chunks_exact(bb))
        .map(|row| row.iter().sum())
        .collect();
    let gene_units = count::column_totals(&means.units, bb);
    record.set_result(
        "tables",
        &serde_json::json!({
            "modality": COUNT_MODALITY,
            "cells": ["count_cell_roles.parquet"],
            "genes": ["count_gene_roles.parquet", "count_blocks.parquet", "count_block_sites.parquet"],
            "trace": "count_trace.parquet",
        }),
    );
    record.set_result(
        "trace",
        &serde_json::json!({
            "sweeps": warmup + options.epochs,
            "burn_in": warmup,
            "quantity": progress.quantity(),
        }),
    );
    record.set_result(
        "count",
        &serde_json::json!({
            "vi": {
                "restarts": (ends.iter().enumerate()).map(|(i, e)| match e {
                    Some(e) => serde_json::json!({
                        "restart": i,
                        "elbo": e.elbo,
                        "elbo_per_count": e.per_count,
                        "epochs": e.epochs,
                        "svi_epochs": e.svi_epochs,
                        "converged": e.converged,
                    }),
                    None => serde_json::json!({"restart": i, "ran": false}),
                }).collect::<Vec<_>>(),
                "best": best,
                "converged": fitted.converged,
                "stopped_early": stopped,
                "epochs": fitted.epochs,
                "asked": options.epochs,
                "spread_per_count": spread,
                "tol": options.tol,
            },
            "empty_roles": {"cell": empty(cell_units), "gene": empty(gene_units)},
        }),
    );
    let verdict = progress.verdict().expect("summarised");
    record.set_result("converged", &verdict.converged());
    Ok(Written {
        modality: COUNT_MODALITY.into(),
        record: RECORD.into(),
    })
}

/// The best restart's fit and index, and each restart's end (see
/// [`restarts_fit`]).
type Restarts = (Option<(usize, vi::Fitted)>, Vec<Option<End>>);

/// How a restart ended, for the record.
struct End {
    elbo: f64,
    per_count: f64,
    epochs: usize,
    svi_epochs: usize,
    converged: bool,
}

/// Each restart's fit: the k-means starts side by side (each on its own
/// thread), then the variational fits one after another, each using every
/// thread. Restart `i`'s random numbers are seeded `seed + i`; its k-means
/// iterations, then its epochs, are counted on `progress` as chain `i` and
/// its ELBO traced there. The best restart's fit (by ELBO) and its index,
/// none when no restart finished an epoch; each restart's end, none for
/// one stopped before an epoch ended or whose start failed.
fn restarts_fit(
    counts: &Counts,
    settings: &count::Settings,
    options: &vi::Options,
    seed: u64,
    restarts: usize,
    progress: &progress::Progress,
) -> anyhow::Result<Restarts> {
    let stop = || progress.stopping();
    let cancelled = || progress.cancelled();
    // The profiles serve the starts alone: dropped once they are made.
    let profiles = vi::Profiles::new(counts);
    let starts: Vec<anyhow::Result<(vi::Start, SmallRng)>> = std::thread::scope(|s| {
        let profiles = &profiles;
        let handles: Vec<_> = (0..restarts)
            .map(|i| {
                s.spawn(move || {
                    let mut rng = SmallRng::seed_from_u64(seed.wrapping_add(i as u64));
                    progress.reach(i, 0);
                    let tick = |n: usize| progress.reach(i, n);
                    vi::kmeans_start(counts, profiles, settings, &mut rng, &tick, &stop)
                        .map(|start| (start, rng))
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("a k-means start does not panic"))
            .collect()
    });
    drop(profiles);
    // Every start done (its iterations may have ended early): a restart
    // waits at its first epoch for its turn.
    let warmup = vi::KMEANS_ITERS;
    (0..restarts).for_each(|i| progress.reach(i, warmup));
    let mut best: Option<(usize, vi::Fitted)> = None;
    let mut ends = Vec::with_capacity(restarts);
    let mut failed = None;
    for (i, start) in starts.into_iter().enumerate() {
        let (start, mut rng) = match start {
            Ok(start) => start,
            Err(e) => {
                log::warn!("restart {i}: {e:#}");
                failed.get_or_insert(e);
                ends.push(None);
                continue;
            }
        };
        if stop() {
            ends.push(None);
            continue;
        }
        let epoch = |e: usize, f: f64| {
            progress.reach(i, warmup + e);
            progress.note(i, warmup + e, f);
        };
        let hooks = vi::Hooks {
            epoch: &epoch,
            stop: &stop,
            cancelled: &cancelled,
        };
        let Some(fitted) = vi::fit(counts, settings, options, start, &mut rng, &hooks) else {
            ends.push(None);
            continue;
        };
        // A restart that settled before its last epoch has nothing left.
        if fitted.converged {
            progress.reach(i, warmup + options.epochs);
        }
        ends.push(Some(End {
            elbo: fitted.elbo,
            per_count: fitted.elbo_per_count(counts),
            epochs: fitted.epochs,
            svi_epochs: fitted.svi_epochs,
            converged: fitted.converged,
        }));
        // Only the best restart's factors are kept.
        if best.as_ref().is_none_or(|(_, b)| fitted.elbo > b.elbo) {
            best = Some((i, fitted));
        }
    }
    if let (None, Some(e), false) = (&best, failed, stop()) {
        return Err(e);
    }
    Ok((best, ends))
}
