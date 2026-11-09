//! `faba cluster --modality count`: the count blockmodel's fit (see
//! [`crate::cluster::count`]) of gene count matrices, from reading them to
//! its tables.

use super::super::count::{self, Counts};
use super::super::{args, progress, read, write, Written};
use crate::run_record::RunRecord;
use data_beans::aux::feature_rows::split_count_row;
use data_beans::hdf5_io::resolve_backend_file;
use data_beans::sparse_io::open_sparse_matrix;
use legume_numeric::mcmc::engine::{McmcConfig, McmcModel};
use log::info;
use rand::rngs::SmallRng;
use rustc_hash::FxHashMap;
use std::path::Path;

/// Cells read from a matrix at a time.
const COLUMN_BLOCK: usize = 4096;

/// The count matrices `args` name: each data file, or each `{batch}_count`
/// matrix of an output directory, with its batch.
fn count_inputs(args: &args::ClusterArgs) -> anyhow::Result<Vec<(Box<str>, Box<str>)>> {
    let mut out = Vec::new();
    for f in &args.data_files {
        match crate::run_record::output_dir_of(f) {
            Some(dir) => {
                let layout = crate::qc::layout::scan_input_dir(&dir.to_string_lossy())?;
                let found: Vec<_> = layout
                    .matrices
                    .iter()
                    .filter(|m| &*m.kind == "count")
                    .map(|m| (m.path.clone(), m.batch.clone()))
                    .collect();
                anyhow::ensure!(!found.is_empty(), "{} has no count matrices", dir.display());
                out.extend(found);
            }
            None => {
                let name = crate::qc::layout::file_name(f);
                let batch = crate::qc::layout::matrix_batch(&name).unwrap_or(name);
                out.push((f.clone(), batch));
            }
        }
    }
    Ok(out)
}

/// The genes' names, the cells, and each cell's `(gene, count)`s.
type CountsRead = (Vec<Box<str>>, Vec<read::Cell>, Vec<Vec<(u32, u32)>>);

/// Every cell's spliced counts over the genes of every input: the genes'
/// names, the cells (barcode and batch), and each cell's `(gene, count)`.
fn read_counts(
    inputs: &[(Box<str>, Box<str>)],
    progress: &progress::Progress,
) -> anyhow::Result<CountsRead> {
    let mut genes: Vec<Box<str>> = Vec::new();
    let mut index: FxHashMap<Box<str>, u32> = FxHashMap::default();
    let mut cells = Vec::new();
    let mut per_cell: Vec<Vec<(u32, u32)>> = Vec::new();
    for (path, batch) in inputs {
        progress.check()?;
        let (backend, resolved) = resolve_backend_file(path, None)?;
        let data = open_sparse_matrix(&resolved, &backend)?;
        let rows: Vec<Option<u32>> = data
            .row_names()?
            .iter()
            .map(|name| {
                let (gene, false) = split_count_row(name)? else {
                    return None;
                };
                let next = genes.len() as u32;
                Some(*index.entry(gene.into()).or_insert_with(|| {
                    genes.push(gene.into());
                    next
                }))
            })
            .collect();
        anyhow::ensure!(
            rows.iter().any(Option::is_some),
            "{path} holds no spliced count rows"
        );
        let barcodes = data.column_names()?;
        let first = cells.len();
        cells.extend(barcodes.iter().map(|b| read::Cell {
            barcode: b.clone(),
            batch: batch.clone(),
            fit: true,
        }));
        per_cell.resize_with(cells.len(), Vec::new);
        for start in (0..barcodes.len()).step_by(COLUMN_BLOCK) {
            progress.check()?;
            let block: Vec<usize> = (start..(start + COLUMN_BLOCK).min(barcodes.len())).collect();
            let (_, _, triplets) = data.read_triplets_by_columns(block)?;
            for (row, col, value) in triplets {
                let (Some(Some(g)), x) = (rows.get(row as usize), value.round().max(0.0) as u32)
                else {
                    continue;
                };
                if x > 0 {
                    per_cell[first + start + col as usize].push((*g, x));
                }
            }
        }
        info!("{path}: {} cells", barcodes.len());
    }
    Ok((genes, cells, per_cell))
}

/// The `n` genes of highest dispersion (variance over mean of the counts,
/// each cell's scaled to the median total) among those with any count,
/// in their order; all of them when there are fewer.
fn variable_genes(per_cell: &[Vec<(u32, u32)>], genes: usize, n: usize) -> Vec<usize> {
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
        for &(g, x) in cell {
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

/// A chain's posterior means: each cell's and gene's role shares, each
/// block's rate (log), each block's units, each cell's background share.
struct Means {
    cell: Vec<f32>,
    gene: Vec<f32>,
    log_rate: Vec<f64>,
    log_rate_sd: Vec<f64>,
    units: Vec<f64>,
    background: Vec<f64>,
}

/// Row-major shares of each row of `x` (`k` per row), as f32.
fn shares(x: &[f64], k: usize) -> impl Iterator<Item = f32> + '_ {
    x.chunks(k).flat_map(|row| {
        let total: f64 = row.iter().sum::<f64>().max(1e-300);
        row.iter().map(move |v| (v / total) as f32)
    })
}

/// The fit's draws: the model's, with each kept draw's shares, rates and
/// allocation summed.
struct Fit<'a> {
    model: count::Model<'a>,
    sizes: (usize, usize, usize, usize),
}

/// One kept draw.
struct Draw {
    log_likelihood: f64,
    cell: Vec<f32>,
    gene: Vec<f32>,
    log_rate: Vec<f64>,
    units: Vec<u32>,
    background: Vec<f64>,
}

impl McmcModel for Fit<'_> {
    type State = count::State;
    type Sample = Draw;
    type Result = (Vec<f64>, Means);

    fn init(&self, rng: &mut SmallRng) -> count::State {
        self.model.init(rng)
    }

    fn sweep(&self, st: &mut count::State, rng: &mut SmallRng) {
        self.model.sweep(st, rng);
    }

    fn collect(&self, st: &count::State) -> Draw {
        let (_, _, aa, bb) = self.sizes;
        let units_of =
            |c: usize| st.n0[c] as f64 + st.n[c * aa..(c + 1) * aa].iter().sum::<u32>() as f64;
        Draw {
            log_likelihood: self.model.log_likelihood(st),
            cell: shares(&st.theta, aa).collect(),
            gene: shares(&st.psi, bb).collect(),
            log_rate: st.lambda.iter().map(|l| l.ln()).collect(),
            units: st.blocks.clone(),
            background: (0..st.n0.len())
                .map(|c| st.n0[c] as f64 / units_of(c).max(1.0))
                .collect(),
        }
    }

    fn summarize(&self, draws: Vec<Draw>) -> (Vec<f64>, Means) {
        let n = draws.len().max(1) as f64;
        let mean32 = |f: &dyn Fn(&Draw) -> &[f32]| -> Vec<f32> {
            let mut out = vec![0.0f64; draws.first().map_or(0, |d| f(d).len())];
            for d in &draws {
                out.iter_mut().zip(f(d)).for_each(|(o, &x)| *o += x as f64);
            }
            out.iter().map(|x| (x / n) as f32).collect()
        };
        let mean64 = |f: &dyn Fn(&Draw) -> Vec<f64>| -> Vec<f64> {
            let mut out: Vec<f64> = Vec::new();
            for d in &draws {
                let x = f(d);
                out.resize(x.len(), 0.0);
                out.iter_mut().zip(&x).for_each(|(o, v)| *o += v);
            }
            out.iter().map(|x| x / n).collect()
        };
        let log_rate = mean64(&|d| d.log_rate.clone());
        let sq = mean64(&|d| d.log_rate.iter().map(|x| x * x).collect());
        let means = Means {
            cell: mean32(&|d| &d.cell),
            gene: mean32(&|d| &d.gene),
            log_rate_sd: sq
                .iter()
                .zip(&log_rate)
                .map(|(s, m)| (s - m * m).max(0.0).sqrt())
                .collect(),
            log_rate,
            units: mean64(&|d| d.units.iter().map(|&u| u as f64).collect()),
            background: mean64(&|d| d.background.clone()),
        };
        (draws.iter().map(|d| d.log_likelihood).collect(), means)
    }
}

impl progress::Traced for Fit<'_> {
    /// The counts' log-likelihood: what a count fit is watched by.
    fn log_joint_of(&self, st: &count::State) -> Option<f64> {
        Some(self.model.log_likelihood(st))
    }
}

/// Fit `args`' count matrices, write the tables in `dir` from the chain of
/// highest mean log-likelihood, and fill the record.
pub(in crate::cluster) fn count_fit(
    args: &args::ClusterArgs,
    config: &McmcConfig,
    dir: &Path,
    record: &mut RunRecord,
    progress: &progress::Progress,
) -> anyhow::Result<Written> {
    let inputs = count_inputs(args)?;
    let (all_genes, cells, per_cell) = read_counts(&inputs, progress)?;
    let chosen = variable_genes(&per_cell, all_genes.len(), args.count_genes);
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
    let kept: Vec<usize> = (0..cells.len())
        .filter(|&c| restricted[c].len() >= args.min_cell_sites)
        .collect();
    anyhow::ensure!(
        kept.len() >= 2 && !names.is_empty(),
        "{} cells and {} genes pass the filters: too few to cluster",
        kept.len(),
        names.len()
    );
    let kept_rows: Vec<Vec<(u32, u32)>> = kept.iter().map(|&c| restricted[c].clone()).collect();
    let counts = Counts::from_cells(names.len(), &kept_rows);
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
        gamma: args.community_gamma,
        ..count::Settings::default()
    };
    record.set_result("model", &settings);
    record.set_result(
        "cells",
        &serde_json::json!({"read": cells.len(), "kept": kept.len(), "fitted": kept.len()}),
    );
    record.set_result(
        "genes",
        &serde_json::json!({"read": all_genes.len(), "kept": names.len()}),
    );
    let (aa, bb) = (settings.cell_roles, settings.gene_roles);
    let fit = Fit {
        model: count::Model::new(&counts, settings)?,
        sizes: (kept.len(), names.len(), aa, bb),
    };
    progress.check()?;
    let chains = progress::run(&fit, config, args.chains, progress)?;
    let lls: Vec<Vec<f64>> = chains.iter().map(|(ll, _)| ll.clone()).collect();
    let rhat = super::super::summary::centred_rhat(&lls);
    progress.summarised(rhat)?;
    let mean = |v: &[f64]| v.iter().sum::<f64>() / v.len().max(1) as f64;
    let chain_ll: Vec<f64> = lls.iter().map(|v| mean(v)).collect();
    let best = (0..chains.len())
        .max_by(|&x, &y| chain_ll[x].total_cmp(&chain_ll[y]))
        .expect("a chain");
    let means = &chains[best].1;
    info!(
        "count: each chain's mean log-likelihood {:?}; split R-hat {rhat:.2}; tables from chain {best}",
        chain_ll.iter().map(|x| x.round()).collect::<Vec<_>>()
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
        // The `level` column: the cell's share of its units in the
        // background.
        rows[c].level = Some(means.background[i] as f32);
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
    for (i, row) in kept_rows.iter().enumerate() {
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
    let main: Vec<usize> = means
        .gene
        .chunks(bb)
        .map(|row| (0..bb).fold(0, |best, b| if row[b] > row[best] { b } else { best }))
        .collect();
    write::count_genes(&names, &main, &umis, &module_cells, &path("block_sites"))?;
    write::trace(&progress.traces(), &path("trace"))?;
    record.set_result(
        "tables",
        &serde_json::json!({
            "modality": "count",
            "cells": ["count_cell_roles.parquet"],
            "genes": ["count_gene_roles.parquet", "count_blocks.parquet", "count_block_sites.parquet"],
            "trace": "count_trace.parquet",
            "measure": "per_cell",
        }),
    );
    record.set_result(
        "trace",
        &serde_json::json!({
            "sweeps": progress::sweeps(config),
            "burn_in": config.warmup,
        }),
    );
    record.set_result(
        "count",
        &serde_json::json!({
            "chain_log_likelihood": chain_ll,
            "rhat_log_likelihood": rhat,
            "best_chain": best,
            "background_share": mean(&means.background),
        }),
    );
    let verdict = progress.verdict().expect("summarised");
    record.set_result("converged", &verdict.converged());
    Ok(Written {
        modality: "count".into(),
        record: "count_blockmodel.run.json".into(),
    })
}
