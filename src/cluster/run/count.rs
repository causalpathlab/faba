//! `faba cluster --modality count`: the count blockmodel's fit (see
//! [`crate::cluster::count`]) of gene count matrices, from reading them to
//! its tables.

use super::super::command::COUNT_MODALITY;
use super::super::count::{self, Counts};
use super::super::read::{walk_matrix, Visit};
use super::super::summary::{chain_scores, main_labels};
use super::super::{args, progress, read, write, Written};
use crate::run_record::RunRecord;
use data_beans::aux::feature_rows::split_count_row;
use legume_numeric::mcmc::engine::{McmcConfig, McmcModel};
use log::info;
use rand::rngs::SmallRng;
use rustc_hash::FxHashMap;
use std::path::Path;
use std::sync::{Arc, Mutex};

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

/// A chain's posterior means: each cell's and gene's role shares, each
/// block's rate (log) with its sd, each block's units, each cell's
/// background share.
#[derive(Default)]
struct Means {
    cell: Vec<f32>,
    gene: Vec<f32>,
    log_rate: Vec<f64>,
    log_rate_sd: Vec<f64>,
    units: Vec<f64>,
    background: Vec<f64>,
}

/// A chain's kept draws summed as they come (none kept whole).
#[derive(Default)]
struct Tally {
    draws: usize,
    cell: Vec<f64>,
    gene: Vec<f64>,
    log_rate: Vec<f64>,
    log_rate_sq: Vec<f64>,
    units: Vec<f64>,
    background: Vec<f64>,
}

/// `x` added into `sums`, grown to its length.
fn add(sums: &mut Vec<f64>, x: impl ExactSizeIterator<Item = f64>) {
    sums.resize(x.len(), 0.0);
    sums.iter_mut().zip(x).for_each(|(s, v)| *s += v);
}

/// Row-major shares of each row of `x` (`k` per row).
fn shares(x: &[f64], k: usize) -> impl ExactSizeIterator<Item = f64> + '_ {
    let totals: Vec<f64> = x
        .chunks(k)
        .map(|r| r.iter().sum::<f64>().max(1e-300))
        .collect();
    x.iter().enumerate().map(move |(i, v)| v / totals[i / k])
}

impl Tally {
    /// The draw at `st` (`aa`, `bb` roles) added.
    fn add(&mut self, st: &count::State, (aa, bb): (usize, usize)) {
        self.draws += 1;
        add(&mut self.cell, shares(&st.theta, aa));
        add(&mut self.gene, shares(&st.psi, bb));
        add(&mut self.log_rate, st.lambda.iter().map(|l| l.ln()));
        add(
            &mut self.log_rate_sq,
            st.lambda.iter().map(|l| l.ln().powi(2)),
        );
        add(&mut self.units, st.blocks.iter().map(|&u| u as f64));
        let units =
            |c: usize| st.n0[c] as f64 + st.n[c * aa..(c + 1) * aa].iter().sum::<u32>() as f64;
        add(
            &mut self.background,
            (0..st.n0.len()).map(|c| st.n0[c] as f64 / units(c).max(1.0)),
        );
    }

    fn means(&self) -> Means {
        let n = self.draws.max(1) as f64;
        let mean = |v: &[f64]| -> Vec<f64> { v.iter().map(|x| x / n).collect() };
        let log_rate = mean(&self.log_rate);
        Means {
            cell: self.cell.iter().map(|x| (x / n) as f32).collect(),
            gene: self.gene.iter().map(|x| (x / n) as f32).collect(),
            log_rate_sd: (self.log_rate_sq.iter().zip(&log_rate))
                .map(|(s, m)| (s / n - m * m).max(0.0).sqrt())
                .collect(),
            log_rate,
            units: mean(&self.units),
            background: mean(&self.background),
        }
    }
}

/// The fit's chains: the model's sweeps, each chain's kept draws tallied.
struct Fit<'a> {
    model: count::Model<'a>,
    roles: (usize, usize),
}

/// A chain: the model's state and its tally.
struct Chain {
    st: count::State,
    tally: Arc<Mutex<Tally>>,
}

impl McmcModel for Fit<'_> {
    type State = Chain;
    /// A kept draw's log-likelihood, and its chain's tally.
    type Sample = (f64, Arc<Mutex<Tally>>);
    type Result = (Vec<f64>, Means);

    fn init(&self, rng: &mut SmallRng) -> Chain {
        Chain {
            st: self.model.start(rng),
            tally: Arc::default(),
        }
    }

    fn sweep(&self, ch: &mut Chain, rng: &mut SmallRng) {
        self.model.sweep(&mut ch.st, rng);
    }

    fn collect(&self, ch: &Chain) -> Self::Sample {
        ch.tally
            .lock()
            .expect("a chain's tally")
            .add(&ch.st, self.roles);
        (self.model.log_likelihood(&ch.st), ch.tally.clone())
    }

    fn summarize(&self, draws: Vec<Self::Sample>) -> (Vec<f64>, Means) {
        let means = draws.first().map_or_else(Means::default, |(_, t)| {
            t.lock().expect("a chain's tally").means()
        });
        (draws.into_iter().map(|(ll, _)| ll).collect(), means)
    }
}

impl progress::Traced for Fit<'_> {
    /// The counts' log-likelihood: what a count fit is watched by.
    fn log_joint_of(&self, ch: &Chain) -> Option<f64> {
        Some(self.model.log_likelihood(&ch.st))
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
        roles: (aa, bb),
    };
    progress.check()?;
    let chains = progress::run(&fit, config, args.chains, progress)?;
    progress.note_stop(config, record)?;
    let lls: Vec<Vec<f64>> = chains.iter().map(|(ll, _)| ll.clone()).collect();
    let (chain_ll, rhat, best) = chain_scores(&lls);
    progress.summarised(rhat)?;
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
            "background_share": means.background.iter().sum::<f64>() / means.background.len().max(1) as f64,
        }),
    );
    let verdict = progress.verdict().expect("summarised");
    record.set_result("converged", &verdict.converged());
    Ok(Written {
        modality: COUNT_MODALITY.into(),
        record: "count_blockmodel.run.json".into(),
    })
}
