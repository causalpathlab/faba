//! A kept draw: the chain's running posterior means, its held-out score and
//! its scalars.

use super::{ln_mean, Alignment, Model, State};
use crate::cluster::crp::{ln_sigmoid, log_sum_exp, logit};
use std::sync::{Arc, Mutex};

/// A chain's running sums over its kept draws, shared by its draws.
#[derive(Debug, Default)]
pub struct Tally {
    cell: Vec<f64>,
    gene: Vec<f64>,
    beta: Vec<f64>,
    shift: Vec<f64>,
    shift_sq: Vec<f64>,
    /// Each cell's level (at infinite strength only).
    level: Vec<f64>,
    /// Links per gene and block, per cell and block (`K` blocks per row),
    /// and reads per site and block.
    gene_block: Vec<f64>,
    cell_block: Vec<f64>,
    site_on: Vec<f64>,
    site_off: Vec<f64>,
    draws: usize,
}

/// A chain's posterior means over its kept draws: each cell's and gene's
/// role weights (row-major), the global weights, each block's shift with
/// its sd, each cell's level (empty at finite strength, where there is
/// none), and each gene's, cell's and site's links or reads in each block
/// (`K` blocks per row; the sites' only when asked for).
pub struct Means {
    pub cell: Vec<f32>,
    pub gene: Vec<f32>,
    pub beta: Vec<f64>,
    pub shift: Vec<f64>,
    pub shift_sd: Vec<f64>,
    pub level: Vec<f64>,
    pub gene_block: Vec<f64>,
    pub cell_block: Vec<f64>,
    pub site_on: Vec<f64>,
    pub site_off: Vec<f64>,
}

/// One kept draw: its scalars and its chain's shared tally.
#[derive(Debug, Clone)]
pub struct Draw {
    /// Held-out reads' log-likelihood and read count, at this draw.
    pub heldout: (f64, f64),
    pub log_joint: f64,
    /// Each block's shift and each cell role's mean membership over the
    /// cells, at this draw: the chains' convergence read from them.
    pub(super) shift: Vec<f64>,
    pub(super) role_mean: Vec<f64>,
    pub(super) tally: Arc<Mutex<Tally>>,
}

/// Each chain's log joint at each kept draw.
pub(super) fn log_joints(chains: &[Vec<Draw>]) -> Vec<Vec<f64>> {
    chains
        .iter()
        .map(|c| c.iter().map(|d| d.log_joint).collect())
        .collect()
}

impl Draw {
    /// The chain's means; with `sites`, also the sites' reads per block,
    /// the largest part of the tally.
    pub fn means(&self, sites: bool) -> Means {
        self.tally.lock().expect("a chain's tally").means(sites)
    }
}

impl Tally {
    /// The means over the tallied draws; with `sites`, also the sites'
    /// reads per block.
    pub(super) fn means(&self, sites: bool) -> Means {
        let n = self.draws.max(1) as f64;
        let mean = |v: &[f64]| -> Vec<f64> { v.iter().map(|x| x / n).collect() };
        let site = |v: &[f64]| if sites { mean(v) } else { Vec::new() };
        let shift = mean(&self.shift);
        Means {
            cell: self.cell.iter().map(|x| (x / n) as f32).collect(),
            gene: self.gene.iter().map(|x| (x / n) as f32).collect(),
            beta: mean(&self.beta),
            shift_sd: self
                .shift_sq
                .iter()
                .zip(&shift)
                .map(|(sq, m)| (sq / n - m * m).max(0.0).sqrt())
                .collect(),
            shift,
            level: mean(&self.level),
            gene_block: mean(&self.gene_block),
            cell_block: mean(&self.cell_block),
            site_on: site(&self.site_on),
            site_off: site(&self.site_off),
        }
    }

    /// `other`'s sums added, its roles renumbered into this tally's by
    /// `to`.
    pub(super) fn add_aligned(&mut self, other: &Tally, to: &Alignment) {
        let blocks = to.blocks();
        add_columns(&mut self.cell, &other.cell, &to.cell);
        add_columns(&mut self.gene, &other.gene, &to.gene);
        add_columns(&mut self.beta, &other.beta, &to.cell);
        add(&mut self.level, &other.level, |x| x);
        for (sums, x) in [
            (&mut self.shift, &other.shift),
            (&mut self.shift_sq, &other.shift_sq),
            (&mut self.gene_block, &other.gene_block),
            (&mut self.cell_block, &other.cell_block),
            (&mut self.site_on, &other.site_on),
            (&mut self.site_off, &other.site_off),
        ] {
            add_columns(sums, x, &blocks);
        }
        self.draws += other.draws;
    }
}

/// `sums` grown to `x`'s length, then `x` added, each row's column `j`
/// (`to.len()` columns per row) into column `to[j]`.
pub(super) fn add_columns(sums: &mut Vec<f64>, x: &[f64], to: &[usize]) {
    sums.resize(x.len(), 0.0);
    for (s, row) in sums.chunks_mut(to.len()).zip(x.chunks(to.len())) {
        row.iter().zip(to).for_each(|(&v, &j)| s[j] += v);
    }
}

/// Each block's counts summed over the rows (`blocks` per row).
fn column_totals(counts: &[f64], blocks: usize) -> Vec<f64> {
    let mut column = vec![0.0f64; blocks];
    for row in counts.chunks(blocks) {
        column.iter_mut().zip(row).for_each(|(c, &x)| *c += x);
    }
    column
}

/// Each row's share of its counts in each block over the block's share of
/// all counts: above 1, the row is enriched in the block (0 for an empty
/// row).
pub(super) fn enrichment(counts: &[f64], blocks: usize) -> Vec<f64> {
    let column = column_totals(counts, blocks);
    let total: f64 = column.iter().sum();
    counts
        .chunks(blocks)
        .flat_map(|row| {
            let n: f64 = row.iter().sum();
            row.iter().zip(&column).map(move |(&x, &c)| {
                if n > 0.0 && c > 0.0 {
                    (x / n) / (c / total)
                } else {
                    0.0
                }
            })
        })
        .collect()
}

impl Means {
    /// Blocks, `A x B`.
    pub fn blocks(&self) -> usize {
        self.shift.len()
    }

    /// Each block's links, over all genes.
    pub fn block_links(&self) -> Vec<f64> {
        column_totals(&self.gene_block, self.blocks())
    }

    /// Each gene's enrichment in each block.
    pub fn gene_enrichment(&self) -> Vec<f64> {
        enrichment(&self.gene_block, self.blocks())
    }

    /// Each cell's enrichment in each block.
    pub fn cell_enrichment(&self) -> Vec<f64> {
        enrichment(&self.cell_block, self.blocks())
    }

    /// Each site's methylation in each block over its baseline `mu`, on the
    /// logit scale (`NaN` where the block holds none of its reads).
    pub fn site_lift(&self, mu: &[f64]) -> Vec<f64> {
        self.site_on
            .iter()
            .zip(&self.site_off)
            .enumerate()
            .map(|(at, (&on, &off))| {
                if on + off <= 0.0 {
                    return f64::NAN;
                }
                let p = (on + 0.5) / (on + off + 1.0);
                logit(p) - mu[at / self.blocks()]
            })
            .collect()
    }
}

/// `sums` grown to `x`'s length, then `f(x)` added.
fn add(sums: &mut Vec<f64>, x: &[f64], f: fn(f64) -> f64) {
    sums.resize(x.len(), 0.0);
    sums.iter_mut().zip(x).for_each(|(s, &x)| *s += f(x));
}

impl Model<'_> {
    /// Row-major role weights given the counts: each row's share of its
    /// links in each role around `conc * base`.
    fn weights(counts: &[u32], conc: f64, base: &[f64]) -> Vec<f64> {
        counts
            .chunks(base.len())
            .flat_map(|row| {
                let total = row.iter().sum::<u32>() as f64 + conc;
                row.iter()
                    .zip(base)
                    .map(move |(&x, &b)| (x as f64 + conc * b) / total)
            })
            .collect()
    }

    /// Site `s`'s log fractions `(ln p, ln (1 - p))` in block `k` for cell
    /// `c`: the block's posterior-mean profile, or at infinite strength the
    /// baseline shifted by the block's shift and the cell's level.
    fn ln_fraction(&self, st: &State, c: usize, s: usize, k: usize) -> (f64, f64) {
        if self.finite() {
            let at = s * self.blocks() + k;
            ln_mean(st.prior[at], st.on[at], st.off[at])
        } else {
            let z = self.mu[s] + st.delta[k] + st.nu[c];
            (ln_sigmoid(z), ln_sigmoid(-z))
        }
    }

    /// Held-out reads' log-likelihood and count: each read at its cell's
    /// and gene's role weights mixed over the blocks, a converted and an
    /// unconverted read each a log-sum-exp, so neither is taken as 1 less
    /// the other. Each block's term is computed once, into a buffer.
    fn score_heldout(&self, st: &State, cell: &[f64], gene: &[f64]) -> (f64, f64) {
        let Some(held) = self.heldout else {
            return (0.0, 0.0);
        };
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let (mut on, mut off) = (vec![0.0f64; aa * bb], vec![0.0f64; aa * bb]);
        super::super::read::score(held, |c, j| {
            let g = self.gene_of[j] as usize;
            for k in 0..aa * bb {
                let (p, q) = self.ln_fraction(st, c, j, k);
                let w = (cell[c * aa + k / bb] * gene[g * bb + k % bb]).ln();
                (on[k], off[k]) = (w + p, w + q);
            }
            (
                log_sum_exp(on.iter().copied()),
                log_sum_exp(off.iter().copied()),
            )
        })
    }

    /// The draw at `st`, its chain's tally brought up to date.
    pub(super) fn draw(&self, st: &State) -> Draw {
        let kk = self.blocks();
        let cell = Self::weights(&st.n, self.settings.alpha, &st.beta);
        let gene = Self::weights(&st.m, self.settings.gene_alpha, &self.gene_base);
        {
            let mut t = st.tally.lock().expect("a chain's tally");
            let same = |x: f64| x;
            add(&mut t.cell, &cell, same);
            add(&mut t.gene, &gene, same);
            add(&mut t.beta, &st.beta, same);
            add(&mut t.shift, &st.delta, same);
            add(&mut t.shift_sq, &st.delta, |x| x * x);
            if !self.finite() {
                add(&mut t.level, &st.nu, same);
            }
            // Links per block counted from the labels, exactly; reads per
            // site and block as the state keeps them at finite strength,
            // else gathered from the links.
            t.gene_block.resize(self.genes() * kk, 0.0);
            t.cell_block.resize(self.links.cells * kk, 0.0);
            t.site_on.resize(self.mu.len() * kk, 0.0);
            t.site_off.resize(self.mu.len() * kk, 0.0);
            let finite = self.finite();
            if finite {
                add(&mut t.site_on, &st.on, same);
                add(&mut t.site_off, &st.off, same);
            }
            for l in 0..self.links.len() {
                let k = self.block(st, l);
                t.gene_block[self.link_gene[l] as usize * kk + k] += 1.0;
                t.cell_block[self.links.cell[l] as usize * kk + k] += 1.0;
                if finite {
                    continue;
                }
                let (s, on, off) = self.links.of(l);
                for i in 0..s.len() {
                    t.site_on[s[i] as usize * kk + k] += on[i];
                    t.site_off[s[i] as usize * kk + k] += off[i];
                }
            }
            t.draws += 1;
        }
        let cells = self.links.cells.max(1) as f64;
        let mut role_mean = vec![0.0f64; self.settings.cell_roles];
        for row in cell.chunks(role_mean.len()) {
            role_mean
                .iter_mut()
                .zip(row)
                .for_each(|(m, &w)| *m += w / cells);
        }
        Draw {
            heldout: self.score_heldout(st, &cell, &gene),
            log_joint: self.log_joint(st),
            shift: st.delta.clone(),
            role_mean,
            tally: Arc::clone(&st.tally),
        }
    }
}
