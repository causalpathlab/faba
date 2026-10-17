//! Overlapping cell and gene modules: a binomial two-role blockmodel over
//! (cell, gene) links, after Ball, Karrer & Newman's link communities and
//! the mixed-membership blockmodel.
//!
//! A cell has membership weights `theta_c` over `A` cell roles, a gene
//! `psi_g` over `B` gene roles. Each link (a cell's reads at one gene's
//! sites) takes a cell role `a ~ theta_c` and a gene role `b ~ psi_g`, and
//! its reads at site `s` are
//!
//! ```text
//! Binomial(n_cs, p),  logit p = mu_s + nu_c + delta[a, b]
//! ```
//!
//! with `mu_s` the site's pooled baseline, `nu_c ~ N(0, cell_sd^2)` the
//! cell's own level (so that roles carry patterns, not levels), and one
//! shift `delta[a, b] ~ N(0, shift_sd^2)` per block, shared by every gene
//! in role `b`: a gene module is the genes whose methylation shifts the
//! same way across the cell roles.
//! With `theta_c ~ Dirichlet(alpha / A)` and `psi_g ~ Dirichlet(gene_alpha
//! / B)` integrated out, a link's pair of roles is drawn from
//!
//! ```text
//! P(a, b | rest) ∝ (n_ca + alpha / A) (m_gb + gene_alpha / B) Π_s Binomial(m_s, u_s | p_sab)
//! ```
//!
//! (counts excluding the link): collapsed Gibbs, one link at a time. Whole
//! cells exchange two cell roles and whole genes two gene roles by
//! Metropolis-Hastings swaps, each its own inverse; each shift and each
//! cell's level is drawn from its conditional by slice sampling. Every move
//! leaves the posterior invariant.

use super::crp::{ln_sigmoid, ln_urn, log_sum_exp, pooled_fractions, sample_log};
use super::link::Links;
use super::read::Reads;
use legume_numeric::mcmc::engine::McmcModel;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::RngExt;
use std::sync::{Arc, Mutex};

#[derive(Debug, Clone)]
pub struct Settings {
    /// Cell roles `A`.
    pub cell_roles: usize,
    /// Gene roles `B`.
    pub gene_roles: usize,
    /// Dirichlet concentration of each cell's role weights.
    pub alpha: f64,
    /// Dirichlet concentration of each gene's role weights.
    pub gene_alpha: f64,
    /// Prior standard deviation of each block's shift.
    pub shift_sd: f64,
    /// Prior standard deviation of each cell's level.
    pub cell_sd: f64,
}

pub struct Model<'a> {
    links: &'a Links,
    /// Each link's gene.
    link_gene: Vec<u32>,
    /// Gene `g`'s links are `gene_link[gene_ptr[g]..gene_ptr[g + 1]]`.
    gene_ptr: Vec<usize>,
    gene_link: Vec<usize>,
    /// Each site's baseline, the logit of its pooled fraction.
    mu: Vec<f64>,
    settings: Settings,
    /// Held-out reads by cell, and each site's gene, scored at every draw.
    heldout: Option<(&'a Reads, &'a [u32])>,
}

/// A chain's labels and the counts they imply: links per cell and cell
/// role (`cells x A`), per gene and gene role (`genes x B`).
pub struct State {
    a: Vec<u16>,
    b: Vec<u16>,
    n: Vec<u32>,
    m: Vec<u32>,
    shift: Vec<f64>,
    /// Each cell's level.
    nu: Vec<f64>,
    tally: Arc<Mutex<Tally>>,
}

/// A chain's running sums over its kept draws, shared by its draws.
#[derive(Debug, Default)]
pub struct Tally {
    cell: Vec<f64>,
    gene: Vec<f64>,
    shift: Vec<f64>,
    shift_sq: Vec<f64>,
    nu: Vec<f64>,
    /// Links per gene and block, per cell and block (`A x B` blocks,
    /// row-major), and reads per site and block.
    gene_block: Vec<f64>,
    cell_block: Vec<f64>,
    site_on: Vec<f64>,
    site_off: Vec<f64>,
    draws: usize,
}

/// A chain's posterior means over its kept draws: each cell's and gene's
/// role weights (row-major), each block's shift with its sd, each cell's
/// level, and each gene's, cell's and site's links or reads in each block
/// (`A x B` blocks per row, row-major; the sites' only when asked for).
pub struct Means {
    pub cell: Vec<f32>,
    pub gene: Vec<f32>,
    pub shift: Vec<f64>,
    pub shift_sd: Vec<f64>,
    pub cell_base: Vec<f64>,
    pub gene_block: Vec<f64>,
    pub cell_block: Vec<f64>,
    pub site_on: Vec<f64>,
    pub site_off: Vec<f64>,
}

/// Each row's share of its counts in each block over the block's share of
/// all counts: above 1, the row is enriched in the block (0 for an empty
/// row).
fn enrichment(counts: &[f64], blocks: usize) -> Vec<f64> {
    let mut column = vec![0.0f64; blocks];
    for row in counts.chunks(blocks) {
        column.iter_mut().zip(row).for_each(|(c, &x)| *c += x);
    }
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
                (p / (1.0 - p)).ln() - mu[at / self.blocks()]
            })
            .collect()
    }
}

/// One kept draw: its scalars and its chain's shared tally.
#[derive(Debug, Clone)]
pub struct Draw {
    /// Held-out reads' log-likelihood and read count, at this draw.
    pub heldout: (f64, f64),
    pub log_joint: f64,
    tally: Arc<Mutex<Tally>>,
}

impl Draw {
    /// The chain's means; with `sites`, also the sites' reads per block,
    /// the largest part of the tally.
    pub fn means(&self, sites: bool) -> Means {
        let t = self.tally.lock().expect("a chain's tally");
        let n = t.draws.max(1) as f64;
        let mean = |v: &[f64]| -> Vec<f64> { v.iter().map(|x| x / n).collect() };
        let shift: Vec<f64> = t.shift.iter().map(|x| x / n).collect();
        Means {
            cell: t.cell.iter().map(|x| (x / n) as f32).collect(),
            gene: t.gene.iter().map(|x| (x / n) as f32).collect(),
            shift_sd: t
                .shift_sq
                .iter()
                .zip(&shift)
                .map(|(sq, m)| (sq / n - m * m).max(0.0).sqrt())
                .collect(),
            shift,
            cell_base: mean(&t.nu),
            gene_block: mean(&t.gene_block),
            cell_block: mean(&t.cell_block),
            site_on: if sites { mean(&t.site_on) } else { Vec::new() },
            site_off: if sites { mean(&t.site_off) } else { Vec::new() },
        }
    }
}

impl<'a> Model<'a> {
    /// `gene_of`: each site's gene, numbered from 0.
    pub fn new(links: &'a Links, gene_of: &[u32], sites: usize, settings: Settings) -> Self {
        let mu = pooled_fractions(&links.site, &links.on, &links.off, sites)
            .into_iter()
            .map(|p| (p / (1.0 - p)).ln())
            .collect();
        let link_gene = links.genes_of(gene_of);
        let genes = gene_of.iter().max().map_or(0, |&g| g as usize + 1);
        let mut gene_ptr = vec![0usize; genes + 1];
        for &g in &link_gene {
            gene_ptr[g as usize + 1] += 1;
        }
        for g in 0..genes {
            gene_ptr[g + 1] += gene_ptr[g];
        }
        let mut next = gene_ptr.clone();
        let mut gene_link = vec![0usize; link_gene.len()];
        for (l, &g) in link_gene.iter().enumerate() {
            gene_link[next[g as usize]] = l;
            next[g as usize] += 1;
        }
        Self {
            links,
            link_gene,
            gene_ptr,
            gene_link,
            mu,
            settings,
            heldout: None,
        }
    }

    /// Score `heldout` (each cell's held-out reads; `gene_of`, each site's
    /// gene) at every kept draw.
    pub fn with_heldout(mut self, heldout: &'a Reads, gene_of: &'a [u32]) -> Self {
        self.heldout = Some((heldout, gene_of));
        self
    }

    fn genes(&self) -> usize {
        self.gene_ptr.len() - 1
    }

    /// The sites' baselines.
    pub fn baseline(&self) -> &[f64] {
        &self.mu
    }

    /// Link `l`'s reads' log-likelihood at shift `shift` on top of its
    /// cell's level.
    fn ln_link(&self, st: &State, l: usize, shift: f64) -> f64 {
        ln_reads(
            &self.mu,
            self.links,
            l,
            shift + st.nu[self.links.cell[l] as usize],
        )
    }

    fn put(&self, st: &mut State, l: usize, sign: i64) {
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let (c, g) = (self.links.cell[l] as usize, self.link_gene[l] as usize);
        let (n, m) = (
            &mut st.n[c * aa + st.a[l] as usize],
            &mut st.m[g * bb + st.b[l] as usize],
        );
        *n = (*n as i64 + sign) as u32;
        *m = (*m as i64 + sign) as u32;
    }

    /// Link `l`'s pair of roles from its conditional.
    fn relabel(&self, st: &mut State, l: usize, w: &mut [f64], rng: &mut SmallRng) {
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let (c, g) = (self.links.cell[l] as usize, self.link_gene[l] as usize);
        self.put(st, l, -1);
        let (ca, gb) = (
            self.settings.alpha / aa as f64,
            self.settings.gene_alpha / bb as f64,
        );
        for a in 0..aa {
            let cell = (st.n[c * aa + a] as f64 + ca).ln();
            for b in 0..bb {
                w[a * bb + b] = cell + self.ln_link(st, l, st.shift[a * bb + b]);
            }
        }
        // The gene's urn, one log per gene role.
        for b in 0..bb {
            let gene = (st.m[g * bb + b] as f64 + gb).ln();
            (0..aa).for_each(|a| w[a * bb + b] += gene);
        }
        let k = sample_log(w, rng);
        st.a[l] = (k / bb) as u16;
        st.b[l] = (k % bb) as u16;
        self.put(st, l, 1);
    }

    /// A swap of two cell roles on all of cell `c`'s links (`cell` true),
    /// or of two gene roles on all of gene `c`'s links: its own inverse, so
    /// the proposal is symmetric. The urn terms are symmetric in the roles,
    /// so the Metropolis-Hastings ratio is the change in the reads'
    /// likelihood alone.
    fn swap(&self, st: &mut State, c: usize, cell: bool, rng: &mut SmallRng) {
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let roles = if cell { aa } else { bb };
        if roles < 2 {
            return;
        }
        let r = rng.random_range(0..roles) as u16;
        let r2 = ((r as usize + rng.random_range(1..roles)) % roles) as u16;
        let links: Vec<usize> = if cell {
            (self.links.cell_ptr[c]..self.links.cell_ptr[c + 1]).collect()
        } else {
            self.gene_link[self.gene_ptr[c]..self.gene_ptr[c + 1]].to_vec()
        };
        let role = |st: &State, l: usize| if cell { st.a[l] } else { st.b[l] };
        let other = |x: u16| if x == r { r2 } else { r };
        let block = |a: u16, b: u16| a as usize * bb + b as usize;
        let mut delta = 0.0;
        let mut moved = Vec::new();
        for &l in &links {
            let x = role(st, l);
            if x != r && x != r2 {
                continue;
            }
            let (from, to) = if cell {
                (block(x, st.b[l]), block(other(x), st.b[l]))
            } else {
                (block(st.a[l], x), block(st.a[l], other(x)))
            };
            delta += self.ln_link(st, l, st.shift[to]) - self.ln_link(st, l, st.shift[from]);
            moved.push(l);
        }
        if moved.is_empty() || rng.random::<f64>().ln() >= delta {
            return;
        }
        for &l in &moved {
            self.put(st, l, -1);
            if cell {
                st.a[l] = other(st.a[l]);
            } else {
                st.b[l] = other(st.b[l]);
            }
            self.put(st, l, 1);
        }
    }

    /// Each block's shift from its conditional (its Normal prior times its
    /// links' Binomial reads), by a slice sampler with stepping out and
    /// shrinkage (Neal 2003), which leaves the conditional invariant.
    fn resample_shifts(&self, st: &mut State, rng: &mut SmallRng) {
        let bb = self.settings.gene_roles;
        let mut members: Vec<Vec<usize>> = vec![Vec::new(); st.shift.len()];
        for l in 0..self.links.len() {
            members[st.a[l] as usize * bb + st.b[l] as usize].push(l);
        }
        let var = self.settings.shift_sd * self.settings.shift_sd;
        for (k, links) in members.iter().enumerate() {
            let target = |d: f64| {
                -d * d / (2.0 * var) + links.iter().map(|&l| self.ln_link(st, l, d)).sum::<f64>()
            };
            st.shift[k] = slice(st.shift[k], &target, 1.0, rng);
        }
    }

    /// Each cell's level from its conditional (its Normal prior times its
    /// links' Binomial reads at their blocks' shifts), by slice sampling.
    fn resample_levels(&self, st: &mut State, rng: &mut SmallRng) {
        let bb = self.settings.gene_roles;
        let var = self.settings.cell_sd * self.settings.cell_sd;
        for c in 0..self.links.cells {
            let own = self.links.cell_ptr[c]..self.links.cell_ptr[c + 1];
            let target = |v: f64| {
                -v * v / (2.0 * var)
                    + own
                        .clone()
                        .map(|l| {
                            let shift = st.shift[st.a[l] as usize * bb + st.b[l] as usize];
                            ln_reads(&self.mu, self.links, l, shift + v)
                        })
                        .sum::<f64>()
            };
            st.nu[c] = slice(st.nu[c], &target, 1.0, rng);
        }
    }

    /// Log joint of the labels, the shifts and the reads (role weights
    /// integrated out).
    fn log_joint(&self, st: &State) -> f64 {
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let (cell_base, gene_base) = (vec![1.0 / aa as f64; aa], vec![1.0 / bb as f64; bb]);
        let cells: f64 =
            st.n.chunks(aa)
                .map(|r| ln_urn(r, self.settings.alpha, &cell_base))
                .sum();
        let genes: f64 =
            st.m.chunks(bb)
                .map(|r| ln_urn(r, self.settings.gene_alpha, &gene_base))
                .sum();
        let reads: f64 = (0..self.links.len())
            .map(|l| self.ln_link(st, l, st.shift[st.a[l] as usize * bb + st.b[l] as usize]))
            .sum();
        let normal =
            |x: &[f64], sd: f64| -> f64 { x.iter().map(|d| -d * d / (2.0 * sd * sd)).sum() };
        cells
            + genes
            + reads
            + normal(&st.shift, self.settings.shift_sd)
            + normal(&st.nu, self.settings.cell_sd)
    }

    /// Row-major role weights, posterior mean given the counts.
    fn weights(counts: &[u32], roles: usize, conc: f64) -> Vec<f64> {
        counts
            .chunks(roles)
            .flat_map(|row| {
                let total = row.iter().sum::<u32>() as f64 + conc;
                row.iter()
                    .map(move |&x| (x as f64 + conc / roles as f64) / total)
            })
            .collect()
    }

    /// Held-out reads' log-likelihood and count: each read at its cell's
    /// and gene's role weights mixed over the blocks.
    /// The log probabilities of a converted and of an unconverted read are
    /// each a log-sum-exp over the blocks, so neither is taken as 1 less the
    /// other.
    fn score_heldout(&self, st: &State, cell: &[f64], gene: &[f64]) -> (f64, f64) {
        let Some((held, gene_of)) = self.heldout else {
            return (0.0, 0.0);
        };
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        super::read::score(held, |c, j| {
            let g = gene_of[j] as usize;
            let terms = |sign: f64| {
                (0..aa * bb).map(move |k| {
                    let (a, b) = (k / bb, k % bb);
                    let x = self.mu[j] + st.nu[c] + st.shift[k];
                    (cell[c * aa + a] * gene[g * bb + b]).ln() + ln_sigmoid(sign * x)
                })
            };
            (log_sum_exp(terms(1.0)), log_sum_exp(terms(-1.0)))
        })
    }

    /// The labels after each of `draws` sweeps of the label moves alone
    /// (link relabels, then cell and gene swaps) after `warmup`, the shifts
    /// held at `shifts` and the cell levels at `levels`.
    #[cfg(test)]
    pub fn label_draws(
        &self,
        shifts: &[f64],
        levels: &[f64],
        warmup: usize,
        draws: usize,
        seed: u64,
    ) -> Vec<(Vec<u16>, Vec<u16>)> {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut st = self.init(&mut rng);
        st.shift = shifts.to_vec();
        st.nu = levels.to_vec();
        let mut w = vec![0.0f64; shifts.len()];
        (0..warmup + draws)
            .filter_map(|sweep| {
                let mut order: Vec<usize> = (0..self.links.len()).collect();
                order.shuffle(&mut rng);
                for l in order {
                    self.relabel(&mut st, l, &mut w, &mut rng);
                }
                for c in 0..self.links.cells {
                    self.swap(&mut st, c, true, &mut rng);
                }
                for g in 0..self.genes() {
                    self.swap(&mut st, g, false, &mut rng);
                }
                (sweep >= warmup).then(|| (st.a.clone(), st.b.clone()))
            })
            .collect()
    }

    /// Cell 0's level after each of `draws` level updates, the labels held
    /// at `a`, `b` and the shifts at `shifts` (to check the update against
    /// its conditional).
    #[cfg(test)]
    pub fn baseline_draws(
        &self,
        a: &[u16],
        b: &[u16],
        shifts: &[f64],
        draws: usize,
        seed: u64,
    ) -> Vec<f64> {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut st = self.state(a.to_vec(), b.to_vec());
        st.shift = shifts.to_vec();
        (0..draws)
            .map(|_| {
                self.resample_levels(&mut st, &mut rng);
                st.nu[0]
            })
            .collect()
    }

    /// `delta[0][0]` after each of `draws` shift updates, the labels held
    /// at `a`, `b` (to check the update against its conditional).
    #[cfg(test)]
    pub fn shift_draws(&self, a: &[u16], b: &[u16], draws: usize, seed: u64) -> Vec<f64> {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut st = self.state(a.to_vec(), b.to_vec());
        (0..draws)
            .map(|_| {
                self.resample_shifts(&mut st, &mut rng);
                st.shift[0]
            })
            .collect()
    }

    fn state(&self, a: Vec<u16>, b: Vec<u16>) -> State {
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let mut st = State {
            a,
            b,
            n: vec![0; self.links.cells * aa],
            m: vec![0; self.genes() * bb],
            shift: vec![0.0; aa * bb],
            nu: vec![0.0; self.links.cells],
            tally: Arc::default(),
        };
        for l in 0..self.links.len() {
            self.put(&mut st, l, 1);
        }
        st
    }
}

/// One slice-sampling update of `x` under log density `f` (Neal 2003):
/// a level under `f(x)`, an interval of width `w` placed at random around
/// `x` and stepped out until both ends fall below the level, then points
/// drawn from it, shrinking it toward `x`, until one is above the level.
fn slice(x: f64, f: &dyn Fn(f64) -> f64, w: f64, rng: &mut SmallRng) -> f64 {
    let level = f(x) + (1.0 - rng.random::<f64>()).ln();
    let mut lo = x - w * rng.random::<f64>();
    let mut hi = lo + w;
    while f(lo) > level {
        lo -= w;
    }
    while f(hi) > level {
        hi += w;
    }
    loop {
        let y = lo + (hi - lo) * rng.random::<f64>();
        if f(y) > level {
            return y;
        }
        if y < x {
            lo = y;
        } else {
            hi = y;
        }
    }
}

impl McmcModel for Model<'_> {
    type State = State;
    type Sample = Draw;
    type Result = Vec<Draw>;

    /// Each cell's links in one cell role, each gene's in one gene role,
    /// drawn at random.
    fn init(&self, rng: &mut SmallRng) -> State {
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let cell: Vec<u16> = (0..self.links.cells)
            .map(|_| rng.random_range(0..aa) as u16)
            .collect();
        let gene: Vec<u16> = (0..self.genes())
            .map(|_| rng.random_range(0..bb) as u16)
            .collect();
        let a = self.links.cell.iter().map(|&c| cell[c as usize]).collect();
        let b = self.link_gene.iter().map(|&g| gene[g as usize]).collect();
        self.state(a, b)
    }

    /// Every link relabelled once in random order, then a swap per cell and
    /// per gene, then the shifts and the cells' levels.
    fn sweep(&self, st: &mut State, rng: &mut SmallRng) {
        let mut order: Vec<usize> = (0..self.links.len()).collect();
        order.shuffle(rng);
        let mut w = vec![0.0f64; self.settings.cell_roles * self.settings.gene_roles];
        for l in order {
            self.relabel(st, l, &mut w, rng);
        }
        for c in 0..self.links.cells {
            self.swap(st, c, true, rng);
        }
        for g in 0..self.genes() {
            self.swap(st, g, false, rng);
        }
        self.resample_shifts(st, rng);
        self.resample_levels(st, rng);
    }

    fn collect(&self, st: &State) -> Draw {
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let cell = Self::weights(&st.n, aa, self.settings.alpha);
        let gene = Self::weights(&st.m, bb, self.settings.gene_alpha);
        {
            let mut t = st.tally.lock().expect("a chain's tally");
            let add = |sums: &mut Vec<f64>, x: &[f64], f: fn(f64) -> f64| {
                sums.resize(x.len(), 0.0);
                sums.iter_mut().zip(x).for_each(|(s, &x)| *s += f(x));
            };
            let same = |x: f64| x;
            add(&mut t.cell, &cell, same);
            add(&mut t.gene, &gene, same);
            add(&mut t.shift, &st.shift, same);
            add(&mut t.shift_sq, &st.shift, |x| x * x);
            add(&mut t.nu, &st.nu, same);
            let k = aa * bb;
            let sites = self.mu.len();
            t.gene_block.resize(self.genes() * k, 0.0);
            t.cell_block.resize(self.links.cells * k, 0.0);
            t.site_on.resize(sites * k, 0.0);
            t.site_off.resize(sites * k, 0.0);
            for l in 0..self.links.len() {
                let block = st.a[l] as usize * bb + st.b[l] as usize;
                t.gene_block[self.link_gene[l] as usize * k + block] += 1.0;
                t.cell_block[self.links.cell[l] as usize * k + block] += 1.0;
                let (s, m, u) = self.links.of(l);
                for i in 0..s.len() {
                    t.site_on[s[i] as usize * k + block] += m[i];
                    t.site_off[s[i] as usize * k + block] += u[i];
                }
            }
            t.draws += 1;
        }
        Draw {
            heldout: self.score_heldout(st, &cell, &gene),
            log_joint: self.log_joint(st),
            tally: Arc::clone(&st.tally),
        }
    }

    fn summarize(&self, draws: Vec<Draw>) -> Vec<Draw> {
        draws
    }
}

/// `links`' link `l` reads' log-likelihood at shift `shift` from baselines
/// `mu`: `m ln sigmoid(x) + u ln sigmoid(-x)`, written with
/// `ln sigmoid(-x) = ln sigmoid(x) - x` for one transcendental per entry.
fn ln_reads(mu: &[f64], links: &Links, l: usize, shift: f64) -> f64 {
    let (s, m, u) = links.of(l);
    (0..s.len())
        .map(|i| {
            let x = mu[s[i] as usize] + shift;
            (m[i] + u[i]) * ln_sigmoid(x) - u[i] * x
        })
        .sum()
}

impl Model<'_> {
    /// Cells outside the fit, placed against a fitted chain's `means`: each
    /// cell's links relabelled by collapsed Gibbs over its cell roles (the
    /// gene roles summed out at the genes' mean weights, the shifts at their
    /// means) alternating with its level by slice sampling, and its role
    /// weights and level averaged over the sweeps after a burn-in. The
    /// cell's own urn and level are the model's; the gene weights and shifts
    /// are plug-in estimates, not integrated over. Each cell's role weights
    /// (`cells x A`, row-major) and level.
    pub fn place(
        &self,
        links: &Links,
        gene_of: &[u32],
        means: &Means,
        seed: u64,
    ) -> (Vec<f32>, Vec<f64>) {
        use super::link::{cell_rng, PLACE_BURN_IN, PLACE_SWEEPS};
        use rayon::prelude::*;
        let (aa, bb) = (self.settings.cell_roles, self.settings.gene_roles);
        let ca = self.settings.alpha / aa as f64;
        let var = self.settings.cell_sd * self.settings.cell_sd;
        let link_gene = links.genes_of(gene_of);
        let ln_gene: Vec<f64> = means
            .gene
            .iter()
            .map(|&w| (w as f64).max(1e-300).ln())
            .collect();
        // Link `l`'s log-likelihood in cell role `a` at level `v`, its gene
        // role summed out.
        let like = |l: usize, a: usize, v: f64| {
            let g = link_gene[l] as usize;
            log_sum_exp((0..bb).map(|b| {
                ln_gene[g * bb + b] + ln_reads(&self.mu, links, l, means.shift[a * bb + b] + v)
            }))
        };
        let (weights, levels): (Vec<Vec<f32>>, Vec<f64>) = (0..links.cells)
            .into_par_iter()
            .map(|c| {
                let mut rng = cell_rng(seed, c);
                let own: Vec<usize> = (links.cell_ptr[c]..links.cell_ptr[c + 1]).collect();
                let mut row = vec![0u32; aa];
                let mut z: Vec<usize> = own.iter().map(|_| rng.random_range(0..aa)).collect();
                z.iter().for_each(|&a| row[a] += 1);
                let total = z.len() as f64 + self.settings.alpha;
                let (mut v, mut level) = (0.0f64, 0.0f64);
                let mut weights = vec![0.0f64; aa];
                let mut w = vec![0.0f64; aa];
                for sweep in 0..PLACE_SWEEPS {
                    for (i, &l) in own.iter().enumerate() {
                        row[z[i]] -= 1;
                        for (a, x) in w.iter_mut().enumerate() {
                            *x = (row[a] as f64 + ca).ln() + like(l, a, v);
                        }
                        z[i] = sample_log(&mut w, &mut rng);
                        row[z[i]] += 1;
                    }
                    let target = |x: f64| {
                        -x * x / (2.0 * var)
                            + own
                                .iter()
                                .zip(&z)
                                .map(|(&l, &a)| like(l, a, x))
                                .sum::<f64>()
                    };
                    v = slice(v, &target, 1.0, &mut rng);
                    if sweep >= PLACE_BURN_IN {
                        for (a, x) in weights.iter_mut().enumerate() {
                            *x += (row[a] as f64 + ca) / total;
                        }
                        level += v;
                    }
                }
                let n = (PLACE_SWEEPS - PLACE_BURN_IN) as f64;
                (
                    weights.into_iter().map(|x| (x / n) as f32).collect(),
                    level / n,
                )
            })
            .unzip();
        (weights.concat(), levels)
    }
}

/// What the chains say, together.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Summary {
    pub cell_roles: usize,
    pub gene_roles: usize,
    pub chains: usize,
    pub draws: usize,
    /// Split R-hat of the log joint across chains.
    pub rhat_log_joint: f64,
    /// Each chain's mean log joint.
    pub chain_log_joint: Vec<f64>,
    /// Adjusted Rand index between each pair of chains' cells' main cell
    /// roles, and genes' main gene roles (posterior-mean weights).
    pub cell_agreement: Vec<f64>,
    pub gene_agreement: Vec<f64>,
    /// The best chain (highest mean log joint) and its posterior-mean block
    /// shifts (`A x B`, row-major) with their sds.
    pub best_chain: usize,
    pub shifts: Vec<f64>,
    pub shift_sd: Vec<f64>,
    /// Held-out reads' log-likelihood per read, over every kept draw.
    pub heldout_per_read: Option<f64>,
    /// Each chain's strong gene roles (whose shift ranges at least
    /// [`STRONG_SHIFT`] across the cell roles), and for each pair of chains
    /// the mean overlap of a strong role's top genes with the best-matching
    /// strong role's in the other chain.
    pub strong_gene_roles: Vec<usize>,
    pub top_gene_overlap: Vec<f64>,
}

/// A gene role is strong when its shift ranges at least this much (logit
/// scale) across the cell roles.
pub const STRONG_SHIFT: f64 = 0.5;
/// Top genes per strong gene role, compared across chains.
const TOP_GENES: usize = 50;
/// Links a gene needs in a chain to be ranked.
const MIN_LINKS: f64 = 5.0;

/// Each strong gene role's top genes in a chain: the [`TOP_GENES`] with at
/// least [`MIN_LINKS`] links most enriched in the role.
fn strong_top_genes(m: &Means, aa: usize, bb: usize) -> Vec<Vec<usize>> {
    let k = aa * bb;
    let genes = m.gene_block.len() / k.max(1);
    // Links per gene and gene role, cell roles summed.
    let mut role = vec![0.0f64; genes * bb];
    for (at, &x) in m.gene_block.iter().enumerate() {
        role[(at / k) * bb + (at % k) % bb] += x;
    }
    let enrich = enrichment(&role, bb);
    (0..bb)
        .filter(|&b| {
            let column: Vec<f64> = (0..aa).map(|a| m.shift[a * bb + b]).collect();
            let (lo, hi) = column
                .iter()
                .fold((f64::INFINITY, f64::NEG_INFINITY), |(l, h), &x| {
                    (l.min(x), h.max(x))
                });
            hi - lo >= STRONG_SHIFT
        })
        .map(|b| {
            let mut ranked: Vec<usize> = (0..genes)
                .filter(|&g| role[g * bb..(g + 1) * bb].iter().sum::<f64>() >= MIN_LINKS)
                .collect();
            ranked.sort_by(|&x, &y| enrich[y * bb + b].total_cmp(&enrich[x * bb + b]));
            ranked.truncate(TOP_GENES);
            ranked
        })
        .collect()
}

/// The mean, over `a`'s top-gene sets, of the best overlap (shared share)
/// with any of `b`'s (0 when either has none).
fn best_overlap(a: &[Vec<usize>], b: &[Vec<usize>]) -> f64 {
    if a.is_empty() || b.is_empty() {
        return 0.0;
    }
    let shared = |x: &[usize], y: &[usize]| {
        x.iter().filter(|g| y.contains(g)).count() as f64 / x.len().max(1) as f64
    };
    a.iter()
        .map(|x| b.iter().map(|y| shared(x, y)).fold(0.0, f64::max))
        .sum::<f64>()
        / a.len() as f64
}

pub fn summarize(model: &Model, chains: &[Vec<Draw>]) -> Summary {
    use super::summary::{chain_scores, main_labels, pairwise, pairwise_ari, per_read};
    let (aa, bb) = (model.settings.cell_roles, model.settings.gene_roles);
    let log_joint: Vec<Vec<f64>> = chains
        .iter()
        .map(|c| c.iter().map(|d| d.log_joint).collect())
        .collect();
    let (chain_log_joint, rhat_log_joint, best_chain) = chain_scores(&log_joint);
    let means: Vec<Means> = chains
        .iter()
        .filter_map(|c| c.last().map(|d| d.means(false)))
        .collect();
    let tops: Vec<Vec<Vec<usize>>> = means.iter().map(|m| strong_top_genes(m, aa, bb)).collect();
    let best = means.get(best_chain);
    Summary {
        cell_roles: aa,
        gene_roles: bb,
        chains: chains.len(),
        draws: chains.iter().map(Vec::len).sum(),
        rhat_log_joint,
        chain_log_joint,
        cell_agreement: pairwise_ari(
            &means
                .iter()
                .map(|m| main_labels(&m.cell, aa))
                .collect::<Vec<_>>(),
        ),
        gene_agreement: pairwise_ari(
            &means
                .iter()
                .map(|m| main_labels(&m.gene, bb))
                .collect::<Vec<_>>(),
        ),
        best_chain,
        shifts: best.map_or_else(Vec::new, |m| m.shift.clone()),
        shift_sd: best.map_or_else(Vec::new, |m| m.shift_sd.clone()),
        heldout_per_read: per_read(chains.iter().flatten().map(|d| d.heldout)),
        strong_gene_roles: tops.iter().map(Vec::len).collect(),
        top_gene_overlap: pairwise(&tops, |x, y| best_overlap(x, y)),
    }
}
