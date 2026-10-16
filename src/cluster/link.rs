//! Link communities: cells' reads grouped by (cell, gene) link, each link in
//! one of `K` communities.
//!
//! A community `k` is a methylation profile over sites, `p_sk`. A cell has
//! membership weights `theta_c` over the communities. Each of its links (its
//! reads at one gene's sites) takes a community `z ~ theta_c`, and its reads
//! at site `s` are `Binomial(n_cs, p_s,z)`. With `theta_c ~ Dirichlet(alpha)`
//! and `p_sk ~ Beta(A p0_s, A (1 - p0_s))` integrated out, a link's
//! community is drawn from
//!
//! ```text
//! P(z = k | rest) ∝ (n_ck + alpha) · Π_{s in link} BetaBinomial(m_s, u_s | site s's reads in k)
//! ```
//!
//! with every count excluding the link itself: collapsed Gibbs, one link at
//! a time, each move in detailed balance.
//!
//! Memberships share global weights (a weak-limit hierarchical Dirichlet
//! process): `beta ~ Dirichlet(gamma / K, ...)`, `theta_c ~ Dirichlet(alpha
//! beta)`, so a link's prior weight is `n_ck + alpha beta_k`. Within a cell
//! that is a Pólya urn; across cells, a community many cells use gets a
//! larger `beta_k` and pulls every cell, so the reads condense (rich get
//! richer) instead of spreading over equal communities. `beta` is drawn from
//! its conditional through each cell's table counts (Teh et al. 2006). A
//! cell swap moves a whole cell between communities, which single-link
//! moves cannot do against the cell's own urn.

use super::crp::{ln_beta, ln_pred, ln_urn, pooled_fractions, sample_log};
use super::read::{Reads, Weighted};
use legume_numeric::mcmc::engine::McmcModel;
use rand::rngs::SmallRng;
use rand::seq::SliceRandom;
use rand::RngExt;
use rand_distr::Distribution;
use std::sync::{Arc, Mutex};

/// Reads grouped into (cell, gene) links: link `l` is cell `cell[l]`'s
/// entries `ptr[l]..ptr[l + 1]` (sites and weighted reads) at one gene.
pub struct Links {
    /// Cell `c`'s links are `cell_ptr[c]..cell_ptr[c + 1]`.
    pub cell_ptr: Vec<usize>,
    pub cell: Vec<u32>,
    pub ptr: Vec<usize>,
    pub site: Vec<u32>,
    pub on: Vec<f64>,
    pub off: Vec<f64>,
    pub cells: usize,
}

impl Links {
    /// `reads` grouped by `gene_of[site]` within each cell.
    pub fn new(reads: &Weighted, gene_of: &[u32]) -> Self {
        let mut links = Links {
            cell_ptr: vec![0],
            cell: Vec::new(),
            ptr: vec![0],
            site: Vec::with_capacity(reads.on.len()),
            on: Vec::with_capacity(reads.on.len()),
            off: Vec::with_capacity(reads.on.len()),
            cells: reads.cells(),
        };
        for c in 0..reads.cells() {
            let (s, m, u) = reads.of(c);
            let mut order: Vec<usize> = (0..s.len()).collect();
            order.sort_by_key(|&i| (gene_of[s[i] as usize], s[i]));
            for (j, &i) in order.iter().enumerate() {
                let gene = gene_of[s[i] as usize];
                if j > 0 && gene_of[s[order[j - 1]] as usize] != gene {
                    links.cell.push(c as u32);
                    links.ptr.push(links.site.len());
                }
                links.site.push(s[i]);
                links.on.push(m[i]);
                links.off.push(u[i]);
            }
            if !order.is_empty() {
                links.cell.push(c as u32);
                links.ptr.push(links.site.len());
            }
            links.cell_ptr.push(links.cell.len());
        }
        links
    }

    /// Each link's gene, by `gene_of[site]` of its first site.
    pub fn genes_of(&self, gene_of: &[u32]) -> Vec<u32> {
        (0..self.len())
            .map(|l| gene_of[self.site[self.ptr[l]] as usize])
            .collect()
    }

    pub fn len(&self) -> usize {
        self.cell.len()
    }

    /// Link `l`'s `(sites, converted, unconverted)`.
    pub(super) fn of(&self, l: usize) -> (&[u32], &[f64], &[f64]) {
        let r = self.ptr[l]..self.ptr[l + 1];
        (&self.site[r.clone()], &self.on[r.clone()], &self.off[r])
    }
}

#[derive(Debug, Clone)]
pub struct Settings {
    /// Communities.
    pub k: usize,
    /// Dirichlet concentration of each cell's membership around the global
    /// weights.
    pub alpha: f64,
    /// Concentration of the global weights (`Dirichlet(gamma / K)`).
    pub gamma: f64,
    /// Prior strength `A` of each site's fraction in a community.
    pub strength: f64,
}

pub struct Model<'a> {
    links: &'a Links,
    sites: usize,
    a: Vec<f64>,
    b: Vec<f64>,
    ln_b0: Vec<f64>,
    settings: Settings,
    /// Each cell's held-out reads, scored at every kept draw.
    heldout: Option<&'a Reads>,
}

/// A chain's link labels and the counts they imply: links per cell and
/// community (`cells x K`), reads per site and community (`sites x K`).
pub struct State {
    z: Vec<u16>,
    n: Vec<u32>,
    on: Vec<f64>,
    off: Vec<f64>,
    /// Global community weights.
    beta: Vec<f64>,
    /// The chain's tally, shared with its draws.
    tally: Arc<Mutex<Tally>>,
}

/// What a chain keeps across its kept draws, updated in place at each:
/// its latest labels and the running sum of each cell's membership. Every
/// draw of the chain holds the same tally, so a chain costs one copy of the
/// labels, not one per draw.
#[derive(Debug, Default)]
pub struct Tally {
    z: Vec<u16>,
    membership: Vec<f64>,
    draws: usize,
}

/// A chain's labels at its last kept draw, and each cell's posterior-mean
/// membership (`cells x K`, row-major) over its kept draws.
pub struct ChainMeans {
    pub z: Vec<u16>,
    pub membership: Vec<f32>,
}

/// One kept draw: its scalars, and its chain's shared tally.
#[derive(Debug, Clone)]
pub struct Draw {
    /// Held-out reads' log-likelihood and read count, at this draw.
    pub heldout: (f64, f64),
    pub beta: Vec<f64>,
    pub log_joint: f64,
    tally: Arc<Mutex<Tally>>,
    /// The labels at this draw, kept only in tests.
    #[cfg(test)]
    pub z: Vec<u16>,
}

impl Draw {
    /// The chain's means and latest labels, as of its last kept draw.
    pub fn means(&self) -> ChainMeans {
        let t = self.tally.lock().expect("a chain's tally");
        let n = t.draws.max(1) as f64;
        ChainMeans {
            z: t.z.clone(),
            membership: t.membership.iter().map(|x| (x / n) as f32).collect(),
        }
    }
}

impl<'a> Model<'a> {
    pub fn new(links: &'a Links, sites: usize, settings: Settings) -> Self {
        let p = pooled_fractions(&links.site, &links.on, &links.off, sites);
        let a: Vec<f64> = p.iter().map(|p| settings.strength * p).collect();
        let b: Vec<f64> = p.iter().map(|p| settings.strength * (1.0 - p)).collect();
        let ln_b0 = a.iter().zip(&b).map(|(&x, &y)| ln_beta(x, y)).collect();
        Self {
            links,
            sites,
            a,
            b,
            ln_b0,
            settings,
            heldout: None,
        }
    }

    /// Score `heldout` (each cell's held-out reads) at every kept draw.
    pub fn with_heldout(mut self, heldout: &'a Reads) -> Self {
        self.heldout = Some(heldout);
        self
    }

    /// Held-out reads' log-likelihood and count, each read at its cell's
    /// posterior-mean fraction: each community's share of the cell times
    /// the community's posterior-mean site fraction.
    fn score_heldout(&self, st: &State) -> (f64, f64) {
        let Some(held) = self.heldout else {
            return (0.0, 0.0);
        };
        let kk = self.settings.k;
        let mut weight = vec![0.0f64; held.cells() * kk];
        for (c, w) in weight.chunks_mut(kk).enumerate() {
            self.label_prior(&st.n[c * kk..(c + 1) * kk], &st.beta, w);
        }
        // Converted and unconverted each summed on its own, so neither is
        // taken as 1 less the other.
        super::read::score(held, |c, j| {
            let w = &weight[c * kk..(c + 1) * kk];
            let p: f64 = (0..kk).map(|k| w[k] * self.fraction(st, j, k)).sum();
            let q: f64 = (0..kk)
                .map(|k| w[k] * (1.0 - self.fraction(st, j, k)))
                .sum();
            (p.ln(), q.ln())
        })
    }

    /// A new link's community probabilities in a cell with links per
    /// community `row`: its share of the cell's links around the global
    /// weights `beta`.
    fn label_prior(&self, row: &[u32], beta: &[f64], out: &mut [f64]) {
        let alpha = self.settings.alpha;
        let total = row.iter().sum::<u32>() as f64 + alpha;
        for (k, w) in out.iter_mut().enumerate() {
            *w = (row[k] as f64 + alpha * beta[k]) / total;
        }
    }

    /// Site `s`'s posterior-mean fraction in community `k`.
    fn fraction(&self, st: &State, s: usize, k: usize) -> f64 {
        let at = s * self.settings.k + k;
        (self.a[s] + st.on[at]) / (self.a[s] + self.b[s] + st.on[at] + st.off[at])
    }

    fn put(&self, st: &mut State, l: usize, k: usize, sign: f64) {
        let kk = self.settings.k;
        let c = self.links.cell[l] as usize;
        st.n[c * kk + k] = (st.n[c * kk + k] as i64 + sign as i64) as u32;
        let (s, m, u) = self.links.of(l);
        for i in 0..s.len() {
            let at = s[i] as usize * kk + k;
            st.on[at] = (st.on[at] + sign * m[i]).max(0.0);
            st.off[at] = (st.off[at] + sign * u[i]).max(0.0);
        }
    }

    /// Log weight of link `l` in each community, given the others: its
    /// cell's community prior times its reads' Beta-binomial predictive.
    fn weights(&self, st: &State, l: usize, out: &mut [f64]) {
        let kk = self.settings.k;
        let c = self.links.cell[l] as usize;
        self.label_prior(&st.n[c * kk..(c + 1) * kk], &st.beta, out);
        let (s, m, u) = self.links.of(l);
        for (k, w) in out.iter_mut().enumerate() {
            *w = w.ln();
            for i in 0..s.len() {
                let (j, at) = (s[i] as usize, s[i] as usize * kk + k);
                *w += ln_pred(self.a[j] + st.on[at], self.b[j] + st.off[at], m[i], u[i]);
            }
        }
    }

    /// A cell swap: for cell `c` and two communities `k`, `k2` drawn
    /// uniformly, propose exchanging the two labels on all of `c`'s links.
    /// A swap is its own inverse, so the proposal is symmetric and the
    /// Metropolis-Hastings ratio is the change in log joint: the cell's own
    /// term and the two communities' sites it touches. Lets a whole cell
    /// move to another community, which single-link moves cannot do against
    /// the cell's own urn.
    fn swap_cell(&self, st: &mut State, c: usize, rng: &mut SmallRng) {
        let kk = self.settings.k;
        if kk < 2 {
            return;
        }
        let k = rng.random_range(0..kk);
        let k2 = (k + rng.random_range(1..kk)) % kk;
        let Some((delta, links)) = self.swap_delta(st, c, k, k2) else {
            return;
        };
        if rng.random::<f64>().ln() < delta {
            for &l in &links {
                self.put(st, l, st.z[l] as usize, -1.0);
            }
            for &l in &links {
                let to = if st.z[l] as usize == k { k2 } else { k };
                st.z[l] = to as u16;
                self.put(st, l, to, 1.0);
            }
        }
    }

    /// The change in log joint of exchanging communities `k` and `k2` on
    /// cell `c`'s links, and the links it moves (none: nothing to swap).
    fn swap_delta(&self, st: &State, c: usize, k: usize, k2: usize) -> Option<(f64, Vec<usize>)> {
        let kk = self.settings.k;
        let links: Vec<usize> = (self.links.cell_ptr[c]..self.links.cell_ptr[c + 1])
            .filter(|&l| st.z[l] as usize == k || st.z[l] as usize == k2)
            .collect();
        if links.is_empty() {
            return None;
        }
        let row = &st.n[c * kk..(c + 1) * kk];
        let mut swapped = row.to_vec();
        swapped.swap(k, k2);
        let mut delta = self.cell_term(&swapped, &st.beta) - self.cell_term(row, &st.beta);
        for &l in &links {
            // A cell covers a site once, so each site here moves one way.
            let (from, to) = if st.z[l] as usize == k {
                (k, k2)
            } else {
                (k2, k)
            };
            // Each site's change in its two communities' marginals: the
            // reads' predictive in `to` less their predictive in `from`
            // (both without the reads themselves).
            let (s, m, u) = self.links.of(l);
            for i in 0..s.len() {
                let j = s[i] as usize;
                let (a, b) = (j * kk + from, j * kk + to);
                let (x, y) = (self.a[j], self.b[j]);
                delta += ln_pred(x + st.on[b], y + st.off[b], m[i], u[i])
                    - ln_pred(x + st.on[a] - m[i], y + st.off[a] - u[i], m[i], u[i]);
            }
        }
        Some((delta, links))
    }

    /// For labels `z` at global weights `beta`: the swap's computed change
    /// in log joint for cell `c` and communities `k`, `k2`, and the change
    /// the log joint itself shows (both none when nothing would move).
    #[cfg(test)]
    pub fn swap_check(
        &self,
        z: &[u16],
        beta: &[f64],
        c: usize,
        k: usize,
        k2: usize,
    ) -> Option<(f64, f64)> {
        let mut st = self.state(z.to_vec());
        st.beta = beta.to_vec();
        let (delta, _) = self.swap_delta(&st, c, k, k2)?;
        let swapped: Vec<u16> = z
            .iter()
            .enumerate()
            .map(|(l, &x)| {
                let mine = self.links.cell[l] as usize == c;
                match x as usize {
                    y if mine && y == k => k2 as u16,
                    y if mine && y == k2 => k as u16,
                    _ => x,
                }
            })
            .collect();
        Some((
            delta,
            self.log_joint_of(&swapped, beta) - self.log_joint_of(z, beta),
        ))
    }

    /// A cell's term of the log joint, from its links per community `row`:
    /// its urn around the global weights (its membership integrated out).
    fn cell_term(&self, row: &[u32], beta: &[f64]) -> f64 {
        ln_urn(row, self.settings.alpha, beta)
    }

    /// Site `s`'s term of the log joint in a community holding reads
    /// `(on, off)` there: the Beta-binomial marginal.
    fn site_term(&self, s: usize, on: f64, off: f64) -> f64 {
        ln_beta(self.a[s] + on, self.b[s] + off) - self.ln_b0[s]
    }

    /// Log joint of the labels and the reads given the global weights
    /// (memberships and site fractions integrated out; the weights' prior
    /// not included).
    fn log_joint(&self, st: &State) -> f64 {
        let kk = self.settings.k;
        let cells: f64 = (0..self.links.cells)
            .map(|c| self.cell_term(&st.n[c * kk..(c + 1) * kk], &st.beta))
            .sum();
        // A site and community without reads adds exactly 0.
        let sites: f64 = (0..self.sites * kk)
            .filter(|&at| st.on[at] + st.off[at] > 0.0)
            .map(|at| self.site_term(at / kk, st.on[at], st.off[at]))
            .sum();
        cells + sites
    }

    /// The counts implied by labels `z`.
    fn state(&self, z: Vec<u16>) -> State {
        let kk = self.settings.k;
        let mut st = State {
            z: Vec::new(),
            n: vec![0; self.links.cells * kk],
            on: vec![0.0; self.sites * kk],
            off: vec![0.0; self.sites * kk],
            beta: vec![1.0 / kk as f64; kk],
            tally: Arc::default(),
        };
        for (l, &k) in z.iter().enumerate() {
            self.put(&mut st, l, k as usize, 1.0);
        }
        st.z = z;
        st
    }

    /// [`Model::log_joint`] of labels `z` at global weights `beta`.
    #[cfg(test)]
    pub fn log_joint_of(&self, z: &[u16], beta: &[f64]) -> f64 {
        let mut st = self.state(z.to_vec());
        st.beta = beta.to_vec();
        self.log_joint(&st)
    }

    /// `beta_1` after each of `draws` updates of the global weights, the
    /// labels held at `z` (to check the update against its conditional).
    #[cfg(test)]
    pub fn beta_draws(&self, z: &[u16], draws: usize, seed: u64) -> Vec<f64> {
        let mut rng = <SmallRng as rand::SeedableRng>::seed_from_u64(seed);
        let mut st = self.state(z.to_vec());
        (0..draws)
            .map(|_| {
                self.resample_beta(&mut st, &mut rng);
                st.beta[0]
            })
            .collect()
    }

    /// The global weights from their conditional: each cell's table counts
    /// in each community (an Antoniak draw), then
    /// `beta ~ Dirichlet(gamma / K + tables)`.
    fn resample_beta(&self, st: &mut State, rng: &mut SmallRng) {
        let kk = self.settings.k;
        let mut tables = vec![0.0f64; kk];
        for c in 0..self.links.cells {
            for (k, t) in tables.iter_mut().enumerate() {
                let ab = self.settings.alpha * st.beta[k];
                for j in 0..st.n[c * kk + k] {
                    if rng.random::<f64>() < ab / (ab + j as f64) {
                        *t += 1.0;
                    }
                }
            }
        }
        let draws: Vec<f64> = tables
            .iter()
            .map(|&t| {
                rand_distr::Gamma::new(self.settings.gamma / kk as f64 + t, 1.0)
                    .unwrap()
                    .sample(rng)
                    .max(1e-300)
            })
            .collect();
        let total: f64 = draws.iter().sum();
        st.beta = draws.iter().map(|d| d / total).collect();
    }

    /// Each cell's posterior-mean membership given the labels, added to
    /// `sums` (`cells x K`).
    fn add_membership(&self, st: &State, sums: &mut [f64]) {
        let kk = self.settings.k;
        let mut w = vec![0.0f64; kk];
        for (c, row) in sums.chunks_mut(kk).enumerate() {
            self.label_prior(&st.n[c * kk..(c + 1) * kk], &st.beta, &mut w);
            row.iter_mut().zip(&w).for_each(|(s, x)| *s += x);
        }
    }
}

impl McmcModel for Model<'_> {
    type State = State;
    type Sample = Draw;
    type Result = Vec<Draw>;

    /// Each cell's links all in one community, drawn per cell.
    fn init(&self, rng: &mut SmallRng) -> State {
        let per_cell: Vec<u16> = (0..self.links.cells)
            .map(|_| rng.random_range(0..self.settings.k) as u16)
            .collect();
        let z = self
            .links
            .cell
            .iter()
            .map(|&c| per_cell[c as usize])
            .collect();
        self.state(z)
    }

    /// Every link relabelled once, in random order; then a swap per cell
    /// and the global weights.
    fn sweep(&self, st: &mut State, rng: &mut SmallRng) {
        let mut order: Vec<usize> = (0..self.links.len()).collect();
        order.shuffle(rng);
        let mut w = vec![0.0f64; self.settings.k];
        for l in order {
            self.put(st, l, st.z[l] as usize, -1.0);
            self.weights(st, l, &mut w);
            let k = sample_log(&mut w, rng);
            st.z[l] = k as u16;
            self.put(st, l, k, 1.0);
        }
        for c in 0..self.links.cells {
            self.swap_cell(st, c, rng);
        }
        self.resample_beta(st, rng);
    }

    fn collect(&self, st: &State) -> Draw {
        {
            let mut t = st.tally.lock().expect("a chain's tally");
            t.z.clone_from(&st.z);
            t.membership.resize(self.links.cells * self.settings.k, 0.0);
            self.add_membership(st, &mut t.membership);
            t.draws += 1;
        }
        Draw {
            heldout: self.score_heldout(st),
            beta: st.beta.clone(),
            log_joint: self.log_joint(st),
            tally: Arc::clone(&st.tally),
            #[cfg(test)]
            z: st.z.clone(),
        }
    }

    fn summarize(&self, draws: Vec<Draw>) -> Vec<Draw> {
        draws
    }
}

/// Each community's share of the reads, under labels `z`.
pub fn read_shares(links: &Links, z: &[u16], k: usize) -> Vec<f64> {
    let mut reads = vec![0.0f64; k];
    for (l, &zl) in z.iter().enumerate() {
        let (_, m, u) = links.of(l);
        reads[zl as usize] += m.iter().sum::<f64>() + u.iter().sum::<f64>();
    }
    let total: f64 = reads.iter().sum::<f64>().max(1e-12);
    reads.iter().map(|r| r / total).collect()
}

/// Each community's reads per site and posterior-mean fraction (`sites x
/// K`), from labels `z`.
pub struct Profile {
    pub k: usize,
    pub on: Vec<f64>,
    pub total: Vec<f64>,
    pub mean: Vec<f64>,
}

impl Model<'_> {
    pub fn profile(&self, z: &[u16]) -> Profile {
        let st = self.state(z.to_vec());
        let kk = self.settings.k;
        let mean = (0..self.sites * kk)
            .map(|at| self.fraction(&st, at / kk, at % kk))
            .collect();
        Profile {
            k: kk,
            total: st.on.iter().zip(&st.off).map(|(a, b)| a + b).collect(),
            on: st.on,
            mean,
        }
    }

    /// Cells outside the fit, each placed against fitted `profile` and
    /// global weights `beta`: its links relabelled by collapsed Gibbs, the
    /// profiles held fixed, and its membership averaged over the sweeps
    /// after a burn-in. The cell's own urn is the model's; the profiles are
    /// plug-in estimates, not integrated over.
    pub fn place(&self, links: &Links, profile: &Profile, beta: &[f64], seed: u64) -> Vec<f32> {
        use rayon::prelude::*;
        let kk = self.settings.k;
        // Beta posterior means, strictly inside (0, 1).
        let ln_p: Vec<(f64, f64)> = profile
            .mean
            .iter()
            .map(|&p| (p.ln(), (1.0 - p).ln()))
            .collect();
        (0..links.cells)
            .into_par_iter()
            .flat_map_iter(|c| {
                let mut rng = cell_rng(seed, c);
                let own = links.cell_ptr[c]..links.cell_ptr[c + 1];
                // Each link's reads' log-likelihood in each community.
                let ln_p = &ln_p;
                let like: Vec<f64> = own
                    .clone()
                    .flat_map(|l| {
                        let (s, m, u) = links.of(l);
                        (0..kk).map(move |k| {
                            (0..s.len())
                                .map(|i| {
                                    let (a, b) = ln_p[s[i] as usize * kk + k];
                                    m[i] * a + u[i] * b
                                })
                                .sum::<f64>()
                        })
                    })
                    .collect();
                let mut row = vec![0u32; kk];
                let mut z: Vec<usize> = own.map(|_| rng.random_range(0..kk)).collect();
                z.iter().for_each(|&k| row[k] += 1);
                let mut membership = vec![0.0f64; kk];
                let mut w = vec![0.0f64; kk];
                for sweep in 0..PLACE_SWEEPS {
                    for (i, zi) in z.iter_mut().enumerate() {
                        row[*zi] -= 1;
                        self.label_prior(&row, beta, &mut w);
                        for (k, x) in w.iter_mut().enumerate() {
                            *x = x.ln() + like[i * kk + k];
                        }
                        *zi = sample_log(&mut w, &mut rng);
                        row[*zi] += 1;
                    }
                    if sweep >= PLACE_BURN_IN {
                        self.label_prior(&row, beta, &mut w);
                        membership.iter_mut().zip(&w).for_each(|(m, x)| *m += x);
                    }
                }
                let n = (PLACE_SWEEPS - PLACE_BURN_IN) as f64;
                membership.into_iter().map(move |x| (x / n) as f32)
            })
            .collect()
    }
}

/// Sweeps of a placed cell's own Gibbs, and those discarded before its
/// averages begin.
pub(super) const PLACE_SWEEPS: usize = 60;
pub(super) const PLACE_BURN_IN: usize = 20;

/// Placed cell `c`'s own random stream, from `seed`.
pub(super) fn cell_rng(seed: u64, c: usize) -> SmallRng {
    <SmallRng as rand::SeedableRng>::seed_from_u64(
        seed ^ (c as u64).wrapping_mul(0x9e37_79b9_7f4a_7c15),
    )
}

/// What the chains say, together.
#[derive(Debug, Clone, serde::Serialize)]
pub struct Summary {
    pub communities: usize,
    pub chains: usize,
    pub draws: usize,
    /// Split R-hat of the log joint across chains.
    pub rhat_log_joint: f64,
    /// Each chain's mean log joint.
    pub chain_log_joint: Vec<f64>,
    /// Adjusted Rand index between each pair of chains' cells' main
    /// communities (posterior-mean memberships).
    pub chain_agreement: Vec<f64>,
    /// The best chain (highest mean log joint): each community's share of
    /// the reads at its last draw, and its global weights there.
    pub best_chain: usize,
    pub read_shares: Vec<f64>,
    pub global_weights: Vec<f64>,
    /// Held-out reads' log-likelihood per read, over every kept draw.
    pub heldout_per_read: Option<f64>,
}

pub fn summarize(model: &Model, links: &Links, chains: &[Vec<Draw>]) -> Summary {
    use super::summary::{chain_scores, main_labels, pairwise_ari, per_read};
    let k = model.settings.k;
    let log_joint: Vec<Vec<f64>> = chains
        .iter()
        .map(|c| c.iter().map(|d| d.log_joint).collect())
        .collect();
    let (chain_log_joint, rhat_log_joint, best_chain) = chain_scores(&log_joint);
    let means: Vec<ChainMeans> = chains
        .iter()
        .filter_map(|c| c.last().map(Draw::means))
        .collect();
    Summary {
        communities: k,
        chains: chains.len(),
        draws: chains.iter().map(Vec::len).sum(),
        rhat_log_joint,
        chain_log_joint,
        chain_agreement: pairwise_ari(
            &means
                .iter()
                .map(|m| main_labels(&m.membership, k))
                .collect::<Vec<_>>(),
        ),
        best_chain,
        read_shares: means
            .get(best_chain)
            .map_or_else(Vec::new, |m| read_shares(links, &m.z, k)),
        global_weights: chains[best_chain]
            .last()
            .map_or_else(Vec::new, |d| d.beta.clone()),
        heldout_per_read: per_read(chains.iter().flatten().map(|d| d.heldout)),
    }
}
