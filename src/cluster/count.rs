//! A Poisson blockmodel of gene counts with a background: cell roles and
//! gene roles over a cell x gene count matrix, after Ball, Karrer and
//! Newman (2011) and the Bayesian Poisson Tucker decomposition (Schein et
//! al. 2016).
//!
//! ```text
//! x_cg ~ Poisson(s_c r_g + sum_ab theta_ca lambda_ab psi_gb)
//! theta_ca ~ Gamma(w_a, rho_c),   psi_gb ~ Gamma(v_b, rho'_g),   lambda_ab ~ Gamma(kappa, kappa)
//! s_c, r_g ~ Gamma(a0, b0),       rho_c, rho'_g ~ Gamma(e0, f0)
//! w_a ~ Gamma(gamma / A, c0),     v_b ~ Gamma(gamma / B, c0)
//! ```
//!
//! (Gamma by shape and rate.) The background `s_c r_g`, a cell's size times
//! a gene's level, belongs to no role: it takes the counts that no block
//! needs, so a block is a departure from it and no role has to hold the
//! bulk of the data (the resolution limit). Each count unit (UMI) is in the
//! background or in one block `(a, b)`; a cell's and a gene's units spread
//! over roles, so cells and genes have mixed memberships. The global
//! weights `w`, `v` are a Gamma-negative-binomial process's (Zhou and
//! Carin 2015): roles not needed empty.
//!
//! Each sweep, uncollapsed Gibbs, every draw from its exact conditional:
//! each nonzero count's units shared out over the background and the
//! blocks, multinomially by their rates (none kept between sweeps); then
//! the global weights through each unit's table counts (a CRT draw, its
//! role's weights integrated out); the cells' weights, their rates, the
//! genes' likewise, the blocks' rates and the background. A zero count
//! enters only through the exposures, `sum_cg` of each rate, which factor
//! into column totals (`sum_ab lambda_ab (sum_c theta_ca)(sum_g psi_gb)`
//! and `(sum_c s_c)(sum_g r_g)`): no zero is visited.

use rand::rngs::SmallRng;
use rand::RngExt;
use rand_distr::{Binomial, Distribution, Gamma};

/// A cell x gene count matrix, its nonzero entries by cell: cell `c`'s are
/// `ptr[c]..ptr[c + 1]` of `gene` and `x`.
#[derive(Debug, Clone, Default)]
pub struct Counts {
    pub genes: usize,
    pub ptr: Vec<usize>,
    pub gene: Vec<u32>,
    pub x: Vec<u32>,
    /// `sum ln x!` over the nonzero counts: the likelihood's constant.
    ln_factorials: f64,
}

impl Counts {
    pub fn cells(&self) -> usize {
        self.ptr.len().saturating_sub(1)
    }

    /// Cell `c`'s `(genes, counts)`.
    pub fn of(&self, c: usize) -> (&[u32], &[u32]) {
        let r = self.ptr[c]..self.ptr[c + 1];
        (&self.gene[r.clone()], &self.x[r])
    }

    /// From each cell's `(gene, count)` entries, zeros dropped.
    pub fn from_cells<'c>(genes: usize, cells: impl IntoIterator<Item = &'c [(u32, u32)]>) -> Self {
        let mut out = Counts {
            genes,
            ptr: vec![0],
            ..Counts::default()
        };
        for cell in cells {
            for &(g, x) in cell.iter().filter(|e| e.1 > 0) {
                out.gene.push(g);
                out.x.push(x);
                out.ln_factorials += statrs::function::gamma::ln_gamma(x as f64 + 1.0);
            }
            out.ptr.push(out.gene.len());
        }
        out
    }
}

/// The model's sizes and priors.
#[derive(Debug, Clone, PartialEq, serde::Serialize, serde::Deserialize)]
pub struct Settings {
    pub cell_roles: usize,
    pub gene_roles: usize,
    /// The global weights' concentration `gamma` and rate `c0`.
    pub gamma: f64,
    pub weight_rate: f64,
    /// The blocks' rates' shape and rate `kappa` (mean 1).
    pub kappa: f64,
    /// The background's sizes' and levels' shape and rate.
    pub background: (f64, f64),
    /// The cells' and genes' rates' shape and rate.
    pub rate_prior: (f64, f64),
}

impl Default for Settings {
    fn default() -> Self {
        Self {
            cell_roles: 8,
            gene_roles: 8,
            gamma: 1.0,
            weight_rate: 1.0,
            kappa: 1.0,
            background: (0.1, 0.1),
            rate_prior: (0.1, 0.1),
        }
    }
}

/// The model over one count matrix.
pub struct Model<'a> {
    counts: &'a Counts,
    settings: Settings,
}

/// A chain's parameters, and the units' allocation of its last sweep.
#[derive(Debug, Clone)]
pub struct State {
    /// `cells x A`, `genes x B`, `A x B`, row-major.
    pub theta: Vec<f64>,
    pub psi: Vec<f64>,
    pub lambda: Vec<f64>,
    /// Each cell's size and each gene's level in the background.
    pub size: Vec<f64>,
    pub level: Vec<f64>,
    /// Each cell's and each gene's rate in its weights' prior.
    pub cell_rate: Vec<f64>,
    pub gene_rate: Vec<f64>,
    /// The global weights over the cell roles and the gene roles.
    pub w: Vec<f64>,
    pub v: Vec<f64>,
    /// Units per cell and role, gene and role, block; and in the
    /// background per cell and per gene.
    pub n: Vec<u32>,
    pub m: Vec<u32>,
    pub blocks: Vec<u32>,
    pub n0: Vec<u32>,
    pub m0: Vec<u32>,
}

/// One draw of a Gamma by shape and rate, kept off 0.
fn gamma(shape: f64, rate: f64, rng: &mut SmallRng) -> f64 {
    Gamma::new(shape, 1.0 / rate)
        .expect("a positive shape and rate")
        .sample(rng)
        .max(1e-300)
}

/// The tables of `n` customers at concentration `r` (a Chinese
/// restaurant table count): the first opens one for certain, the `j`-th
/// after it one with probability `r / (r + j)`.
fn tables(n: u32, r: f64, rng: &mut SmallRng) -> u32 {
    if n == 0 {
        return 0;
    }
    1 + (1..n)
        .filter(|&j| rng.random::<f64>() * (r + j as f64) < r)
        .count() as u32
}

/// Units at or under which [`multinomial`] draws each unit's bucket
/// (one uniform each) rather than a binomial per bucket: most counts are
/// a few, and a binomial costs a power or more to set up.
const ONE_BY_ONE: u32 = 12;

/// `n` units shared out over `weights` (all finite and at least 0, not all
/// 0), multinomially, into `out` (added to): a few units one by one, more
/// as a binomial per bucket.
fn multinomial(n: u32, weights: &[f64], out: &mut [u32], rng: &mut SmallRng) {
    if n <= ONE_BY_ONE {
        let total: f64 = weights.iter().sum();
        for _ in 0..n {
            let mut u = rng.random::<f64>() * total;
            let mut k = weights.len() - 1;
            for (i, &w) in weights.iter().enumerate() {
                if u < w {
                    k = i;
                    break;
                }
                u -= w;
            }
            out[k] += 1;
        }
        return;
    }
    let mut left = n;
    let mut rest: f64 = weights.iter().sum();
    // The last bucket of positive weight takes what is left, so no unit
    // is lost to rounding in `rest`.
    let last = weights.iter().rposition(|&w| w > 0.0).unwrap_or(0);
    for (i, (o, &w)) in out.iter_mut().zip(weights).enumerate() {
        if left == 0 {
            break;
        }
        let k = if i >= last {
            left
        } else {
            let p = (w / rest).clamp(0.0, 1.0);
            if p >= 1.0 {
                left
            } else {
                Binomial::new(left as u64, p)
                    .expect("a probability")
                    .sample(rng) as u32
            }
        };
        *o += k;
        left -= k;
        rest -= w;
    }
}

impl<'a> Model<'a> {
    pub fn new(counts: &'a Counts, settings: Settings) -> anyhow::Result<Self> {
        anyhow::ensure!(
            settings.cell_roles >= 1 && settings.gene_roles >= 1,
            "at least one cell role and one gene role"
        );
        Ok(Self { counts, settings })
    }

    fn sizes(&self) -> (usize, usize, usize, usize) {
        let s = &self.settings;
        (
            self.counts.cells(),
            self.counts.genes,
            s.cell_roles,
            s.gene_roles,
        )
    }

    /// A state of these sizes, every parameter and count 0.
    fn zeroed(&self) -> State {
        let (cc, gg, aa, bb) = self.sizes();
        State {
            theta: vec![0.0; cc * aa],
            psi: vec![0.0; gg * bb],
            lambda: vec![0.0; aa * bb],
            size: vec![0.0; cc],
            level: vec![0.0; gg],
            cell_rate: vec![0.0; cc],
            gene_rate: vec![0.0; gg],
            w: vec![0.0; aa],
            v: vec![0.0; bb],
            n: vec![0; cc * aa],
            m: vec![0; gg * bb],
            blocks: vec![0; aa * bb],
            n0: vec![0; cc],
            m0: vec![0; gg],
        }
    }

    /// Each gene's `sum_b lambda_ab psi_gb` (`genes x A`): what a count of
    /// the gene in cell role `a` is rated by, per unit of theta.
    fn gene_rate_sums(&self, st: &State) -> Vec<f64> {
        let (_, gg, aa, bb) = self.sizes();
        let mut sums = vec![0.0; gg * aa];
        for (g, row) in sums.chunks_exact_mut(aa).enumerate() {
            let psi = &st.psi[g * bb..(g + 1) * bb];
            for (a, x) in row.iter_mut().enumerate() {
                let lambda = &st.lambda[a * bb..(a + 1) * bb];
                *x = lambda.iter().zip(psi).map(|(l, p)| l * p).sum();
            }
        }
        sums
    }

    /// Each gene's `lambda_ab psi_gb` (`genes x A x B`) and their sums over
    /// `b` (`genes x A`): what a count of the gene in cell role `a` and
    /// block `(a, b)` is rated by, per unit of theta.
    fn gene_rates(&self, st: &State) -> (Vec<f64>, Vec<f64>) {
        let (_, gg, aa, bb) = self.sizes();
        let mut within = vec![0.0; gg * aa * bb];
        let mut sums = vec![0.0; gg * aa];
        for g in 0..gg {
            let psi = &st.psi[g * bb..(g + 1) * bb];
            for a in 0..aa {
                let row = &mut within[(g * aa + a) * bb..(g * aa + a + 1) * bb];
                for ((x, &l), &p) in row
                    .iter_mut()
                    .zip(&st.lambda[a * bb..(a + 1) * bb])
                    .zip(psi)
                {
                    *x = l * p;
                }
                sums[g * aa + a] = row.iter().sum();
            }
        }
        (within, sums)
    }

    /// One sweep: the allocation, then the parameters given it.
    pub fn sweep(&self, st: &mut State, rng: &mut SmallRng) {
        self.allocate(st, rng);
        self.draw_parameters(st, rng);
    }

    /// Each nonzero count's units shared out over the background and the
    /// blocks by their rates, the counts per cell and role, gene and role,
    /// block and background tallied afresh.
    pub fn allocate(&self, st: &mut State, rng: &mut SmallRng) {
        let (cc, _, aa, bb) = self.sizes();
        for t in [&mut st.n, &mut st.m, &mut st.blocks, &mut st.n0, &mut st.m0] {
            t.fill(0);
        }
        let (within, sums) = self.gene_rates(st);
        // Per link: the background's rate then each cell role's (its gene
        // roles summed), and within a role each gene role's.
        let mut role = vec![0.0; 1 + aa];
        let mut units = vec![0u32; 1 + aa];
        let mut per_b = vec![0u32; bb];
        for c in 0..cc {
            let (genes, xs) = self.counts.of(c);
            let theta = &st.theta[c * aa..(c + 1) * aa];
            for (&g, &x) in genes.iter().zip(xs) {
                let g = g as usize;
                role[0] = st.size[c] * st.level[g];
                for (r, (&t, &l)) in role[1..]
                    .iter_mut()
                    .zip(theta.iter().zip(&sums[g * aa..(g + 1) * aa]))
                {
                    *r = t * l;
                }
                units.fill(0);
                multinomial(x, &role, &mut units, rng);
                st.n0[c] += units[0];
                st.m0[g] += units[0];
                for a in 0..aa {
                    let k = units[1 + a];
                    if k == 0 {
                        continue;
                    }
                    st.n[c * aa + a] += k;
                    per_b.fill(0);
                    let rates = &within[(g * aa + a) * bb..(g * aa + a + 1) * bb];
                    multinomial(k, rates, &mut per_b, rng);
                    for (b, &kb) in per_b.iter().enumerate() {
                        st.m[g * bb + b] += kb;
                        st.blocks[a * bb + b] += kb;
                    }
                }
            }
        }
    }

    /// Everything but the allocation from its conditional, the allocation
    /// held: the cell side, the gene side, the blocks, the background.
    pub fn draw_parameters(&self, st: &mut State, rng: &mut SmallRng) {
        let (cc, gg, aa, bb) = self.sizes();
        let s = &self.settings;
        let (e0, f0) = s.rate_prior;
        // The cell side: each cell role's exposure per unit of theta.
        let psi_total = column_totals(&st.psi, bb);
        let cell_exposure: Vec<f64> = (0..aa)
            .map(|a| (0..bb).map(|b| st.lambda[a * bb + b] * psi_total[b]).sum())
            .collect();
        side(
            (&st.n, &mut st.theta, &mut st.cell_rate, &mut st.w),
            &cell_exposure,
            (cc, aa),
            (s.gamma, s.weight_rate, e0, f0),
            rng,
        );
        let theta_total = column_totals(&st.theta, aa);
        let gene_exposure: Vec<f64> = (0..bb)
            .map(|b| {
                (0..aa)
                    .map(|a| st.lambda[a * bb + b] * theta_total[a])
                    .sum()
            })
            .collect();
        side(
            (&st.m, &mut st.psi, &mut st.gene_rate, &mut st.v),
            &gene_exposure,
            (gg, bb),
            (s.gamma, s.weight_rate, e0, f0),
            rng,
        );
        let psi_total = column_totals(&st.psi, bb);
        for (a, &t) in theta_total.iter().enumerate() {
            for (b, &p) in psi_total.iter().enumerate() {
                let k = a * bb + b;
                st.lambda[k] = gamma(s.kappa + st.blocks[k] as f64, s.kappa + t * p, rng);
            }
        }
        let (a0, b0) = s.background;
        let levels: f64 = st.level.iter().sum();
        for c in 0..cc {
            st.size[c] = gamma(a0 + st.n0[c] as f64, b0 + levels, rng);
        }
        let sizes: f64 = st.size.iter().sum();
        for g in 0..gg {
            st.level[g] = gamma(a0 + st.m0[g] as f64, b0 + sizes, rng);
        }
    }

    /// Each count's log-likelihood under `st`, summed: every nonzero
    /// count's `x ln mu - ln x!` less every rate's exposure.
    pub fn log_likelihood(&self, st: &State) -> f64 {
        let (cc, _, aa, bb) = self.sizes();
        let sums = self.gene_rate_sums(st);
        let mut out = -self.counts.ln_factorials;
        for c in 0..cc {
            let (genes, xs) = self.counts.of(c);
            let theta = &st.theta[c * aa..(c + 1) * aa];
            for (&g, &x) in genes.iter().zip(xs) {
                let g = g as usize;
                let rates = &sums[g * aa..(g + 1) * aa];
                let mu = st.size[c] * st.level[g]
                    + theta.iter().zip(rates).map(|(t, l)| t * l).sum::<f64>();
                out += x as f64 * mu.ln();
            }
        }
        let (theta_total, psi_total) = (column_totals(&st.theta, aa), column_totals(&st.psi, bb));
        for (a, &t) in theta_total.iter().enumerate() {
            for (b, &p) in psi_total.iter().enumerate() {
                out -= st.lambda[a * bb + b] * t * p;
            }
        }
        out - st.size.iter().sum::<f64>() * st.level.iter().sum::<f64>()
    }
}

/// Each column's total of `x`, `k` columns per row.
fn column_totals(x: &[f64], k: usize) -> Vec<f64> {
    let mut out = vec![0.0; k];
    for row in x.chunks_exact(k) {
        out.iter_mut().zip(row).for_each(|(o, v)| *o += v);
    }
    out
}

/// One side's draws, the cells' or the genes' (`units x roles`): the
/// global weights through table counts, each unit's weights integrated out
/// (a role's units in a unit are negative binomial, `NB(w_a, E_a / (rho +
/// E_a))`), then each unit's weights, then its rate.
fn side(
    (counts, weights, rate, global): (&[u32], &mut [f64], &mut [f64], &mut [f64]),
    exposure: &[f64],
    (units, roles): (usize, usize),
    (conc, conc_rate, e0, f0): (f64, f64, f64, f64),
    rng: &mut SmallRng,
) {
    for (a, w) in global.iter_mut().enumerate() {
        let (mut t, mut ln_q) = (0u64, 0.0);
        for u in 0..units {
            t += tables(counts[u * roles + a], *w, rng) as u64;
            ln_q += (1.0 + exposure[a] / rate[u]).ln();
        }
        *w = gamma(conc / roles as f64 + t as f64, conc_rate + ln_q, rng);
    }
    for u in 0..units {
        let row = &mut weights[u * roles..(u + 1) * roles];
        for (a, x) in row.iter_mut().enumerate() {
            *x = gamma(
                global[a] + counts[u * roles + a] as f64,
                rate[u] + exposure[a],
                rng,
            );
        }
        let shape: f64 = global.iter().sum();
        rate[u] = gamma(e0 + shape, f0 + row.iter().sum::<f64>(), rng);
    }
}

impl Model<'_> {
    /// A start from the data: the background at half of each cell's and
    /// gene's total (`s_c r_g` half the counts' product over their sum),
    /// each cell's and gene's other half spread over its roles at random
    /// (each role's share a unit Gamma's), the blocks' rates and the global
    /// weights at 1, the rates at their priors' means.
    pub fn start(&self, rng: &mut SmallRng) -> State {
        let (cc, gg, aa, bb) = self.sizes();
        let mut st = self.zeroed();
        let mut cell_total = vec![0.0f64; cc];
        let mut gene_total = vec![0.0f64; gg];
        for (c, n) in cell_total.iter_mut().enumerate() {
            let (genes, xs) = self.counts.of(c);
            for (&g, &x) in genes.iter().zip(xs) {
                *n += x as f64;
                gene_total[g as usize] += x as f64;
            }
        }
        let total: f64 = cell_total.iter().sum::<f64>().max(1.0);
        // s_c r_g = half of N_c N_g / N, split evenly between the two.
        let root = (2.0 * total).sqrt();
        for (s, &n) in st.size.iter_mut().zip(&cell_total) {
            *s = (n / root).max(1e-3);
        }
        for (r, &n) in st.level.iter_mut().zip(&gene_total) {
            *r = (n / root).max(1e-3);
        }
        // Blocks at rate 1: theta_ca psi_gb summed is the other half when
        // theta sums to N_c / root over a cell's roles and psi to N_g /
        // root over a gene's.
        let spread = |out: &mut [f64], totals: &[f64], k: usize, rng: &mut SmallRng| {
            for (row, &t) in out.chunks_mut(k).zip(totals) {
                let shares: Vec<f64> = (0..k).map(|_| gamma(1.0, 1.0, rng)).collect();
                let sum: f64 = shares.iter().sum();
                for (x, w) in row.iter_mut().zip(&shares) {
                    *x = (t / root * w / sum).max(1e-6);
                }
            }
        };
        spread(&mut st.theta, &cell_total, aa, rng);
        spread(&mut st.psi, &gene_total, bb, rng);
        st.lambda.fill(1.0);
        st.w.fill(1.0);
        st.v.fill(1.0);
        let (e0, f0) = self.settings.rate_prior;
        st.cell_rate.fill(e0 / f0);
        st.gene_rate.fill(e0 / f0);
        st
    }
}

#[cfg(test)]
mod tests {
    use super::multinomial;
    use rand::rngs::SmallRng;
    use rand::SeedableRng;

    /// Every unit lands in a bucket of positive weight, few or many, with
    /// weights that leave rounding in their running sum and trailing
    /// buckets of weight 0.
    #[test]
    fn a_multinomial_places_every_unit() {
        let mut rng = SmallRng::seed_from_u64(1);
        let weights = [0.1, 1e-17, 0.3, 0.2, 1e-300, 0.0, 0.0];
        for n in [1, 5, 12, 13, 40, 1000, 100_000] {
            let mut out = [0u32; 7];
            multinomial(n, &weights, &mut out, &mut rng);
            assert_eq!(out.iter().sum::<u32>(), n, "{n} units: {out:?}");
            assert_eq!(out[5] + out[6], 0, "none in a bucket of weight 0");
        }
    }
}
