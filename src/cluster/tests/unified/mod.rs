//! The blockmodel: every switch's moves against posteriors computed from
//! the model's definition, its swaps against its log joint, its presets
//! against the link and blocks models they replaced, and what the presets
//! recover, place and write. The fixtures and references are here; the
//! exact-posterior checks, the conditionals and the presets' tests in the
//! submodules.

use super::{config, draw_cells};
use crate::cluster::crp::ln_beta;
use crate::cluster::fast;
use crate::cluster::read::{self, Reads, Weighted};
use crate::cluster::summary::{self, ari};
use crate::cluster::unified::{self, Held, Links, Model, Settings, Weights};
use legume_numeric::mcmc::engine::run_mcmc_parallel;
use rand::rngs::SmallRng;
use rand::{RngExt, SeedableRng};
use rustc_hash::FxHashMap;
use statrs::function::gamma::ln_gamma;

mod conditionals;
mod exact;
mod golden;
mod pool;
mod presets;

fn sigmoid(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

/// A link labelling under the blockmodel: cell roles and gene roles.
type Labels = (Vec<u16>, Vec<u16>);

/// The tiny link fixture: two cells, two genes (sites 0-1 and site 2), four
/// links.
fn tiny_links() -> Links {
    let reads = Reads::from_cells(vec![
        vec![(0, 3, 0), (1, 0, 2), (2, 1, 1)],
        vec![(0, 0, 3), (1, 2, 0), (2, 2, 0)],
    ]);
    Links::new(&Weighted::new(&reads, None), &TINY_GENES)
}

/// Each tiny-fixture site's gene.
const TINY_GENES: [u32; 3] = [0, 0, 1];

/// The tiny link fixture with its reads fractionally weighted (each
/// entry's converted and unconverted reads their own weight), and the
/// fixture as read: the two the exact checks run over.
fn tiny_fixtures() -> [Links; 2] {
    let reads = Reads::from_cells(vec![
        vec![(0, 3, 0), (1, 0, 2), (2, 1, 1)],
        vec![(0, 0, 3), (1, 2, 0), (2, 2, 0)],
    ]);
    let (on, off) = (
        [0.6, 1.0, 1.4, 1.0, 0.8, 1.0],
        [1.0, 1.7, 0.5, 0.9, 1.0, 1.2],
    );
    [
        tiny_links(),
        Links::new(&Weighted::new(&reads, Some((&on, &off))), &TINY_GENES),
    ]
}

/// Today's link model on the tiny fixture (two communities, alpha 0.7,
/// gamma 2, prior strength 2), enumerated by that model before the
/// blockmodel replaced it: each labelling of the four links (link `l`'s
/// community bit `l` of the index) with its posterior probability, and the
/// posterior mean of the global weight of link 0's community. A fixed
/// reference, so the link preset is held to the model it replaced.
const LINK_REFERENCE: ([f64; 16], f64) = (
    [
        0.021108756186692907,
        0.06527071045966835,
        0.003071562845161041,
        0.3241003416579008,
        0.06527071045966835,
        0.00044126670890401844,
        0.017665088836843514,
        0.0030715628451610352,
        0.003071562845161041,
        0.01766508883684353,
        0.00044126670890401844,
        0.06527071045966823,
        0.3241003416579011,
        0.0030715628451610352,
        0.06527071045966823,
        0.021108756186692907,
    ],
    0.513412103495515,
);

/// Every labelling of `links` under the blocks model the blockmodel
/// replaced (two cell roles, two gene roles), with its exact posterior
/// probability, the block shifts held at `shifts` and the cells' levels at
/// `levels`: each cell's and gene's Dirichlet-multinomial urn around equal
/// weights (`alpha`, `gene_alpha`), and Binomial reads at `logit p = site
/// baseline mu + cell level + shift`. Written from that model's
/// definition, a fixed reference for the blocks preset.
fn exact_block_labels(
    links: &Links,
    gene_of: &[u32],
    mu: &[f64],
    (alpha, gene_alpha): (f64, f64),
    (shifts, levels): (&[f64], &[f64]),
) -> Vec<(Labels, f64)> {
    let link_gene = links.genes_of(gene_of);
    let n = links.len();
    // Each labelling, (cell roles, gene roles), with its log probability.
    let mut exact: Vec<(Labels, f64)> = Vec::new();
    for code in 0..1usize << (2 * n) {
        let a: Vec<u16> = (0..n).map(|l| ((code >> (2 * l)) & 1) as u16).collect();
        let b: Vec<u16> = (0..n).map(|l| ((code >> (2 * l + 1)) & 1) as u16).collect();
        let mut v = 0.0;
        for c in 0..2 {
            let mut counts = [0.0; 2];
            (0..n)
                .filter(|&l| links.cell[l] as usize == c)
                .for_each(|l| counts[a[l] as usize] += 1.0);
            v += urn(&counts, alpha, &[0.5; 2]);
        }
        for g in 0..2 {
            let mut counts = [0.0; 2];
            (0..n)
                .filter(|&l| link_gene[l] as usize == g)
                .for_each(|l| counts[b[l] as usize] += 1.0);
            v += urn(&counts, gene_alpha, &[0.5; 2]);
        }
        for l in 0..n {
            let shift = shifts[a[l] as usize * 2 + b[l] as usize] + levels[links.cell[l] as usize];
            let (sites, m, u) = links.of(l);
            for i in 0..sites.len() {
                let p = sigmoid(mu[sites[i] as usize] + shift);
                v += m[i] * p.ln() + u[i] * (1.0 - p).ln();
            }
        }
        exact.push(((a, b), v));
    }
    let mx = exact.iter().map(|e| e.1).fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = exact.iter().map(|e| (e.1 - mx).exp()).sum();
    exact
        .into_iter()
        .map(|(z, l)| (z, (l - mx).exp() / total))
        .collect()
}

/// Total variation between the frequencies of labellings `draws` and
/// `exact`.
fn label_tv(draws: Vec<Labels>, exact: &[(Labels, f64)]) -> f64 {
    let total = draws.len() as f64;
    let mut freq: FxHashMap<Labels, f64> = FxHashMap::default();
    for z in draws {
        *freq.entry(z).or_default() += 1.0 / total;
    }
    exact
        .iter()
        .map(|(z, p)| (freq.get(z).copied().unwrap_or(0.0) - p).abs())
        .sum::<f64>()
        / 2.0
}

/// Cells of three planted within-gene patterns, genes of four sites each.
fn planted_links(cells: usize, genes: usize, seed: u64) -> (Reads, Vec<u32>, Vec<u16>) {
    let sites = genes * 4;
    let mut rng = SmallRng::seed_from_u64(seed);
    let base: Vec<f64> = (0..sites).map(|_| rng.random_range(0.1..0.4)).collect();
    let frac: Vec<Vec<f64>> = (0..3)
        .map(|_| {
            base.iter()
                .map(|&p| {
                    if rng.random::<f64>() < 0.4 {
                        rng.random_range(0.0..0.9)
                    } else {
                        p
                    }
                })
                .collect()
        })
        .collect();
    let truth: Vec<u16> = (0..cells).map(|c| (c % 3) as u16).collect();
    let fracs: Vec<&[f64]> = truth.iter().map(|&t| frac[t as usize].as_slice()).collect();
    let gene_of: Vec<u32> = (0..sites as u32).map(|s| s / 4).collect();
    (draw_cells(&fracs, 0.3, &mut rng), gene_of, truth)
}

/// Cells of three planted cell roles, genes of two planted gene roles,
/// each (cell role, gene role) block shifting its sites' methylation, and
/// each cell its own level within `±level`: the reads, each site's gene,
/// the cells' and genes' roles, and the cells' levels.
pub(super) fn planted_blocks(
    cells: usize,
    genes: usize,
    level: f64,
    seed: u64,
) -> (Reads, Vec<u32>, Vec<u16>, Vec<u16>, Vec<f64>) {
    let mut rng = SmallRng::seed_from_u64(seed);
    let sites = genes * 4;
    let gene_of: Vec<u32> = (0..sites as u32).map(|s| s / 4).collect();
    let mu: Vec<f64> = (0..sites).map(|_| rng.random_range(-1.5..0.5)).collect();
    let shift = [[1.5, -1.5], [-1.5, 1.5], [0.0, 0.0]];
    let cell_role: Vec<u16> = (0..cells).map(|c| (c % 3) as u16).collect();
    let gene_role: Vec<u16> = (0..genes).map(|g| (g % 2) as u16).collect();
    // Each cell's own level, whatever its role.
    let level: Vec<f64> = (0..cells)
        .map(|_| level * rng.random_range(-1.0..1.0))
        .collect();
    let fracs: Vec<Vec<f64>> = cell_role
        .iter()
        .zip(&level)
        .map(|(&r, &v)| {
            (0..sites)
                .map(|s| {
                    sigmoid(mu[s] + v + shift[r as usize][gene_role[gene_of[s] as usize] as usize])
                })
                .collect()
        })
        .collect();
    let refs: Vec<&[f64]> = fracs.iter().map(Vec::as_slice).collect();
    (
        draw_cells(&refs, 0.5, &mut rng),
        gene_of,
        cell_role,
        gene_role,
        level,
    )
}

fn correlation(x: &[f64], y: &[f64]) -> f64 {
    let n = x.len() as f64;
    let (mx, my) = (x.iter().sum::<f64>() / n, y.iter().sum::<f64>() / n);
    let cov: f64 = x.iter().zip(y).map(|(a, b)| (a - mx) * (b - my)).sum();
    let vx: f64 = x.iter().map(|a| (a - mx).powi(2)).sum();
    let vy: f64 = y.iter().map(|b| (b - my).powi(2)).sum();
    cov / (vx * vy).sqrt()
}

/// [`planted_blocks`] (each cell's level within `±1`) as `faba cluster`
/// reads them: `cells` cells, the first `fitted` of them fitted and the
/// rest controls, over `genes` genes; the data and the fitted cells' reads.
pub(super) fn planted_site_data(
    cells: usize,
    fitted: usize,
    genes: usize,
    seed: u64,
) -> (read::SiteData, Reads) {
    let (reads, gene_of, _, _, _) = planted_blocks(cells, genes, 1.0, seed);
    let sites = gene_of
        .iter()
        .enumerate()
        .map(|(s, &g)| read::Site {
            gene: format!("GENE{}", g + 1).into(),
            chr: "chr1".into(),
            pos: s as i64 * 10,
        })
        .collect();
    let cells = (0..cells)
        .map(|c| read::Cell {
            barcode: format!("B{c}").into(),
            batch: "b".into(),
            fit: c < fitted,
        })
        .collect();
    let fit_reads = reads.subset(&(0..fitted).collect::<Vec<_>>());
    let data = read::SiteData {
        modality: "m6A".into(),
        sites,
        cells,
        reads,
    };
    (data, fit_reads)
}

/// infinite strength, each cell's own level too).
pub(super) fn settings(
    weights: Weights,
    shift: bool,
    strength: f64,
    gene_roles: usize,
) -> Settings {
    Settings {
        cell_roles: 2,
        gene_roles,
        alpha: 0.7,
        gene_alpha: 1.3,
        weights,
        shift_sd: shift.then_some(1.0),
        strength,
        cell_sd: 0.8,
    }
}

/// Every allowed switch combination: HDP or symmetric weights, the shift on
/// or off, finite strength or infinite (with the cells' levels), one gene
/// role or two (two, or infinite strength, only with the shift).
fn combinations() -> Vec<Settings> {
    let mut out = Vec::new();
    for weights in [Weights::Hdp(1.0), Weights::Symmetric] {
        for shift in [true, false] {
            for strength in [2.0, f64::INFINITY] {
                for gene_roles in [1, 2] {
                    if shift || (gene_roles == 1 && strength.is_finite()) {
                        out.push(settings(weights, shift, strength, gene_roles));
                    }
                }
            }
        }
    }
    out
}

/// Shifts to hold the labels at: one per block, all 0 without the shift.
fn fixed_shifts(s: &Settings) -> Vec<f64> {
    let k = s.cell_roles * s.gene_roles;
    match s.shift_sd {
        Some(_) => [1.2, -0.4, -1.0, 0.7][..k].to_vec(),
        None => vec![0.0; k],
    }
}

/// Cell levels to hold the labels at: unequal at infinite strength, all 0
/// at finite strength, where there are none.
fn fixed_levels(s: &Settings, cells: usize) -> Vec<f64> {
    (0..cells)
        .map(|c| {
            if s.strength.is_finite() {
                0.0
            } else {
                [0.5, -0.3, 0.2][c % 3]
            }
        })
        .collect()
}

/// The cells' levels' log prior, `N(0, cell_sd^2)` each (its normaliser
/// left out), at infinite strength; none at finite strength.
fn level_prior(s: &Settings, nu: &[f64]) -> f64 {
    if s.strength.is_finite() {
        return 0.0;
    }
    nu.iter()
        .map(|v| -v * v / (2.0 * s.cell_sd * s.cell_sd))
        .sum()
}

/// The Dirichlet-multinomial log marginal of `counts` under
/// `Dirichlet(conc * base)`.
fn urn(counts: &[f64], conc: f64, base: &[f64]) -> f64 {
    ln_gamma(conc) - ln_gamma(conc + counts.iter().sum::<f64>())
        + counts
            .iter()
            .zip(base)
            .map(|(&x, &b)| ln_gamma(conc * b + x) - ln_gamma(conc * b))
            .sum::<f64>()
}

/// Each of `sites` sites: whether the model takes its Beta functions by
/// the fast `ln Gamma`, which it does where some link holds fractional
/// reads or more than 16.
fn fast_sites(links: &Links, sites: usize) -> Vec<bool> {
    let mut out = vec![false; sites];
    for (i, &s) in links.site.iter().enumerate() {
        let (m, u) = (links.on[i], links.off[i]);
        out[s as usize] |= m.fract() != 0.0 || u.fract() != 0.0 || m + u > 16.0;
    }
    out
}

/// A block's log marginal at one site holding `on` and `off` reads, centred
/// on `p`: Beta-binomial at strength `strength` (its Beta functions by the
/// fast `ln Gamma` if `fast`), Binomial at infinity.
fn site_marginal(p: f64, on: f64, off: f64, strength: f64, fast: bool) -> f64 {
    if strength.is_finite() {
        let (x, y) = (strength * p, strength * (1.0 - p));
        let beta = if fast { fast::ln_beta } else { ln_beta };
        beta(x + on, y + off) - beta(x, y)
    } else {
        on * p.ln() + off * (1.0 - p).ln()
    }
}

/// The labelling's terms that do not hang on the global weights, from the
/// model's definition: each gene's Dirichlet-multinomial urn and each
/// block's sites' marginals around `sigmoid(mu + delta)` (see
/// [`fast_sites`]); at infinite strength each link's Binomial reads at
/// `sigmoid(mu + nu_c + delta)`, its cell's level `nu` in.
fn genes_and_reads(
    links: &Links,
    gene_of: &[u32],
    mu: &[f64],
    s: &Settings,
    (a, b): (&[u16], &[u16]),
    (delta, nu): (&[f64], &[f64]),
) -> f64 {
    let (bb, n, sites) = (s.gene_roles, links.len(), mu.len());
    let kk = s.cell_roles * bb;
    let link_gene = links.genes_of(gene_of);
    let genes = *gene_of.iter().max().unwrap() as usize + 1;
    let gene_base = vec![1.0 / bb as f64; bb];
    let mut out = 0.0;
    for g in 0..genes {
        let mut counts = vec![0.0; bb];
        (0..n)
            .filter(|&l| link_gene[l] as usize == g)
            .for_each(|l| counts[b[l] as usize] += 1.0);
        out += urn(&counts, s.gene_alpha, &gene_base);
    }
    if s.strength.is_infinite() {
        for l in 0..n {
            let z = nu[links.cell[l] as usize] + delta[a[l] as usize * bb + b[l] as usize];
            let (ss, m, u) = links.of(l);
            for i in 0..ss.len() {
                let p = sigmoid(mu[ss[i] as usize] + z);
                out += site_marginal(p, m[i], u[i], s.strength, false);
            }
        }
        return out;
    }
    let (on, off) = block_reads(links, s, (a, b), sites);
    let fast = fast_sites(links, sites);
    for at in 0..sites * kk {
        if on[at] + off[at] > 0.0 {
            let p = sigmoid(mu[at / kk] + delta[at % kk]);
            out += site_marginal(p, on[at], off[at], s.strength, fast[at / kk]);
        }
    }
    out
}

/// Reads per site and block (`sites x K`, block `a B + b`) under labels
/// `a`, `b`: converted, unconverted.
fn block_reads(
    links: &Links,
    s: &Settings,
    (a, b): (&[u16], &[u16]),
    sites: usize,
) -> (Vec<f64>, Vec<f64>) {
    let (bb, kk) = (s.gene_roles, s.cell_roles * s.gene_roles);
    let (mut on, mut off) = (vec![0.0; sites * kk], vec![0.0; sites * kk]);
    for l in 0..links.len() {
        let k = a[l] as usize * bb + b[l] as usize;
        let (ss, m, u) = links.of(l);
        for i in 0..ss.len() {
            on[ss[i] as usize * kk + k] += m[i];
            off[ss[i] as usize * kk + k] += u[i];
        }
    }
    (on, off)
}

/// Each cell's links per cell role under labels `a`.
fn cell_counts(links: &Links, roles: usize, a: &[u16]) -> Vec<Vec<f64>> {
    let mut out = vec![vec![0.0; roles]; links.cells];
    for (l, &r) in a.iter().enumerate() {
        out[links.cell[l] as usize][r as usize] += 1.0;
    }
    out
}

/// Every labelling of two-cell-role `links` under `s` with its exact
/// posterior probability, the shifts held at `delta` and the cells' levels
/// at `nu`, from the model's
/// definition: each cell's and gene's Dirichlet-multinomial urn, each
/// block's sites' marginals around `sigmoid(mu + delta)`, and under HDP the
/// global weight `beta_1 ~ Beta(gamma / 2, gamma / 2)` integrated on a
/// grid. Also `E[beta]` of link 0's cell role.
fn exact(
    links: &Links,
    gene_of: &[u32],
    mu: &[f64],
    s: &Settings,
    held: (&[f64], &[f64]),
) -> (Vec<(Labels, f64)>, f64) {
    assert_eq!(s.cell_roles, 2);
    let (bb, kk, n) = (s.gene_roles, 2 * s.gene_roles, links.len());
    // The global weight's grid and log prior there (one point, weight 1,
    // when symmetric).
    let (grid, ln_prior): (Vec<f64>, Box<dyn Fn(f64) -> f64>) = match s.weights {
        Weights::Hdp(gamma) => {
            let h = gamma / 2.0;
            (
                (0..2000).map(|i| (i as f64 + 0.5) / 2000.0).collect(),
                Box::new(move |b: f64| (h - 1.0) * (b.ln() + (1.0 - b).ln()) - ln_beta(h, h)),
            )
        }
        Weights::Symmetric => (vec![0.5], Box::new(|_| 0.0)),
    };
    let mut out = Vec::new();
    let mut weight_mean = Vec::new();
    for code in 0..kk.pow(n as u32) {
        let k: Vec<usize> = (0..n).map(|l| code / kk.pow(l as u32) % kk).collect();
        let a: Vec<u16> = k.iter().map(|&k| (k / bb) as u16).collect();
        let b: Vec<u16> = k.iter().map(|&k| (k % bb) as u16).collect();
        let fixed = genes_and_reads(links, gene_of, mu, s, (&a, &b), held);
        let cells = cell_counts(links, 2, &a);
        let terms: Vec<f64> = grid
            .iter()
            .map(|&g| {
                ln_prior(g)
                    + cells
                        .iter()
                        .map(|c| urn(c, s.alpha, &[g, 1.0 - g]))
                        .sum::<f64>()
            })
            .collect();
        let mx = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
        let w: Vec<f64> = terms.iter().map(|t| (t - mx).exp()).collect();
        let mass: f64 = w.iter().sum();
        let own = |g: f64| if a[0] == 0 { g } else { 1.0 - g };
        weight_mean.push(grid.iter().zip(&w).map(|(&g, &x)| own(g) * x).sum::<f64>() / mass);
        out.push(((a, b), fixed + mx + mass.ln()));
    }
    let mx = out.iter().map(|x| x.1).fold(f64::NEG_INFINITY, f64::max);
    let total: f64 = out.iter().map(|x| (x.1 - mx).exp()).sum();
    let probs: Vec<(Labels, f64)> = out
        .into_iter()
        .map(|(z, l)| (z, (l - mx).exp() / total))
        .collect();
    let mean = probs
        .iter()
        .zip(&weight_mean)
        .map(|((_, p), m)| p * m)
        .sum();
    (probs, mean)
}

/// The mean over `draws` of the global weight of link 0's cell role.
fn own_weight(draws: &[(Vec<u16>, Vec<u16>, Vec<f64>)]) -> f64 {
    draws
        .iter()
        .map(|(a, _, beta)| beta[a[0] as usize])
        .sum::<f64>()
        / draws.len() as f64
}

/// The link preset with `cell_roles` roles and concentration `alpha`.
fn link_preset(cell_roles: usize, alpha: f64) -> Settings {
    Settings {
        cell_roles,
        gene_roles: 1,
        alpha,
        gene_alpha: 0.5,
        weights: Weights::Hdp(1.0),
        shift_sd: None,
        strength: 2.0,
        cell_sd: 1.0,
    }
}

/// The blocks preset with three cell roles and two gene roles.
fn blocks_preset() -> Settings {
    Settings {
        cell_roles: 3,
        gene_roles: 2,
        alpha: 0.5,
        gene_alpha: 0.5,
        weights: Weights::Symmetric,
        shift_sd: Some(1.5),
        strength: f64::INFINITY,
        cell_sd: 1.0,
    }
}
