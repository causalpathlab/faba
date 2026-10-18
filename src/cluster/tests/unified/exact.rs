//! The moves against posteriors enumerated from the model's definition:
//! the labels given shifts and levels, the swaps' ratios, the log joint
//! and the held-out score.

use super::*;

/// The label moves (link relabels, cell and gene swaps, the global weights
/// under HDP), the shifts and the cells' levels held, match the enumerated
/// posterior on the tiny fixture, as read and fractionally weighted, for
/// every allowed switch combination (the levels unequal at infinite
/// strength); under HDP, so does the mean global weight of link 0's role.
#[test]
fn every_switch_matches_the_exact_posterior_given_shifts_and_levels() {
    for (links, s) in tiny_fixtures()
        .iter()
        .flat_map(|links| combinations().into_iter().map(move |s| (links, s)))
    {
        let model = Model::new(links, &TINY_GENES, 3, s.clone()).expect("an allowed setting");
        let (delta, nu) = (fixed_shifts(&s), fixed_levels(&s, links.cells));
        let (exact, weight) = exact(links, &TINY_GENES, model.baseline(), &s, (&delta, &nu));
        let draws: Vec<_> = [131, 132]
            .iter()
            .flat_map(|&seed| model.label_draws((&delta, &nu), 500, 40_000, seed))
            .collect();
        let labels: Vec<Labels> = draws
            .iter()
            .map(|(a, b, _)| (a.clone(), b.clone()))
            .collect();
        let tv = label_tv(labels, &exact);
        assert!(tv < 0.02, "{s:?}: total variation {tv:.4}");
        if let Weights::Hdp(_) = s.weights {
            let sampled = own_weight(&draws);
            assert!(
                (sampled - weight).abs() < 0.01,
                "{s:?}: E[beta of link 0's role] sampled {sampled:.4}, exact {weight:.4}"
            );
        }
    }
}

/// Every pair of `roles` roles.
fn pairs(roles: usize) -> Vec<(u16, u16)> {
    (0..roles as u16)
        .flat_map(|r| (r + 1..roles as u16).map(move |r2| (r, r2)))
        .collect()
}

/// The cell swap's and the gene swap's Metropolis-Hastings ratio is exactly
/// the change in log joint, for every labelling, unit and pair of roles,
/// under HDP and symmetric weights, finite strength and infinite (the
/// cells' levels unequal), over reads as read and fractionally weighted.
#[test]
fn swap_ratios_are_the_change_in_log_joint() {
    let reads = Reads::from_cells(vec![
        vec![(0, 3, 0), (1, 0, 2), (2, 1, 1), (3, 2, 2)],
        vec![(0, 0, 3), (1, 2, 0), (2, 2, 0), (3, 0, 1)],
    ]);
    let gene_of = [0, 0, 1, 2];
    let (on, off) = (
        [0.6, 1.0, 1.4, 0.7, 1.0, 0.8, 1.0, 1.3],
        [1.0, 1.7, 0.5, 1.1, 0.9, 1.0, 1.2, 0.6],
    );
    let fixtures = [
        Links::new(&Weighted::new(&reads, None), &gene_of),
        Links::new(&Weighted::new(&reads, Some((&on, &off))), &gene_of),
    ];
    let n = fixtures[0].len();
    let mut checked = [0usize; 2];
    for (links, weights) in fixtures
        .iter()
        .flat_map(|links| [Weights::Hdp(1.0), Weights::Symmetric].map(|w| (links, w)))
    {
        for strength in [2.0, f64::INFINITY] {
            for (cell_roles, gene_roles, shift) in [(2, 2, true), (3, 1, strength.is_infinite())] {
                let s = Settings {
                    cell_roles,
                    gene_roles,
                    alpha: 0.7,
                    gene_alpha: 1.3,
                    weights,
                    shift_sd: shift.then_some(1.0),
                    strength,
                    cell_sd: 0.8,
                };
                let nu = fixed_levels(&s, links.cells);
                let model = Model::new(links, &gene_of, 4, s).expect("an allowed setting");
                let beta: Vec<f64> = match (weights, cell_roles) {
                    (Weights::Hdp(_), 2) => vec![0.7, 0.3],
                    (Weights::Hdp(_), _) => vec![0.5, 0.3, 0.2],
                    (Weights::Symmetric, r) => vec![1.0 / r as f64; r],
                };
                let kk = cell_roles * gene_roles;
                let delta: Vec<f64> = if shift {
                    vec![0.8, -0.5, -1.1, 0.4]
                } else {
                    vec![0.0; kk]
                };
                for code in 0..kk.pow(n as u32) {
                    let k: Vec<usize> = (0..n).map(|l| code / kk.pow(l as u32) % kk).collect();
                    let a: Vec<u16> = k.iter().map(|&k| (k / gene_roles) as u16).collect();
                    let b: Vec<u16> = k.iter().map(|&k| (k % gene_roles) as u16).collect();
                    let units = [(2, cell_roles, true), (3, gene_roles, false)];
                    for (count, roles, cell) in units {
                        for unit in 0..count {
                            for pair in pairs(roles) {
                                let held = Held {
                                    beta: &beta,
                                    delta: &delta,
                                    nu: &nu,
                                };
                                let Some((ratio, change)) =
                                    model.swap_check((&a, &b), &held, unit, cell, pair)
                                else {
                                    continue;
                                };
                                assert!(
                                    (ratio - change).abs() < 1e-9,
                                    "{weights:?} strength {strength} a {a:?} b {b:?} \
                                     {} {unit} swap {pair:?}: ratio {ratio}, log joint {change}",
                                    if cell { "cell" } else { "gene" }
                                );
                                checked[usize::from(!cell)] += 1;
                            }
                        }
                    }
                }
            }
        }
    }
    assert!(
        checked.iter().all(|&c| c > 1000),
        "swaps checked {checked:?}"
    );
}

/// The labels of `n` links in `K` blocks numbered `code` (link `l`'s block
/// digit `l` of `code` base `K`, block `a B + b`).
fn labelling(code: usize, n: usize, s: &Settings) -> Labels {
    let (bb, kk) = (s.gene_roles, s.cell_roles * s.gene_roles);
    let k: Vec<usize> = (0..n).map(|l| code / kk.pow(l as u32) % kk).collect();
    (
        k.iter().map(|&k| (k / bb) as u16).collect(),
        k.iter().map(|&k| (k % bb) as u16).collect(),
    )
}

/// Global weights to hold the labels at: unequal under HDP.
fn fixed_weights(s: &Settings) -> Vec<f64> {
    match s.weights {
        Weights::Hdp(_) => vec![0.7, 0.3],
        Weights::Symmetric => vec![0.5, 0.5],
    }
}

/// The log joint the chains report (which R-hat and the choice of the best
/// chain read) is the model's definition, for every labelling of the tiny
/// fixture (as read and fractionally weighted) at two sets of shifts under
/// every allowed switch combination:
/// each cell's urn at the global weights, each gene's, each block's sites'
/// marginals (at infinite strength each link's reads, its cell's level
/// in), the shifts' Normal prior and the levels' (binomial coefficients
/// and the Normals' normalisers left out of both).
#[test]
fn the_log_joint_is_the_models_definition() {
    for (links, s) in tiny_fixtures()
        .iter()
        .flat_map(|links| combinations().into_iter().map(move |s| (links, s)))
    {
        let model = Model::new(links, &TINY_GENES, 3, s.clone()).expect("an allowed setting");
        let mu = model.baseline().to_vec();
        let beta = fixed_weights(&s);
        let (fixed, nu) = (fixed_shifts(&s), fixed_levels(&s, links.cells));
        let held = [
            (fixed.iter().map(|d| -0.6 * d).collect(), nu.clone()),
            (fixed, nu.iter().map(|v| -1.5 * v).collect::<Vec<f64>>()),
        ];
        let kk = s.cell_roles * s.gene_roles;
        for code in 0..kk.pow(links.len() as u32) {
            let (a, b) = labelling(code, links.len(), &s);
            for (delta, nu) in &held {
                let prior = s.shift_sd.map_or(0.0, |sd| {
                    delta.iter().map(|d| -d * d / (2.0 * sd * sd)).sum()
                });
                let defined = genes_and_reads(links, &TINY_GENES, &mu, &s, (&a, &b), (delta, nu))
                    + cell_counts(links, 2, &a)
                        .iter()
                        .map(|c| urn(c, s.alpha, &beta))
                        .sum::<f64>()
                    + prior
                    + level_prior(&s, nu);
                let at = Held {
                    beta: &beta,
                    delta,
                    nu,
                };
                let reported = model.log_joint_at((&a, &b), &at);
                assert!(
                    (reported - defined).abs() < 1e-9,
                    "{s:?} a {a:?} b {b:?} delta {delta:?} nu {nu:?}: reported {reported}, \
                     defined {defined}"
                );
            }
        }
    }
}

/// Held-out reads are scored at the model's predictive, computed here from
/// its definition, for every labelling of the tiny fixture under every
/// allowed switch combination: a read at site `s` of gene `g` in cell `c`
/// converted with probability `sum_ab theta_ca psi_gb p_s,ab`, the role
/// weights their posterior means given the labels and `p` the block's
/// posterior-mean profile (at infinite strength `sigmoid(mu_s + nu_c +
/// delta_ab)`, the cell's own level in); an unconverted read with
/// `sum_ab theta_ca psi_gb (1 - p_s,ab)`.
#[test]
fn held_out_reads_are_scored_at_the_predictive() {
    let links = tiny_links();
    let held = Reads::from_cells(vec![vec![(0, 2, 1), (2, 0, 3)], vec![(1, 1, 0), (2, 4, 1)]]);
    let held_reads = 12.0;
    let link_gene = links.genes_of(&TINY_GENES);
    for s in combinations() {
        let model = Model::new(&links, &TINY_GENES, 3, s.clone())
            .expect("an allowed setting")
            .with_heldout(&held);
        let mu = model.baseline().to_vec();
        let (beta, delta) = (fixed_weights(&s), fixed_shifts(&s));
        let nu = fixed_levels(&s, links.cells);
        let (aa, bb) = (s.cell_roles, s.gene_roles);
        let kk = aa * bb;
        for code in 0..kk.pow(links.len() as u32) {
            let (a, b) = labelling(code, links.len(), &s);
            let theta: Vec<Vec<f64>> = cell_counts(&links, aa, &a)
                .iter()
                .map(|n| {
                    let total = n.iter().sum::<f64>() + s.alpha;
                    (0..aa)
                        .map(|r| (n[r] + s.alpha * beta[r]) / total)
                        .collect()
                })
                .collect();
            let psi: Vec<Vec<f64>> = (0..2)
                .map(|g| {
                    let mut m = vec![0.0; bb];
                    (0..links.len())
                        .filter(|&l| link_gene[l] as usize == g)
                        .for_each(|l| m[b[l] as usize] += 1.0);
                    let total = m.iter().sum::<f64>() + s.gene_alpha;
                    m.iter()
                        .map(|x| (x + s.gene_alpha / bb as f64) / total)
                        .collect()
                })
                .collect();
            let (on, off) = block_reads(&links, &s, (&a, &b), mu.len());
            let p = |c: usize, j: usize, k: usize| {
                if s.strength.is_finite() {
                    let at = j * kk + k;
                    let centre = sigmoid(mu[j] + delta[k]);
                    (s.strength * centre + on[at]) / (s.strength + on[at] + off[at])
                } else {
                    sigmoid(mu[j] + nu[c] + delta[k])
                }
            };
            let mut defined = 0.0;
            for (c, theta) in theta.iter().enumerate() {
                let (sites, m, u) = held.of(c);
                for i in 0..sites.len() {
                    let (j, g) = (sites[i] as usize, TINY_GENES[sites[i] as usize] as usize);
                    let mix = |f: &dyn Fn(f64) -> f64| -> f64 {
                        (0..kk)
                            .map(|k| theta[k / bb] * psi[g][k % bb] * f(p(c, j, k)))
                            .sum()
                    };
                    defined +=
                        m[i] as f64 * mix(&|p| p).ln() + u[i] as f64 * mix(&|p| 1.0 - p).ln();
                }
            }
            let at = Held {
                beta: &beta,
                delta: &delta,
                nu: &nu,
            };
            let (scored, count) = model.heldout_at((&a, &b), &at);
            assert_eq!(count, held_reads);
            assert!(
                (scored - defined).abs() < 1e-9,
                "{s:?} a {a:?} b {b:?}: scored {scored}, defined {defined}"
            );
        }
    }
}
