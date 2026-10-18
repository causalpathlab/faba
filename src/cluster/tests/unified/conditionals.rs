//! The shifts, the cells' levels and the global weights, each drawn with
//! the rest held, against its conditional.

use super::*;

/// Each block's shift update, labels held, leaves its conditional invariant
/// at finite and infinite strength: draws of `delta[0]` match the
/// conditional's mean and spread (numerical integral of its Normal prior
/// times its block's site marginals, at infinite strength each of its
/// links' Binomial reads at its cell's level), with many reads in the
/// block and with few, where the prior counts.
#[test]
fn shifts_follow_their_conditional() {
    for strength in [2.0, f64::INFINITY] {
        check_shift_conditional(90, 20, 1.0, strength);
        check_shift_conditional(6, 2, 0.5, strength);
    }
}

/// The conditional's mean and sd on a grid over `[-4, 4]` of log density
/// `f`.
fn grid_moments(f: impl Fn(f64) -> f64) -> (f64, f64) {
    let grid: Vec<f64> = (0..8000)
        .map(|i| -4.0 + 8.0 * (i as f64 + 0.5) / 8000.0)
        .collect();
    let terms: Vec<f64> = grid.iter().map(|&x| f(x)).collect();
    let mx = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let w: Vec<f64> = terms.iter().map(|t| (t - mx).exp()).collect();
    let mass: f64 = w.iter().sum();
    let mean = grid.iter().zip(&w).map(|(x, w)| x * w).sum::<f64>() / mass;
    let var = grid
        .iter()
        .zip(&w)
        .map(|(x, w)| (x - mean).powi(2) * w)
        .sum::<f64>()
        / mass;
    (mean, var.sqrt())
}

/// The mean and sd of `draws` after the first 200.
fn draw_moments(draws: Vec<f64>) -> (f64, f64) {
    let kept = &draws[200..];
    let n = kept.len() as f64;
    let mean = kept.iter().sum::<f64>() / n;
    let var = kept.iter().map(|x| (x - mean).powi(2)).sum::<f64>() / n;
    (mean, var.sqrt())
}

fn check_shift_conditional(cells: usize, genes: usize, shift_sd: f64, strength: f64) {
    let (reads, gene_of, cell_role, gene_role, level) = planted_blocks(cells, genes, 1.0, 117);
    let links = Links::new(&Weighted::new(&reads, None), &gene_of);
    let s = Settings {
        cell_roles: 3,
        gene_roles: 2,
        alpha: 1.0,
        gene_alpha: 1.0,
        weights: Weights::Symmetric,
        shift_sd: Some(shift_sd),
        strength,
        cell_sd: 1.0,
    };
    let model = Model::new(&links, &gene_of, gene_of.len(), s).expect("an allowed setting");
    let mu = model.baseline().to_vec();
    let link_gene = links.genes_of(&gene_of);
    let a: Vec<u16> = links.cell.iter().map(|&c| cell_role[c as usize]).collect();
    let b: Vec<u16> = link_gene.iter().map(|&g| gene_role[g as usize]).collect();
    // The planted levels held at infinite strength; none at finite.
    let nu: Vec<f64> = level
        .iter()
        .map(|&v| if strength.is_finite() { 0.0 } else { v })
        .collect();
    // Block 0's reads: per site at finite strength, per link entry (with
    // its cell's level) at infinite.
    let (mut on, mut off) = (vec![0.0; mu.len()], vec![0.0; mu.len()]);
    let mut entries: Vec<(f64, f64, f64)> = Vec::new();
    for l in (0..links.len()).filter(|&l| a[l] == 0 && b[l] == 0) {
        let (sites, m, u) = links.of(l);
        for i in 0..sites.len() {
            on[sites[i] as usize] += m[i];
            off[sites[i] as usize] += u[i];
            entries.push((
                mu[sites[i] as usize] + nu[links.cell[l] as usize],
                m[i],
                u[i],
            ));
        }
    }
    let fast = fast_sites(&links, mu.len());
    let (exact, exact_sd) = grid_moments(|d| {
        let reads: f64 = if strength.is_finite() {
            (0..mu.len())
                .filter(|&j| on[j] + off[j] > 0.0)
                .map(|j| site_marginal(sigmoid(mu[j] + d), on[j], off[j], strength, fast[j]))
                .sum()
        } else {
            entries
                .iter()
                .map(|&(x, m, u)| site_marginal(sigmoid(x + d), m, u, strength, false))
                .sum()
        };
        -d * d / (2.0 * shift_sd * shift_sd) + reads
    });
    let (sampled, sd) = draw_moments(model.shift_draws((&a, &b), &nu, 4000, 119));
    assert!(
        (sampled - exact).abs() < 0.1 * exact_sd && (sd / exact_sd - 1.0).abs() < 0.1,
        "{cells} cells, strength {strength}: delta[0] | labels: sampled mean {sampled:.4} \
         sd {sd:.4}, exact mean {exact:.4} sd {exact_sd:.4}"
    );
}

/// Each cell's level update, labels and shifts held, leaves its conditional
/// invariant: draws of cell 0's level match the conditional's mean and
/// spread (numerical integral of its Normal prior times its links'
/// Binomial reads at their blocks' shifts), with many reads and with few,
/// where the prior counts.
#[test]
fn cell_levels_follow_their_conditional() {
    check_level_conditional(0.5, 1.0);
    check_level_conditional(0.05, 0.5);
}

fn check_level_conditional(cover: f64, cell_sd: f64) {
    let mut rng = SmallRng::seed_from_u64(125);
    let fracs: Vec<Vec<f64>> = (0..4)
        .map(|_| (0..40).map(|_| rng.random_range(0.1..0.6)).collect())
        .collect();
    let refs: Vec<&[f64]> = fracs.iter().map(Vec::as_slice).collect();
    let reads = draw_cells(&refs, cover, &mut rng);
    let gene_of: Vec<u32> = (0..40u32).map(|s| s / 4).collect();
    let links = Links::new(&Weighted::new(&reads, None), &gene_of);
    let shifts = vec![0.8, -0.5, 0.3, -1.1];
    let s = Settings {
        cell_roles: 2,
        gene_roles: 2,
        alpha: 1.0,
        gene_alpha: 1.0,
        weights: Weights::Symmetric,
        shift_sd: Some(1.0),
        strength: f64::INFINITY,
        cell_sd,
    };
    let model = Model::new(&links, &gene_of, 40, s).expect("an allowed setting");
    let mu = model.baseline().to_vec();
    let link_gene = links.genes_of(&gene_of);
    let a: Vec<u16> = (0..links.len()).map(|l| (l % 2) as u16).collect();
    let b: Vec<u16> = link_gene.iter().map(|&g| (g % 2) as u16).collect();
    let entries: Vec<(f64, f64, f64)> = (0..links.len())
        .filter(|&l| links.cell[l] == 0)
        .flat_map(|l| {
            let shift = shifts[a[l] as usize * 2 + b[l] as usize];
            let (sites, m, u) = links.of(l);
            (0..sites.len())
                .map(|i| (mu[sites[i] as usize] + shift, m[i], u[i]))
                .collect::<Vec<_>>()
        })
        .collect();
    let (exact, exact_sd) = grid_moments(|v| {
        -v * v / (2.0 * cell_sd * cell_sd)
            + entries
                .iter()
                .map(|&(x, m, u)| {
                    let p = sigmoid(x + v);
                    m * p.ln() + u * (1.0 - p).ln()
                })
                .sum::<f64>()
    });
    let (sampled, sd) = draw_moments(model.level_draws((&a, &b), &shifts, 4000, 127));
    assert!(
        (sampled - exact).abs() < 0.1 * exact_sd && (sd / exact_sd - 1.0).abs() < 0.1,
        "cover {cover}: cell 0's level | rest: sampled mean {sampled:.4} sd {sd:.4}, \
         exact mean {exact:.4} sd {exact_sd:.4}"
    );
}

/// The global weights' update, labels held fixed, leaves `p(beta | labels)`
/// invariant: draws of `beta_1` match its exact conditional mean and spread
/// (numerical integral over the log joint, under a flat prior at gamma 2)
/// for a large labelling, where a wrong table count would show.
#[test]
fn global_weights_follow_their_conditional() {
    let (reads, gene_of, truth) = planted_links(300, 40, 91);
    let links = Links::new(&Weighted::new(&reads, None), &gene_of);
    let s = Settings {
        weights: Weights::Hdp(2.0),
        ..link_preset(2, 0.5)
    };
    let model = Model::new(&links, &gene_of, gene_of.len(), s).unwrap();
    // Cells of pattern 0 in role 0, the rest in 1: unequal sizes.
    let a: Vec<u16> = links
        .cell
        .iter()
        .map(|&c| u16::from(truth[c as usize] != 0))
        .collect();
    let b = vec![0u16; a.len()];
    let grid: Vec<f64> = (1..4000).map(|i| i as f64 / 4000.0).collect();
    let terms: Vec<f64> = grid
        .iter()
        .map(|&x| {
            let at = Held {
                beta: &[x, 1.0 - x],
                delta: &[0.0, 0.0],
                nu: &vec![0.0; links.cells],
            };
            model.log_joint_at((&a, &b), &at)
        })
        .collect();
    let mx = terms.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let w: Vec<f64> = terms.iter().map(|t| (t - mx).exp()).collect();
    let mass: f64 = w.iter().sum();
    let exact = grid.iter().zip(&w).map(|(x, w)| x * w).sum::<f64>() / mass;
    let exact_sd = (grid
        .iter()
        .zip(&w)
        .map(|(x, w)| (x - exact).powi(2) * w)
        .sum::<f64>()
        / mass)
        .sqrt();
    let draws: Vec<f64> = model
        .beta_draws(&a, &b, 4000, 93)
        .into_iter()
        .skip(200)
        .collect();
    let sampled = draws.iter().sum::<f64>() / draws.len() as f64;
    let sd = (draws.iter().map(|x| (x - sampled).powi(2)).sum::<f64>() / draws.len() as f64).sqrt();
    // The mean alone misses an update off by a common factor: the spread
    // shows it.
    assert!(
        (sampled - exact).abs() < 0.01 && (sd / exact_sd - 1.0).abs() < 0.15,
        "beta_1 | labels: sampled mean {sampled:.4} sd {sd:.4}, exact mean {exact:.4} sd {exact_sd:.4}"
    );
}
