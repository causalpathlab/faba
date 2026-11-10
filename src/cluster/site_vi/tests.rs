//! The site blockmodel's variational fit on small made-up reads: its
//! bound against the exact likelihood, its priors' terms, the ELBO that
//! coordinate ascent never lowers, and the moves that leave it as it is.

use super::*;
use crate::cluster::crp::ln_sigmoid;
use crate::cluster::read::{Reads, Weighted};
use rand::SeedableRng;

/// Two kinds of 30 cells over 6 genes of 3 sites each: at the first three
/// genes the first kind is mostly methylated and the second mostly not, at
/// the last three the other way round; 8 reads a site.
fn two_kinds() -> (Reads, Vec<u32>) {
    let cells: Vec<Vec<(u32, u32, u32)>> = (0..60u32)
        .map(|c| {
            (0..18u32)
                .map(|s| {
                    let high = (c < 30) == (s < 9);
                    let on = if high {
                        6 + (c + s) % 2
                    } else {
                        1 + (c * 3 + s) % 2
                    };
                    (s, on, 8 - on)
                })
                .collect()
        })
        .collect();
    let gene_of = (0..18u32).map(|s| s / 3).collect();
    (Reads::from_cells(cells), gene_of)
}

fn prior() -> Prior {
    Prior {
        cell_roles: 2,
        gene_roles: 2,
        alpha: 1.0,
        gene_alpha: 1.0,
        shift_sd: 1.0,
        cell_sd: 1.0,
    }
}

/// The baselines: each site's pooled fraction's logit.
fn baselines(links: &Links, sites: usize) -> Vec<f64> {
    let (mut on, mut all) = (vec![0.0; sites], vec![0.0; sites]);
    for l in 0..links.len() {
        let (s, m, u) = links.of(l);
        for i in 0..s.len() {
            on[s[i] as usize] += m[i];
            all[s[i] as usize] += m[i] + u[i];
        }
    }
    (0..sites)
        .map(|s| {
            let p = (on[s] + 0.5) / (all[s] + 1.0);
            (p / (1.0 - p)).ln()
        })
        .collect()
}

/// Each epoch's ELBO and the fit, from a k-means start.
fn fitted(options: &Options) -> (Vec<(usize, f64)>, Fitted, f64) {
    let (reads, gene_of) = two_kinds();
    let weighted = Weighted::new(&reads, None);
    let links = Links::new(&weighted, &gene_of);
    let mu = baselines(&links, 18);
    let data = Data::new(&links, &gene_of, 6, mu);
    let fit_of = Fit {
        data: &data,
        prior: prior(),
    };
    let mut rng = SmallRng::seed_from_u64(5);
    let start = kmeans_start(&fit_of, &Profiles::new(&fit_of), &mut rng, &|_| (), &|| {
        false
    })
    .unwrap();
    let seen = std::sync::Mutex::new(Vec::<(usize, f64)>::new());
    let epoch = |e: usize, f: f64| {
        let mut seen = seen.lock().unwrap();
        match seen.last_mut() {
            Some(last) if last.0 == e => last.1 = f,
            _ => seen.push((e, f)),
        }
    };
    let hooks = Hooks {
        epoch: &epoch,
        stop: &|| false,
        cancelled: &|| false,
    };
    let f = fit(&fit_of, options, start, &mut rng, &hooks).unwrap();
    let per = f.elbo_per_read(&data);
    (seen.into_inner().unwrap(), f, per)
}

/// The bound is below the exact log-likelihood (less its binomial
/// coefficient) at every `xi`, and equal at `xi = |eta|`.
#[test]
fn the_bound_is_below_the_likelihood_and_tight_at_its_xi() {
    for (m, u) in [(3.0, 5.0), (0.0, 2.0), (7.5, 0.5), (1.0, 1.0)] {
        let (n, kappa) = (m + u, m - (m + u) / 2.0);
        for eta in [-4.0, -1.0, 0.0, 0.3, 2.5] {
            let exact = m * ln_sigmoid(eta) + u * ln_sigmoid(-eta);
            for xi in [0.0, 0.5, 1.0, 3.0, 6.0] {
                let w = omega(n, xi);
                let bound =
                    kappa * eta - w / 2.0 * (eta * eta - xi * xi) + n * (ln_sigmoid(xi) - xi / 2.0);
                assert!(bound <= exact + 1e-12, "m {m} u {u} eta {eta} xi {xi}");
            }
            let xi = f64::abs(eta);
            let tight = kappa * eta + n * (ln_sigmoid(xi) - xi / 2.0);
            assert!((tight - exact).abs() < 1e-9, "tight at xi = |eta| = {xi}");
            let (w, bound) = jj(n, xi);
            assert!((w - omega(n, xi)).abs() < 1e-12 && (kappa * eta + bound - tight).abs() < 1e-9);
            assert!((w - n * (xi / 2.0).tanh() / (2.0 * xi.max(1e-300))).abs() < 1e-9 || xi < 1e-4);
        }
    }
    // E omega is continuous at 0.
    assert!((omega(2.0, 1e-5) - omega(2.0, 1.01e-4)).abs() < 1e-8);
}

/// A factor whose `q` is its prior adds nothing; one away from it, less.
#[test]
fn a_factor_at_its_prior_adds_nothing() {
    assert!(gaussian_term(0.0, 2.0, 2.0).abs() < 1e-12);
    assert!(gaussian_term(0.5, 2.0, 2.0) < 0.0 && gaussian_term(0.0, 1.0, 2.0) < 0.0);
    assert!(dirichlet_term(&[0.5, 0.5, 0.5], 0.5).0.abs() < 1e-12);
    assert!(dirichlet_term(&[2.5, 0.5, 0.5], 0.5).0 < 0.0);
}

/// Coordinate ascent, every cell in one batch, never lowers the ELBO.
#[test]
fn coordinate_ascent_never_lowers_the_elbo() {
    let options = Options {
        batch_size: 1000,
        epochs: 30,
        tol: 0.0,
        ..Options::default()
    };
    let (trace, f, _) = fitted(&options);
    assert_eq!((f.svi_epochs, f.epochs), (0, 30));
    for w in trace.windows(2) {
        let (a, b) = (w[0].1, w[1].1);
        assert!(
            b >= a - 1e-9 * a.abs(),
            "the ELBO fell from {a} to {b} at epoch {}",
            w[1].0
        );
    }
}

/// Stochastic epochs, then coordinate ascent to the tolerance; the two kinds
/// of cells end in different cell roles, and the shifts carry the pattern,
/// not the levels.
#[test]
fn the_fit_converges_and_tells_the_kinds_apart() {
    let options = Options {
        batch_size: 16,
        svi_epochs: 5,
        epochs: 200,
        tol: 1e-6,
        ..Options::default()
    };
    let (_, f, per) = fitted(&options);
    assert!(f.converged, "not converged in {} epochs", f.epochs);
    assert!(per.is_finite() && per < 0.0);
    let main = |c: usize| {
        let t = &f.locals.t[c * 2..c * 2 + 2];
        usize::from(t[1] > t[0])
    };
    let first: std::collections::BTreeSet<usize> = (0..30).map(main).collect();
    let second: std::collections::BTreeSet<usize> = (30..60).map(main).collect();
    assert_eq!((first.len(), second.len()), (1, 1), "{first:?} {second:?}");
    assert!(first.is_disjoint(&second));
}

/// Moving the shifts and the levels together by the prior's constant leaves
/// every cell's terms as they are and the priors' no worse.
#[test]
fn the_translation_never_lowers_the_elbo() {
    let (reads, gene_of) = two_kinds();
    let weighted = Weighted::new(&reads, None);
    let links = Links::new(&weighted, &gene_of);
    let data = Data::new(&links, &gene_of, 6, baselines(&links, 18));
    let fit_of = Fit {
        data: &data,
        prior: prior(),
    };
    let mut rng = SmallRng::seed_from_u64(1);
    let Start {
        mut globals,
        mut locals,
    } = kmeans_start(&fit_of, &Profiles::new(&fit_of), &mut rng, &|_| (), &|| {
        false
    })
    .unwrap();
    globals.d.iter_mut().for_each(|d| *d += 0.7);
    let all: Vec<usize> = (0..60).collect();
    let before = score(&fit_of, &all, &mut locals, &globals, &|| false).unwrap();
    fit_of.translate(&mut globals, &mut locals);
    let after = score(&fit_of, &all, &mut locals, &globals, &|| false).unwrap();
    assert!(after >= before - 1e-9 * before.abs(), "{before} to {after}");
}
