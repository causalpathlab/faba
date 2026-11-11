//! The count model's variational fit on small made-up counts: a factor at
//! its prior, the ELBO that coordinate ascent never lowers, and the tables'
//! shares.

use super::*;
use rand::SeedableRng;

/// Two kinds of 40 cells over 24 genes: each kind's half of the genes at
/// 6 to 9 counts, the other half at 0 to 2.
fn two_kinds() -> Counts {
    let rows: Vec<Vec<(u32, u32)>> = (0..80)
        .map(|c| {
            (0..24u32)
                .map(|g| {
                    let high = (c < 40) == (g < 12);
                    let x = if high {
                        6 + (c as u32 * 7 + g * 3) % 4
                    } else {
                        (c as u32 + g) % 3
                    };
                    (g, x)
                })
                .collect()
        })
        .collect();
    Counts::from_cells(24, rows.iter().map(Vec::as_slice))
}

fn settings() -> Settings {
    Settings {
        cell_roles: 3,
        gene_roles: 3,
        ..Settings::default()
    }
}

/// A fit of `counts` from a k-means start: its ELBO by epoch and the fit.
fn fitted(counts: &Counts, options: &Options) -> (Vec<(usize, f64)>, Fitted) {
    let mut rng = SmallRng::seed_from_u64(3);
    let profiles = Profiles::new(counts);
    let start = kmeans_start(counts, &profiles, &settings(), &mut rng, &|_| (), &|| false).unwrap();
    // Each epoch's ELBO, an epoch told again replacing its estimate.
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
    let f = fit(counts, &settings(), options, start, &mut rng, &hooks).unwrap();
    (seen.into_inner().unwrap(), f)
}

/// A factor whose `q` is its prior adds nothing to the ELBO; one away from
/// it, less.
#[test]
fn a_factor_at_its_prior_adds_nothing() {
    for (a, b) in [(0.3, 1.0), (1.0, 1.0), (0.1, 0.1), (5.0, 2.0)] {
        assert!(kl_fixed(a, b, a, b).abs() < 1e-12, "Gamma({a}, {b})");
        assert!(kl_fixed(a, b, a + 1.0, b) < 0.0);
        assert!(kl_fixed(a, b, a, 2.0 * b) < 0.0);
    }
}

/// Coordinate ascent, every cell in one batch, never lowers the ELBO.
#[test]
fn coordinate_ascent_never_lowers_the_elbo() {
    let counts = two_kinds();
    let options = Options {
        batch_size: 1000,
        epochs: 40,
        tol: 0.0,
        ..Options::default()
    };
    let (trace, f) = fitted(&counts, &options);
    assert_eq!(f.ran.svi_epochs, 0);
    assert_eq!(f.ran.epochs, 40);
    assert!(trace.len() >= 40, "{trace:?}");
    for w in trace.windows(2) {
        let (a, b) = (w[0].1, w[1].1);
        assert!(
            b >= a - 1e-9 * a.abs(),
            "the ELBO fell from {a} to {b} at epoch {}",
            w[1].0
        );
    }
    assert!((trace.last().unwrap().1 - f.elbo).abs() <= 1e-9 * f.elbo.abs());
}

/// With batches fewer than every cell, stochastic epochs come first, then
/// coordinate ascent to the tolerance; the fit ends close to coordinate
/// ascent's alone.
#[test]
fn stochastic_epochs_then_coordinate_ascent_converge() {
    let counts = two_kinds();
    let options = Options {
        batch_size: 16,
        svi_epochs: 5,
        epochs: 300,
        tol: 1e-4,
        ..Options::default()
    };
    let (trace, f) = fitted(&counts, &options);
    assert_eq!(f.ran.svi_epochs, 5);
    assert!(f.ran.converged, "not converged in {} epochs", f.ran.epochs);
    assert!(trace.iter().all(|p| p.1.is_finite()));
    let cavi = fitted(
        &counts,
        &Options {
            batch_size: 1000,
            ..options
        },
    )
    .1;
    let per = |x: &Fitted| x.elbo_per_count(&counts);
    assert!(
        (per(&f) - per(&cavi)).abs() < 0.05,
        "{} against {}",
        per(&f),
        per(&cavi)
    );
}

/// The tables' shares are shares, and the two kinds of cells come apart:
/// no main role holds cells of both kinds.
#[test]
fn the_summary_gives_shares_and_tells_the_kinds_apart() {
    let counts = two_kinds();
    let options = Options {
        batch_size: 1000,
        epochs: 200,
        ..Options::default()
    };
    let (_, f) = fitted(&counts, &options);
    let s = summarize(&f);
    for row in s.cell.chunks(3).chain(s.gene.chunks(3)) {
        assert!((row.iter().sum::<f32>() - 1.0).abs() < 1e-4, "{row:?}");
    }
    assert_eq!((s.log_rate.len(), s.units.len()), (9, 9));
    let main = |c: usize| {
        let row = &s.cell[c * 3..(c + 1) * 3];
        (0..3).max_by(|&a, &b| row[a].total_cmp(&row[b])).unwrap()
    };
    let first: std::collections::BTreeSet<usize> = (0..40).map(main).collect();
    let second: std::collections::BTreeSet<usize> = (40..80).map(main).collect();
    assert!(first.is_disjoint(&second), "{first:?} {second:?}");
}
