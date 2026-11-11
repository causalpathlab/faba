//! The shared schedule on a stub fit whose ELBO is set in advance: its
//! stochastic batches and step sizes, its epochs told, and where it stops.

use super::*;
use rand::SeedableRng;

/// A stub fit: each coordinate-ascent epoch's ELBO from `elbos`, each step
/// noted.
#[derive(Default)]
struct Stub {
    elbos: Vec<f64>,
    epochs: usize,
    /// Each batch's cells, scale and rho; the stochastic epochs ended; the
    /// global steps of coordinate ascent; whether the last pass ran.
    batches: Vec<(Vec<usize>, f64, f64)>,
    svi_ended: usize,
    steps: usize,
    last: bool,
}

impl Steps for Stub {
    type Pass = f64;
    type End = f64;
    fn units(&self) -> f64 {
        10.0
    }
    fn first(&mut self) -> Option<()> {
        Some(())
    }
    fn batch(&mut self, cells: &[usize], scale: f64, rho: f64) -> Option<()> {
        self.batches.push((cells.to_vec(), scale, rho));
        Some(())
    }
    fn end_svi_epoch(&mut self) {
        self.svi_ended += 1;
    }
    fn score(&mut self) -> Option<f64> {
        Some(-100.0)
    }
    fn epoch(&mut self) -> Option<(f64, f64)> {
        let f = self.elbos[self.epochs.min(self.elbos.len() - 1)];
        self.epochs += 1;
        Some((f, f + 0.5))
    }
    fn converged(&mut self, pass: f64) -> f64 {
        pass
    }
    fn step(&mut self, _: f64) {
        self.steps += 1;
    }
    fn last(&mut self) -> Option<f64> {
        self.last = true;
        Some(-1.0)
    }
    fn elbo(end: &f64) -> f64 {
        *end
    }
}

/// The schedule `options` give `stub` over `cells` cells: its end, how it
/// ran, and every epoch told.
fn ran(stub: &mut Stub, cells: usize, options: &Options) -> (f64, Ran, Vec<(usize, f64)>) {
    let told = std::sync::Mutex::new(Vec::new());
    let epoch = |e: usize, f: f64| told.lock().unwrap().push((e, f));
    let no = || false;
    let hooks = Hooks {
        epoch: &epoch,
        stop: &no,
        cancelled: &no,
    };
    let mut rng = SmallRng::seed_from_u64(1);
    let (end, ran) = run(stub, cells, options, &mut rng, &hooks).unwrap();
    (end, ran, told.into_inner().unwrap())
}

/// Ten cells in batches of four: three batches (4, 4, 2) an epoch, each
/// scaled to every cell, the step size `(tau0 + t)^-forget` with
/// `tau0` the batches; then coordinate ascent until an epoch moves the
/// ELBO by less than the tolerance per unit, its end the pass it stopped
/// on, no last pass.
#[test]
fn stochastic_epochs_then_ascent_to_the_tolerance() {
    let options = Options {
        svi_epochs: 2,
        batch_size: 4,
        epochs: 20,
        tol: 0.01,
        forget: 0.7,
    };
    let mut stub = Stub {
        elbos: vec![-50.0, -20.0, -19.95],
        ..Stub::default()
    };
    let (end, how, told) = ran(&mut stub, 10, &options);
    let sizes: Vec<usize> = stub.batches.iter().map(|b| b.0.len()).collect();
    assert_eq!(sizes, [4, 4, 2, 4, 4, 2]);
    for (t, (cells, scale, rho)) in stub.batches.iter().enumerate() {
        assert_eq!(*scale, 10.0 / cells.len() as f64);
        assert!((rho - (3.0 + t as f64).powf(-0.7)).abs() < 1e-12);
    }
    let mut seen: Vec<usize> = stub.batches[..3].iter().flat_map(|b| b.0.clone()).collect();
    seen.sort_unstable();
    assert_eq!(seen, (0..10).collect::<Vec<_>>(), "every cell an epoch");
    assert_eq!(stub.svi_ended, 2);
    // -20 to -19.95 is 0.005 a unit: converged on the third epoch of
    // ascent, two global steps taken before it.
    assert_eq!((stub.steps, stub.last, end), (2, false, -19.45));
    assert_eq!(
        how,
        Ran {
            svi_epochs: 2,
            epochs: 4,
            converged: true
        }
    );
    // The first stochastic epoch scored, the second by the ascent's first
    // epoch; each ascent epoch told at its start; the end told last.
    assert_eq!(
        told,
        [
            (1, -100.0),
            (2, -50.0),
            (3, -20.0),
            (4, -19.95),
            (4, -19.45)
        ]
    );
}

/// One batch holding every cell is plain coordinate ascent; not converged
/// within the epochs asked, the fit ends with a last pass.
#[test]
fn one_batch_is_ascent_alone_and_unsettled_ends_with_a_last_pass() {
    let options = Options {
        svi_epochs: 5,
        batch_size: 100,
        epochs: 3,
        tol: 1e-9,
        forget: 0.7,
    };
    let mut stub = Stub {
        elbos: vec![-50.0, -40.0, -30.0],
        ..Stub::default()
    };
    let (end, how, told) = ran(&mut stub, 10, &options);
    assert!(stub.batches.is_empty());
    assert_eq!((stub.steps, stub.last, end), (3, true, -1.0));
    assert_eq!(
        how,
        Ran {
            svi_epochs: 0,
            epochs: 3,
            converged: false
        }
    );
    assert_eq!(told, [(0, -50.0), (1, -40.0), (2, -30.0), (3, -1.0)]);
}
