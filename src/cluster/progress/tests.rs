//! The wrapper around each chain: its bar is the chain's own whatever
//! order the chains start in, and a chain that starts after a cancel
//! skips the model's start.

use super::*;
use rand::SeedableRng;

/// A model that counts its starts, and draws nothing.
#[derive(Default)]
struct Starts(AtomicUsize);

impl McmcModel for Starts {
    type State = ();
    type Sample = ();
    type Result = ();

    fn init(&self, _: &mut SmallRng) {
        self.0.fetch_add(1, Ordering::SeqCst);
    }

    fn sweep(&self, _: &mut (), _: &mut SmallRng) {}

    fn collect(&self, _: &()) {}

    fn summarize(&self, _: Vec<()>) {}
}

/// The engine's random numbers for chain `i` of a run seeded `seed`.
fn chain_rng(seed: u64, i: usize) -> SmallRng {
    SmallRng::seed_from_u64(seed.wrapping_add(i as u64))
}

/// Chains started out of order each count on their own bar, as the
/// summary numbers them; the random numbers are left as they were.
#[test]
fn each_chain_counts_on_its_own_bar() {
    let p = Progress::new(3);
    let model = Starts::default();
    let watched = Watched {
        model: &model,
        progress: &p,
        seed: u64::MAX,
    };
    for i in [2, 0, 1] {
        let mut rng = chain_rng(u64::MAX, i);
        let mut st = watched.init(&mut rng);
        assert!(
            rng == chain_rng(u64::MAX, i),
            "chain {i}: the draws untouched"
        );
        for _ in 0..=i {
            watched.sweep(&mut st, &mut rng);
        }
    }
    assert_eq!((0..3).map(|i| p.done(i)).collect::<Vec<_>>(), [1, 2, 3]);
}

/// A chain that starts once the run is cancelled does not start its model.
#[test]
fn a_chain_started_after_a_cancel_skips_the_models_start() {
    let p = Progress::new(2);
    let model = Starts::default();
    let watched = Watched {
        model: &model,
        progress: &p,
        seed: 0,
    };
    let mut first = watched.init(&mut chain_rng(0, 0));
    p.cancel();
    let mut second = watched.init(&mut chain_rng(0, 1));
    assert_eq!(model.0.load(Ordering::SeqCst), 1);
    for st in [&mut first, &mut second] {
        watched.sweep(st, &mut chain_rng(0, 0));
        assert!(watched.collect(st).is_none());
    }
}

/// A schedule too long to count panics as the engine's does, not wrapping.
#[test]
#[should_panic(expected = "overflows usize")]
fn an_overflowing_schedule_panics_as_the_engines_does() {
    sweeps(&McmcConfig {
        n_samples: usize::MAX,
        warmup: 0,
        thin: 2,
        seed: 0,
    });
}
