//! Sampling progress: counting the sweeps and tracing the log joint leave
//! every chain's draws as they were, the count is the engine's schedule, a
//! cancel turns the sweeps left into no-ops, and a cancelled fit writes
//! nothing.

use super::unified::{planted_blocks, planted_site_data, settings};
use super::{config, small_cells};
use crate::cluster::progress::{self, Cancelled, Progress};
use crate::cluster::read::Weighted;
use crate::cluster::unified::{Links, Model, Settings, Weights};
use crate::cluster::{crp, FitData};
use legume_numeric::mcmc::engine::{run_mcmc_parallel, McmcConfig, McmcModel};
use rand::rngs::SmallRng;
use std::sync::atomic::{AtomicUsize, Ordering};

/// A model that only counts its sweeps, and cancels `progress` at sweep
/// `at` (counted over every chain).
struct Counting<'p> {
    sweeps: AtomicUsize,
    at: usize,
    progress: &'p Progress,
}

impl<'p> Counting<'p> {
    fn new(at: usize, progress: &'p Progress) -> Self {
        Self {
            sweeps: AtomicUsize::new(0),
            at,
            progress,
        }
    }
}

impl McmcModel for Counting<'_> {
    type State = ();
    type Sample = ();
    type Result = usize;

    fn init(&self, _: &mut SmallRng) {}

    fn sweep(&self, _: &mut (), _: &mut SmallRng) {
        if self.sweeps.fetch_add(1, Ordering::SeqCst) + 1 == self.at {
            self.progress.cancel();
        }
    }

    fn collect(&self, _: &()) {}

    fn summarize(&self, kept: Vec<()>) -> usize {
        kept.len()
    }
}

impl progress::Traced for Counting<'_> {}

/// Each chain counts the sweeps the engine runs, and keeps the draws it
/// would keep, whatever the schedule.
#[test]
fn each_chain_counts_the_engines_sweeps() {
    for (n, w, thin) in [(0, 5, 1), (4, 3, 2), (3, 2, 0), (1, 0, 1)] {
        let cfg = McmcConfig {
            thin,
            ..config(n, w, 1)
        };
        let p = Progress::new(3);
        let model = Counting::new(usize::MAX, &p);
        let kept = progress::run(&model, &cfg, 3, &p).unwrap();
        assert_eq!(kept, run_mcmc_parallel(&model, &cfg, 3), "{cfg:?}");
        let per_chain = model.sweeps.load(Ordering::SeqCst) / 6;
        assert_eq!(progress::sweeps(&cfg), per_chain, "{cfg:?}");
        for i in 0..3 {
            assert_eq!(p.done(i), per_chain, "chain {i}, {cfg:?}");
            assert!(p.running(i).is_some(), "chain {i} started");
        }
        let (total, warmup) = p.schedule().unwrap();
        assert_eq!((total, warmup), (per_chain, w));
    }
}

/// A cancel mid-run: every chain stops at its next sweep, the run is
/// [`progress::Cancelled`], and no draw is kept.
#[test]
fn a_cancel_stops_every_chain_at_its_next_sweep() {
    let p = Progress::new(3);
    let model = Counting::new(5, &p);
    let cfg = config(100_000, 100_000, 1);
    let out = progress::run(&model, &cfg, 3, &p);
    assert!(out.unwrap_err().is::<Cancelled>());
    let swept = model.sweeps.load(Ordering::SeqCst);
    assert!(swept <= 5 + 3, "{swept} sweeps after a cancel at the 5th");
    // A run cancelled before it starts sweeps none.
    let model = Counting::new(usize::MAX, &p);
    assert!(progress::run(&model, &cfg, 3, &p)
        .unwrap_err()
        .is::<Cancelled>());
    assert_eq!(model.sweeps.load(Ordering::SeqCst), 0);
}

/// Once the fit commits to writing, a cancel is refused; a cancel first
/// makes the commit fail.
#[test]
fn a_cancel_and_a_commit_exclude_each_other() {
    let p = Progress::new(1);
    p.summarised(1.0).unwrap();
    assert!(!p.cancel() && !p.cancelled());
    let p = Progress::new(1);
    assert!(p.cancel() && p.cancelled());
    assert!(p.summarised(1.0).unwrap_err().is::<Cancelled>());
}

/// The blockmodel's settings this file fits: the blocks model (infinite
/// strength, the cells' levels) and the link model (finite strength).
fn both_settings() -> [Settings; 2] {
    [
        Settings {
            cell_roles: 3,
            ..settings(Weights::Symmetric, true, f64::INFINITY, 2)
        },
        Settings {
            cell_roles: 3,
            ..settings(Weights::Hdp(1.0), false, 2.0, 1)
        },
    ]
}

/// The same seed gives the same draws with the sweeps counted as without:
/// every kept draw's log joint and each chain's posterior means, bit for
/// bit, in both blockmodel settings and in the mixture.
#[test]
fn counting_the_sweeps_leaves_the_draws_as_they_were() {
    let (reads, gene_of, _, _, _) = planted_blocks(30, 6, 1.0, 161);
    let links = Links::new(&Weighted::new(&reads, None), &gene_of);
    let cfg = config(6, 8, 163);
    for s in both_settings() {
        let model = Model::new(&links, &gene_of, gene_of.len(), s).unwrap();
        let bare = run_mcmc_parallel(&model, &cfg, 2);
        let p = Progress::new(2);
        let watched = progress::run(&model, &cfg, 2, &p).unwrap();
        let joint = |c: &Vec<crate::cluster::unified::Draw>| -> Vec<u64> {
            c.iter().map(|d| d.log_joint.to_bits()).collect()
        };
        for (a, b) in bare.iter().zip(&watched) {
            assert_eq!(joint(a), joint(b));
            let (a, b) = (a.last().unwrap().means(true), b.last().unwrap().means(true));
            assert_eq!(a.cell, b.cell);
            assert_eq!(a.shift, b.shift);
            assert_eq!(a.level, b.level);
        }
        assert_eq!(p.done(0) + p.done(1), 2 * progress::sweeps(&cfg));
    }
    let reads = small_cells(6);
    let weighted = Weighted::new(&reads, None);
    let model = crp::Model::new(&weighted, 2, super::settings(true, 2, None));
    let bare = run_mcmc_parallel(&model, &cfg, 2);
    let watched = progress::run(&model, &cfg, 2, &Progress::new(2)).unwrap();
    for (a, b) in bare.iter().zip(&watched) {
        let draws = |c: &crp::Chain| -> Vec<(Vec<u16>, u64)> {
            c.samples
                .iter()
                .map(|s| (s.z.clone(), s.log_post.to_bits()))
                .collect()
        };
        assert_eq!(draws(a), draws(b));
        assert_eq!((a.splits, a.merges), (b.splits, b.merges));
    }
}

/// The same seed gives the same draws with the log joint traced as
/// without, and the trace is the chain's log joint: at each kept draw's
/// sweep, that draw's, bit for bit, in both blockmodel settings; a short
/// run traced every sweep, a long one every other.
#[test]
fn tracing_the_log_joint_leaves_the_draws_as_they_were() {
    let (reads, gene_of, _, _, _) = planted_blocks(30, 6, 1.0, 161);
    let links = Links::new(&Weighted::new(&reads, None), &gene_of);
    // Draw j is kept after sweep warmup + j * thin + 1: even in the long
    // run, so traced.
    for (cfg, every) in [(config(6, 8, 163), 1), (config(6, 199, 163), 2)] {
        let cfg = McmcConfig { thin: 2, ..cfg };
        let total = progress::sweeps(&cfg);
        for s in both_settings() {
            let model = Model::new(&links, &gene_of, gene_of.len(), s).unwrap();
            let bare = run_mcmc_parallel(&model, &cfg, 2);
            let p = Progress::traced(2);
            let traced = progress::run(&model, &cfg, 2, &p).unwrap();
            for (i, (a, b)) in bare.iter().zip(&traced).enumerate() {
                let joint = |c: &Vec<crate::cluster::unified::Draw>| -> Vec<u64> {
                    c.iter().map(|d| d.log_joint.to_bits()).collect()
                };
                assert_eq!(joint(a), joint(b), "chain {i}");
                assert_eq!(
                    a.last().unwrap().means(true).cell,
                    b.last().unwrap().means(true).cell
                );
                let trace = p.trace(i);
                assert_eq!(
                    trace.iter().map(|t| t.0).collect::<Vec<_>>(),
                    (every..=total).step_by(every).collect::<Vec<_>>(),
                    "traced every {every} of {total} sweeps"
                );
                for (j, d) in b.iter().enumerate() {
                    let at = cfg.warmup + j * cfg.thin + 1;
                    let v = trace.iter().find(|t| t.0 == at).map(|t| t.1);
                    assert_eq!(
                        v.map(f64::to_bits),
                        Some(d.log_joint.to_bits()),
                        "draw {j} at sweep {at}"
                    );
                }
            }
        }
    }
}

/// Each model's fit under `progress`, into `dir`: what it returned.
fn fit_each(
    model: &str,
    progress: &Progress,
    dir: &std::path::Path,
) -> anyhow::Result<crate::cluster::Written> {
    let (data, fit_reads) = planted_site_data(30, 25, 6, 171);
    let weighted = Weighted::new(&fit_reads, None);
    let fitted: Vec<usize> = (0..25).collect();
    let kept: Vec<usize> = (0..30).collect();
    let controls: Vec<usize> = (25..30).collect();
    let args = super::parse_args(&["cluster", "x.zarr", "--model", model, "--chains", "2"]);
    let fit = FitData {
        data: &data,
        fitted: &fitted,
        weighted: &weighted,
        held: None,
        all_cells: &data.cells,
        kept: &kept,
        controls: &controls,
        progress,
    };
    let mut record = crate::run_record::RunRecord::start("cluster", "unused");
    let cfg = config(4, 6, 173);
    match model {
        "crp" => crate::cluster::run::crp::crp_fit(&args, &cfg, dir, &mut record, fit, &fit_reads),
        _ => crate::cluster::run::blockmodel::blockmodel_fit(&args, &cfg, dir, &mut record, fit),
    }
}

/// A cancelled fit returns [`progress::Cancelled`] and writes nothing, in
/// either model; one left to run counts every sweep, notes its R-hat, and
/// writes its tables.
#[test]
fn a_cancelled_fit_writes_nothing() {
    for model in ["blocks", "link", "crp"] {
        let tmp = tempfile::tempdir().unwrap();
        let dir = tmp.path().join("out");
        let p = Progress::new(2);
        p.cancel();
        let out = fit_each(model, &p, &dir);
        assert!(out.err().unwrap().is::<Cancelled>(), "{model}");
        assert!(!dir.exists(), "{model}: nothing written");

        let p = Progress::new(2);
        fit_each(model, &p, &dir).unwrap();
        let total = progress::sweeps(&config(4, 6, 173));
        assert_eq!((p.done(0), p.done(1)), (total, total), "{model}");
        assert!(p.rhat().is_some_and(f64::is_finite), "{model}");
        assert!(!p.cancel(), "{model}: too late to cancel");
        assert!(std::fs::read_dir(&dir).unwrap().count() > 0, "{model}");
    }
}

/// A cancel that comes while the fit reads stops it at the next step,
/// before the next input is read: here before the first, which does not
/// exist.
#[test]
fn a_cancel_while_reading_stops_before_the_next_read() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().join("out");
    let missing = tmp.path().join("missing.zarr");
    let args = super::parse_args(&["cluster", missing.to_str().unwrap()]);
    let mut record = crate::run_record::RunRecord::start("cluster", "unused");
    let p = Progress::new(2);
    p.cancel();
    let out = crate::cluster::cluster(&args, &dir, &mut record, &p);
    assert!(out.err().unwrap().is::<Cancelled>());
    assert!(!dir.exists());
}

/// A progress on stderr draws its bar (none off a terminal) and returns
/// as the engine does once the chains end.
#[test]
fn on_stderr_it_runs_as_the_engine_does() {
    let unused = Progress::new(0);
    let model = Counting::new(usize::MAX, &unused);
    let cfg = config(3, 2, 1);
    let kept = progress::run(&model, &cfg, 2, &Progress::on_stderr(2)).unwrap();
    assert_eq!(kept, vec![3, 3]);
    assert_eq!(
        model.sweeps.load(Ordering::SeqCst),
        2 * progress::sweeps(&cfg)
    );
}

/// The progress never moves the draws: the schedule is the burn-in then
/// `(sweeps - burn-in) / thin` draws every `thin`-th sweep, the engine
/// stopping on the last of them, and the bar's total and burn-in mark are
/// the sweeps the engine runs and where its burn-in ends.
#[test]
fn the_progress_shows_the_engines_schedule() {
    // 3000 sweeps, 1000 of burn-in, thin 2: the engine stops on sweep 2999.
    for (sweeps, burn_in, thin) in [(3000, 1000, 2), (150, 50, 2), (20, 3, 3), (9, 1, 1)] {
        let (s, b, t) = (sweeps.to_string(), burn_in.to_string(), thin.to_string());
        let args = super::parse_args(&[
            "cluster",
            "x.zarr",
            "--sweeps",
            &s,
            "--burn-in",
            &b,
            "--thin",
            &t,
        ]);
        let cfg = crate::cluster::mcmc_config(&args);
        let n_samples = (sweeps - burn_in) / thin;
        assert_eq!(
            (cfg.n_samples, cfg.warmup, cfg.thin, cfg.seed),
            (n_samples, burn_in, thin, args.seed),
            "{cfg:?}"
        );
        let run = burn_in + (n_samples - 1) * thin + 1;
        assert_eq!(progress::sweeps(&cfg), run, "{cfg:?}");

        // The chains, counted, run that many sweeps; the bar says so.
        let shown = Progress::new(1);
        let model = Counting::new(usize::MAX, &shown);
        let kept = progress::run(&model, &cfg, 1, &shown).unwrap();
        assert_eq!(kept, vec![n_samples]);
        assert_eq!(model.sweeps.load(Ordering::SeqCst), run);
        assert_eq!(shown.schedule(), Some((run, burn_in)));
    }
}
