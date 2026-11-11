//! What the fits by variational inference (the count fit's, the site
//! blockmodel's) run alike: their restarts, each from a k-means start of its
//! own, the best by ELBO kept; and their verdict and record.

use super::super::args::ClusterArgs;
use super::super::progress::{self, advice};
use super::super::vi::{Hooks, Options, Ran};
use crate::run_record::RunRecord;
use rand::rngs::SmallRng;
use rand::SeedableRng;

/// The variational fit's options from `args`.
pub(in crate::cluster) fn options(args: &ClusterArgs) -> Options {
    Options {
        svi_epochs: args.svi_epochs,
        batch_size: args.batch_size,
        epochs: args.epochs,
        tol: args.tol,
        ..Options::default()
    }
}

/// How a restart ended, for the record and the verdict: its ELBO, that per
/// unit (count or read), and how it ran.
#[derive(Debug, Clone, Copy)]
pub(in crate::cluster) struct End {
    pub elbo: f64,
    pub per_unit: f64,
    pub ran: Ran,
}

/// The best restart's fit and index (none when no restart finished an
/// epoch), and each restart's end (none for one stopped before an epoch
/// ended or whose start failed).
pub(in crate::cluster) type Restarts<F> = (Option<(usize, F)>, Vec<Option<End>>);

/// `restarts` restarts: the starts side by side (each on its own thread,
/// told its Lloyd iterations and whether to stop), then the fits one after
/// another, each using every thread. Restart `i`'s random numbers are
/// seeded `seed + i`; its start's `warmup` steps, then its epochs (of
/// `epochs` at most), are counted on `progress` as chain `i` and its ELBO
/// traced there. Only the best fit by ELBO is kept. Refused when no restart
/// could start and none was stopped.
pub(in crate::cluster) fn restarts<S: Send, F>(
    progress: &progress::Progress,
    (seed, restarts): (u64, usize),
    (warmup, epochs): (usize, usize),
    start: impl Fn(&mut SmallRng, &(dyn Fn(usize) + Sync), &(dyn Fn() -> bool + Sync)) -> anyhow::Result<S>
        + Sync,
    fit: impl Fn(S, &mut SmallRng, &Hooks) -> Option<F>,
    end: impl Fn(&F) -> End,
) -> anyhow::Result<Restarts<F>> {
    let stop = || progress.stopping();
    let cancelled = || progress.cancelled();
    let starts: Vec<anyhow::Result<(S, SmallRng)>> = std::thread::scope(|s| {
        let (start, stop) = (&start, &stop);
        let handles: Vec<_> = (0..restarts)
            .map(|i| {
                s.spawn(move || {
                    let mut rng = SmallRng::seed_from_u64(seed.wrapping_add(i as u64));
                    progress.reach(i, 0);
                    let tick = |n: usize| progress.reach(i, n);
                    start(&mut rng, &tick, stop).map(|st| (st, rng))
                })
            })
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("a k-means start does not panic"))
            .collect()
    });
    // What the starts alone needed goes with them.
    drop(start);
    // Every start done (its iterations may have ended early): a restart
    // waits at its first epoch for its turn.
    (0..restarts).for_each(|i| progress.reach(i, warmup));
    let mut best: Option<(usize, F, End)> = None;
    let mut ends = Vec::with_capacity(restarts);
    let mut failed = None;
    for (i, start) in starts.into_iter().enumerate() {
        let (start, mut rng) = match start {
            Ok(start) => start,
            Err(e) => {
                log::warn!("restart {i}: {e:#}");
                failed.get_or_insert(e);
                ends.push(None);
                continue;
            }
        };
        if stop() {
            ends.push(None);
            continue;
        }
        let epoch = |e: usize, f: f64| {
            progress.reach(i, warmup + e);
            progress.note(i, warmup + e, f);
        };
        let hooks = Hooks {
            epoch: &epoch,
            stop: &stop,
            cancelled: &cancelled,
        };
        let Some(fitted) = fit(start, &mut rng, &hooks) else {
            ends.push(None);
            continue;
        };
        let this = end(&fitted);
        // A restart that settled before its last epoch has nothing left.
        if this.ran.converged {
            progress.reach(i, warmup + epochs);
        }
        ends.push(Some(this));
        if best.as_ref().is_none_or(|(_, _, b)| this.elbo > b.elbo) {
            best = Some((i, fitted, this));
        }
    }
    if let (None, Some(e), false) = (&best, failed, stop()) {
        return Err(e);
    }
    Ok((best.map(|(i, f, _)| (i, f)), ends))
}

/// The restarts' outcome a driver writes from: the best restart's index,
/// fit and end, each restart's end, and whether the run was stopped early
/// (as it stood when the fits ended: a stop during the tables' passes
/// after does not make the fit one stopped early).
pub(in crate::cluster) struct Best<F> {
    pub index: usize,
    pub fit: F,
    pub end: End,
    pub ends: Vec<Option<End>>,
    pub stopped: bool,
}

/// The restarts of `work` (see [`restarts`]), their `warmup` steps of
/// k-means then `epochs` at most counted on `progress` (drawn on its bar);
/// the best of them. [`progress::Cancelled`] when the run was, or was
/// stopped before any restart finished an epoch (nothing to write).
pub(in crate::cluster) fn best<F: Send>(
    progress: &progress::Progress,
    (warmup, epochs): (usize, usize),
    work: impl FnOnce() -> anyhow::Result<Restarts<F>> + Send,
) -> anyhow::Result<Best<F>> {
    progress.begin_vi(warmup + epochs, warmup);
    let (best, ends) = progress::barred(progress, work)?;
    if progress.cancelled() {
        return Err(progress::Cancelled.into());
    }
    let Some((index, fit)) = best else {
        log::warn!("stopped before any restart finished an epoch: nothing to write");
        return Err(progress::Cancelled.into());
    };
    let end = ends[index].expect("the best restart ended");
    let stopped = progress.stopped_early().is_some();
    Ok(Best {
        index,
        fit,
        end,
        ends,
        stopped,
    })
}

/// Roles a fit by variational inference keeps but gives (nearly) nothing,
/// from its blocks' units (`A x B`, `bb` gene roles): each side's, from 1.
pub(in crate::cluster) fn empty_roles(units: &[f64], bb: usize) -> serde_json::Value {
    let empty = |totals: Vec<f64>| -> Vec<usize> {
        (totals.iter().enumerate())
            .filter(|(_, &u)| u < 1.0)
            .map(|(k, _)| k + 1)
            .collect()
    };
    let cell = units.chunks_exact(bb).map(|row| row.iter().sum()).collect();
    let gene = super::super::count::column_totals(units, bb);
    serde_json::json!({"cell": empty(cell), "gene": empty(gene)})
}

/// The verdict on the restarts `best` came of, out of `asked` epochs at
/// most, noted on `progress`: the run then commits to writing (refused as
/// cancelled when it was). What the record
/// keeps of them, under `results.fit` (with `empty_roles`, see
/// [`empty_roles`]), and the trace's schedule. The ELBO per `unit` (count
/// or read) compares the restarts.
pub(in crate::cluster) fn conclude<F>(
    progress: &progress::Progress,
    record: &mut RunRecord,
    unit: &'static str,
    best: &Best<F>,
    (warmup, asked, tol): (usize, usize, f64),
    empty_roles: serde_json::Value,
) -> anyhow::Result<()> {
    let (ends, kept, stopped) = (&best.ends, best.end, best.stopped);
    let best = best.index;
    let per: Vec<f64> = ends.iter().flatten().map(|e| e.per_unit).collect();
    let spread = per.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        - per.iter().copied().fold(f64::INFINITY, f64::min);
    if stopped {
        let epochs: Vec<Option<usize>> = ends.iter().map(|e| e.map(|e| e.ran.epochs)).collect();
        log::warn!(
            "stopped early: the tables are from the epochs run so far ({epochs:?} by restart)"
        );
        // As 0.18 recorded it, for whoever reads it there.
        record.set_result("stopped_early", &serde_json::json!({"epochs": epochs}));
    }
    let verdict = advice::Verdict::vi(
        kept.ran.converged,
        stopped,
        (kept.ran.epochs, asked),
        per.len(),
        (spread, unit),
    );
    progress.optimised(verdict.clone())?;
    log::info!(
        "each restart's ELBO per {unit} {:?}; tables from restart {best} ({} epochs: {})",
        per.iter()
            .map(|x| (x * 1e4).round() / 1e4)
            .collect::<Vec<_>>(),
        kept.ran.epochs,
        verdict.said(),
    );
    record.set_result(
        "trace",
        &serde_json::json!({
            "sweeps": warmup + asked,
            "burn_in": warmup,
            "quantity": progress::ELBO,
        }),
    );
    // (A restart that never ran says so: readers count those with an ELBO.)
    let restarts: Vec<serde_json::Value> = (ends.iter().enumerate())
        .map(|(i, e)| match e {
            Some(e) => serde_json::json!({
                "restart": i,
                "elbo": e.elbo,
                "elbo_per_unit": e.per_unit,
                "epochs": e.ran.epochs,
                "svi_epochs": e.ran.svi_epochs,
                "converged": e.ran.converged,
            }),
            None => serde_json::json!({"restart": i, "ran": false}),
        })
        .collect();
    record.set_result(
        "fit",
        &serde_json::json!({
            "method": "vi",
            "unit": unit,
            "restarts": restarts,
            "best": best,
            "converged": kept.ran.converged,
            "stopped_early": stopped,
            "epochs": kept.ran.epochs,
            "asked": asked,
            "spread_per_unit": spread,
            "tol": tol,
            "empty_roles": empty_roles,
        }),
    );
    record.set_result("converged", &verdict.converged());
    Ok(())
}
