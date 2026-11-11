//! The progress a fit's restarts count on: each restart's steps and trace
//! its own, the verdict once the fit ends, and a cancel or a stop asked
//! for too late noted.

use super::*;

/// The fit's verdict, as it ends.
fn verdict() -> advice::Verdict {
    advice::Verdict::vi(true, false, (3, 10), 2, (0.0, "read"))
}

/// Each restart counts on its own bar: a step reached is never undone by
/// an earlier one, and the restarts' steps add up.
#[test]
fn each_restart_counts_on_its_own_bar() {
    let p = Progress::new(2);
    assert_eq!(p.running(0), None, "not started");
    p.reach(0, 5);
    p.reach(0, 3);
    p.reach(1, 2);
    assert_eq!((p.done(0), p.done(1), p.sum()), (5, 2, 7));
    assert!(p.running(0).is_some());
}

/// A restart notes its ELBO at each epoch, a point at the same step
/// replacing the one before (an estimate made exact).
#[test]
fn a_restart_notes_its_elbo_by_step() {
    let p = Progress::new(1);
    p.note(0, 3, -2.0);
    p.note(0, 4, -1.5);
    p.note(0, 4, -1.25);
    assert_eq!(p.trace(0), vec![(3, -2.0), (4, -1.25)]);
}

/// The verdict is none until the fit ends, then the one it ended with;
/// the schedule is the k-means start's steps then the epochs.
#[test]
fn the_verdict_comes_when_the_fit_ends() {
    let p = Progress::new(2);
    assert_eq!((p.verdict(), p.schedule()), (None, None));
    p.begin_vi(40, 30);
    assert_eq!(p.schedule(), Some((40, 30)));
    p.optimised(verdict()).unwrap();
    // (Its R-hat is NaN, never equal: compared by what it says.)
    assert_eq!(p.verdict().map(|v| v.said()), Some(verdict().said()));
}

/// A cancel or a stop asked for once the fit is writing is too late for
/// it, and is noted; one in time is not.
#[test]
fn a_cancel_too_late_for_the_fit_is_noted() {
    let p = Progress::new(1);
    p.optimised(verdict()).unwrap();
    assert!(!p.asked_late());
    assert!(!p.cancel() && p.asked_late() && !p.cancelled());
    let p = Progress::new(1);
    p.optimised(verdict()).unwrap();
    p.stop_early();
    assert!(p.asked_late() && p.stopped_early().is_none());
    let p = Progress::new(1);
    assert!(p.cancel() && !p.asked_late());
}
