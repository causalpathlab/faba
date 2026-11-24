//! The verdict on a fit by variational inference: each outcome, its words
//! and the flags its advice names.

use super::Verdict;
use crate::cluster::vi::Unit;

#[test]
fn a_settled_fit_can_be_read() {
    let v = Verdict::vi(true, false, (90, 200), 3, (0.0123, Unit::Read));
    assert!(v.converged());
    assert_eq!(
        v.said(),
        "ELBO settled in 90 epochs; 3 restarts, the worst 0.0123 per read below the best"
    );
    assert_eq!(v.badge(), "ELBO ok");
    assert!(!v.advice().contains("--"), "{}", v.advice());
}

/// One restart has no spread to say.
#[test]
fn one_restart_says_no_spread() {
    let v = Verdict::vi(true, false, (9, 20), 1, (0.0, Unit::Count));
    assert_eq!(v.said(), "ELBO settled in 9 epochs");
}

/// Still rising: twice the epochs asked; stopped early: the same again.
#[test]
fn an_unsettled_fit_names_the_epochs_to_run() {
    let v = Verdict::vi(false, false, (200, 200), 2, (0.5, Unit::Count));
    assert!(!v.converged());
    assert_eq!(v.badge(), "ELBO check");
    assert!(v
        .said()
        .starts_with("ELBO still rising after 200 epochs; 2 restarts"));
    assert!(
        v.said().ends_with("per count below the best"),
        "{}",
        v.said()
    );
    assert_eq!(v.advice(), "run again with --epochs 400");
    let v = Verdict::vi(false, true, (12, 200), 2, (0.5, Unit::Read));
    assert!(v.said().starts_with("stopped early after 12 epochs"));
    assert!(v.advice().ends_with("(--epochs 200)"), "{}", v.advice());
}
