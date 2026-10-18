//! The verdict on synthetic traces: each outcome, its words and the flags
//! its advice names.

use super::{Status, Verdict};

/// A chain's trace, a point every sweep to 100, wobbling about `level`
/// past a climb over its first 20 sweeps.
fn chain(level: f64) -> Vec<(usize, f64)> {
    (1..=100)
        .map(|n| {
            let climb = if n <= 20 { -50.0 } else { 0.0 };
            (n, level + climb + (n as f64 * 1.7).sin())
        })
        .collect()
}

/// The verdict on chains at `levels`, burn-in 20 of 1000 sweeps, R-hat `rhat`.
fn verdict(levels: &[f64], rhat: f64) -> Verdict {
    let traces: Vec<_> = levels.iter().map(|&l| chain(l)).collect();
    Verdict::of(&traces, 20, 1000, rhat)
}

#[test]
fn chains_that_agree_have_converged() {
    let v = verdict(&[0.0, 0.2, -0.1], 1.02);
    assert_eq!(v.status, Status::Agree);
    assert!(v.converged());
    assert_eq!(v.said(), "chains agree (R-hat 1.02)");
    assert_eq!(v.badge(), "R-hat 1.02 ok");
    assert!(!v.advice().contains("--"), "{}", v.advice());
}

#[test]
fn a_chain_far_below_lags_and_more_sweeps_are_advised() {
    let v = verdict(&[0.0, 0.2, -30.0], 1.04);
    assert_eq!(v.status, Status::Lags(vec![2]));
    assert!(!v.converged());
    assert_eq!(v.said(), "chain 2 lags (R-hat 1.04)");
    assert_eq!(v.advice(), "run again with --sweeps 2000 or more --chains");
    assert_eq!(v.badge(), "R-hat 1.04 check");
    let v = verdict(&[0.0, -30.0, -30.0], 1.04);
    assert_eq!(v.said(), "chains 1, 2 lag (R-hat 1.04)");
}

/// Past 1.05 the chains disagree, whatever their levels; at 1.05 they agree.
#[test]
fn an_r_hat_past_the_bound_says_the_chains_disagree() {
    let v = verdict(&[0.0, 0.2, -0.1], 1.31);
    assert_eq!(v.status, Status::Disagree(Vec::new()));
    assert!(!v.converged());
    assert_eq!(v.said(), "R-hat 1.31: chains disagree");
    assert_eq!(
        v.advice(),
        "run longer (--sweeps 2000) and check the traces"
    );
    assert_eq!(v.badge(), "R-hat 1.31 check");
    let v = verdict(&[0.0, -30.0], 1.31);
    assert_eq!(v.said(), "R-hat 1.31: chains disagree, chain 1 lags");
    assert!(verdict(&[0.0, 0.2], 1.05).converged());
    assert!(!verdict(&[0.0, 0.2], 1.06).converged());
}

#[test]
fn one_chain_has_nothing_to_hold_it_against() {
    let v = verdict(&[0.0], 1.01);
    assert_eq!(v.status, Status::OneChain);
    assert!(v.converged());
    assert_eq!(v.said(), "one chain (R-hat 1.01)");
    assert!(v.advice().contains("--chains"), "{}", v.advice());
}

/// One chain whose halves disagree: its halves are what disagree, and
/// other chains are what would show it.
#[test]
fn one_chain_past_the_bound_says_its_halves_disagree() {
    let v = verdict(&[0.0], 1.3);
    assert_eq!(v.status, Status::Disagree(Vec::new()));
    assert!(!v.converged());
    assert_eq!(v.said(), "R-hat 1.30: one chain, its halves disagree");
    assert_eq!(
        v.advice(),
        "run longer (--sweeps 2000) with --chains 4 and check the trace"
    );
}

/// No R-hat (too few draws) is not converged: run longer.
#[test]
fn no_r_hat_asks_for_more_sweeps() {
    let v = verdict(&[0.0, 0.2], f64::NAN);
    assert_eq!(v.status, Status::Unknown);
    assert!(!v.converged());
    assert_eq!(v.said(), "R-hat not known");
    assert!(v.advice().contains("--sweeps 2000"), "{}", v.advice());
    assert_eq!(v.badge(), "R-hat - check");
}

/// Without a trace (an older fit) the R-hat alone decides, the chains
/// counted from their empty traces: one is one chain.
#[test]
fn without_a_trace_the_r_hat_alone_decides() {
    let none = vec![Vec::new(); 2];
    assert_eq!(Verdict::of(&none, 20, 1000, 1.02).status, Status::Agree);
    assert_eq!(
        Verdict::of(&none, 20, 1000, 1.5).status,
        Status::Disagree(Vec::new())
    );
    assert_eq!(
        Verdict::of(&none[..1], 20, 1000, 1.02).status,
        Status::OneChain
    );
}
