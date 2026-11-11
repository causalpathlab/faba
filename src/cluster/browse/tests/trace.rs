//! The fit's convergence in the browser: the badge in the header from the
//! run record, `t` opening the saved trace full pane with the verdict and
//! its advice under it, and older tables, which kept no trace, saying so.

use super::{drawn, press, screen, tables};
use crate::cluster::browse::{Modules, ModulesView};
use crate::cluster::{record_name, write};
use data_beans::interactive::ui::Screen;
use ratatui::crossterm::event::KeyCode;
use std::path::Path;

/// Two chains' traces, a point every 5 of 200 sweeps, wobbling about
/// `levels` past a climb over the first 20.
pub(super) fn chains(levels: [f64; 2]) -> Vec<Vec<(usize, f64)>> {
    levels
        .iter()
        .map(|&l| {
            (1..=40)
                .map(|i| {
                    let n = 5 * i;
                    let climb = if n <= 20 { -40.0 } else { 0.0 };
                    (n, l + climb + (n as f64).sin())
                })
                .collect()
        })
        .collect()
}

/// The fixture's tables in `dir` beside a Gibbs-sampled blocks run's
/// record (faba before 0.19), split R-hat `rhat`, burn-in 20 of 200
/// sweeps, and with `traced` its trace.
pub(super) fn fitted(dir: &Path, rhat: f64, traced: bool) {
    fitted_with(dir, rhat, traced.then(|| chains([0.0, 0.5])), 2, &[]);
}

/// [`fitted`], the trace `traced` (none: no trace), the run's `--chains`
/// and what its pooling found `unsettled`, the record written as faba
/// before 0.19 wrote it: no quantity in its trace (the log joint).
fn fitted_with(
    dir: &Path,
    rhat: f64,
    traced: Option<Vec<Vec<(usize, f64)>>>,
    chains: usize,
    unsettled: &[&str],
) {
    tables(dir);
    let record = serde_json::json!({
        "job": "cluster", "status": "ok",
        "options": {"model": "blocks", "chains": chains, "sweeps": 200, "burn_in": 20},
        "results": {
            "model": {"preset": "blocks", "settings": {
                "cell_roles": 8, "gene_roles": 8, "alpha": 0.5, "gene_alpha": 0.5,
                "weights": {"hdp": 1.0}, "shift_sd": 1.0, "strength": "inf", "cell_sd": 1.0,
            }},
            "blockmodel": {"rhat_log_joint": rhat},
            "trace": {"sweeps": 200, "burn_in": 20},
            "pooling": {"unsettled": unsettled},
        },
    });
    if let Some(t) = traced {
        write::trace(&t, &dir.join("m6a_trace.parquet")).unwrap();
    }
    std::fs::write(dir.join(record_name("m6a")), record.to_string()).unwrap();
}

/// The modules of [`fitted`].
pub(super) fn load(rhat: f64, traced: bool) -> Modules {
    let tmp = tempfile::tempdir().unwrap();
    fitted(tmp.path(), rhat, traced);
    Modules::load(tmp.path(), "m6a", None).unwrap()
}

#[test]
fn the_trace_is_read_back_as_written() {
    let m = load(1.02, true);
    let t = m.fit.trace.as_ref().unwrap();
    assert_eq!(t.chains, chains([0.0, 0.5]));
    assert_eq!((m.fit.sweeps, t.burn_in), (200, 20));
    assert_eq!(m.fit.verdict().map(|v| v.rhat), Some(1.02));
}

#[test]
fn the_header_says_whether_the_fit_converged() {
    for (rhat, badge) in [(1.02, "R-hat 1.02 ok"), (1.31, "R-hat 1.31 check")] {
        let m = load(rhat, true);
        let mut v = ModulesView::new("test", &m);
        let s = screen(&mut v, 120, 30);
        assert!(s.lines().next().unwrap().contains(badge), "{s}");
    }
    // No record, no badge.
    let m = super::loaded();
    let mut v = ModulesView::new("test", &m);
    assert!(!screen(&mut v, 120, 30).contains("R-hat"));
}

/// `t` opens the trace full pane, a line per chain under the verdict's
/// words and advice; `t` or Esc closes it, and the footer says so.
#[test]
fn t_opens_the_trace_with_the_verdict_and_esc_or_t_closes_it() {
    let m = load(1.31, true);
    let mut v = ModulesView::new("test", &m);
    assert!(
        screen(&mut v, 120, 30).contains("t trace"),
        "listed in the footer"
    );
    for close in [KeyCode::Esc, KeyCode::Char('t')] {
        press(&mut v, KeyCode::Char('t'));
        let s = screen(&mut v, 120, 30);
        assert!(s.contains("R-hat 1.31: chains disagree"), "{s}");
        assert!(s.contains("fit again with `faba cluster`"), "{s}");
        assert!(s.contains("close") && !s.contains("gene modules ·"), "{s}");
        let buf = drawn(&mut v, 120, 30);
        let braille = (0..30)
            .flat_map(|y| (0..120).map(move |x| (x, y)))
            .filter(|&at| {
                ('\u{2801}'..='\u{28ff}').contains(&buf[at].symbol().chars().next().unwrap())
            })
            .count();
        assert!(braille > 20, "the chains drawn: {braille} cells");
        press(&mut v, close);
        assert!(
            screen(&mut v, 120, 30).contains("gene modules ·"),
            "{close:?} closes it"
        );
        assert!(!v.done());
    }
    press(&mut v, KeyCode::Char('t'));
    press(&mut v, KeyCode::Char('q'));
    assert!(v.done(), "q quits from the trace too");
}

/// Older tables kept no trace: `t` says so in one line and opens nothing;
/// the badge still comes from the record's R-hat.
#[test]
fn older_tables_without_a_trace_say_why() {
    let m = load(1.02, false);
    let mut v = ModulesView::new("test", &m);
    press(&mut v, KeyCode::Char('t'));
    let s = screen(&mut v, 120, 30);
    assert!(s.contains("gene modules ·"), "nothing opened: {s}");
    assert!(s.contains("no m6a_trace.parquet"), "{s}");
    assert!(s.lines().next().unwrap().contains("R-hat 1.02 ok"), "{s}");
}

/// An older single chain's fit, with no trace: one chain, by the record's
/// `--chains`, with the advice to run more.
#[test]
fn an_older_single_chain_is_one_chain() {
    let tmp = tempfile::tempdir().unwrap();
    fitted_with(tmp.path(), 1.01, None, 1, &[]);
    let m = Modules::load(tmp.path(), "m6a", None).unwrap();
    let v = m.fit.verdict().unwrap();
    assert_eq!(v.said(), "one chain (R-hat 1.01)");
    assert!(v.advice().contains("compare restarts"), "{}", v.advice());
}

/// A trace whose chains skip one (no rows for it) reads it back empty;
/// the page leaves it out, the others keeping their names.
#[test]
fn a_chain_with_no_points_is_left_off_the_page() {
    let tmp = tempfile::tempdir().unwrap();
    let [a, b] = <[_; 2]>::try_from(chains([0.0, 0.5])).unwrap();
    fitted_with(tmp.path(), 1.01, Some(vec![a, Vec::new(), b]), 3, &[]);
    let m = Modules::load(tmp.path(), "m6a", None).unwrap();
    assert_eq!(m.fit.trace.as_ref().unwrap().chains[1], Vec::new());
    let data = serde_json::to_value(crate::cluster::browse::report::data(
        &m,
        &Default::default(),
    ))
    .unwrap();
    let drawn = data["trace"]["chains"].as_array().unwrap();
    let names: Vec<_> = drawn.iter().map(|c| c["name"].as_str().unwrap()).collect();
    assert_eq!(names, ["chain 0", "chain 2"]);
    assert!(drawn
        .iter()
        .all(|c| !c["points"].as_array().unwrap().is_empty()));
}

/// Chains that agree on their traces, but whose pooling found quantities
/// short of convergence, are not ok: the header says to check, the trace
/// names them.
#[test]
fn what_the_pooling_found_unsettled_is_checked() {
    let tmp = tempfile::tempdir().unwrap();
    fitted_with(
        tmp.path(),
        1.02,
        Some(chains([0.0, 0.5])),
        2,
        &["shift (0, 1)"],
    );
    let m = Modules::load(tmp.path(), "m6a", None).unwrap();
    let mut v = ModulesView::new("test", &m);
    let s = screen(&mut v, 120, 30);
    assert!(
        s.lines().next().unwrap().contains("R-hat 1.02 check"),
        "{s}"
    );
    press(&mut v, KeyCode::Char('t'));
    let s = screen(&mut v, 120, 30);
    assert!(s.contains("shift (0, 1)"), "{s}");
}

/// A fit by variational inference's record: its verdict under `fit` with
/// its unit, or as faba 0.18 kept it, under `count` or `blocks` with the
/// unit in the field's name; its trace the restarts' ELBO.
fn vi_fitted(dir: &Path, results: serde_json::Value) {
    tables(dir);
    let mut record = serde_json::json!({
        "job": "cluster", "status": "ok", "options": {"chains": 2},
        "results": {"trace": {"sweeps": 60, "burn_in": 20, "quantity": "ELBO"}},
    });
    if let (Some(r), serde_json::Value::Object(more)) = (record["results"].as_object_mut(), results)
    {
        r.extend(more);
    }
    write::trace(&chains([-2.0, -2.1]), &dir.join("m6a_trace.parquet")).unwrap();
    std::fs::write(dir.join(record_name("m6a")), record.to_string()).unwrap();
}

/// The verdict of a fit by variational inference reads from its record in
/// either shape, the badge and the pane saying the ELBO, the restarts named
/// as such.
#[test]
fn a_vi_fit_reads_its_verdict_in_either_shape() {
    let restarts = serde_json::json!([{"restart": 0, "elbo": -9.0}, {"restart": 1, "ran": false}]);
    let verdict = |unit: Option<&str>| {
        let mut v = serde_json::json!({"restarts": restarts, "best": 0, "converged": true,
            "stopped_early": false, "epochs": 12, "asked": 40, "tol": 1e-5});
        match unit {
            Some(u) => {
                v["unit"] = u.into();
                v["spread_per_unit"] = 0.0.into();
            }
            None => v["spread_per_read"] = 0.0.into(),
        }
        v
    };
    for results in [
        serde_json::json!({"fit": verdict(Some("read"))}),
        serde_json::json!({"blocks": {"vi": verdict(None)}}),
        serde_json::json!({"count": {"vi": verdict(Some("count"))}}),
    ] {
        let tmp = tempfile::tempdir().unwrap();
        vi_fitted(tmp.path(), results.clone());
        let m = Modules::load(tmp.path(), "m6a", None).unwrap();
        let v = m.fit.verdict().unwrap_or_else(|| panic!("{results}"));
        assert!(
            v.converged() && v.said().contains("ELBO settled in 12 epochs"),
            "{results}"
        );
        let mut view = ModulesView::new("test", &m);
        assert!(screen(&mut view, 120, 30)
            .lines()
            .next()
            .unwrap()
            .contains("ELBO ok"));
        press(&mut view, KeyCode::Char('t'));
        let s = screen(&mut view, 120, 30);
        assert!(s.contains("restart 0") && !s.contains("chain 0"), "{s}");
    }
}
