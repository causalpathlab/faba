//! The fit's convergence in the browser: the badge in the header from the
//! run record, `t` opening the saved trace full pane with the verdict and
//! its advice under it, tables without a trace saying so, and the MCMC
//! fits of faba before 0.19 refused.

use super::{drawn, press, screen, tables};
use crate::cluster::browse::{Modules, ModulesView};
use crate::cluster::{record_name, write};
use crate::tui::ui::Screen;
use ratatui::crossterm::event::KeyCode;
use std::path::Path;

/// Two restarts' traces, a point every 5 of 200 steps, wobbling about
/// `levels` past a climb over the first 20.
pub(super) fn restarts(levels: [f64; 2]) -> Vec<Vec<(usize, f64)>> {
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

/// The fixture's tables in `dir` beside a blocks fit's record, its best
/// restart `converged` (else still rising) after 30 of 180 epochs, the
/// k-means start 20 steps, and with `traced` its trace.
pub(super) fn fitted(dir: &Path, converged: bool, traced: bool) {
    fitted_with(dir, converged, traced.then(|| restarts([0.0, 0.5])));
}

/// [`fitted`], the trace `traced` (none: no trace).
fn fitted_with(dir: &Path, converged: bool, traced: Option<Vec<Vec<(usize, f64)>>>) {
    tables(dir);
    let n = traced.as_ref().map_or(2, Vec::len);
    let restarts: Vec<_> = (0..n)
        .map(|i| serde_json::json!({"restart": i, "elbo": -9.0}))
        .collect();
    let record = serde_json::json!({
        "job": "cluster", "status": "ok",
        "options": {"chains": n, "epochs": 180},
        "results": {
            "model": {"preset": "blocks", "settings": {
                "cell_roles": 8, "gene_roles": 8, "alpha": 0.5, "gene_alpha": 0.5,
                "weights": "symmetric", "shift_sd": 1.0, "strength": "inf", "cell_sd": 1.0,
            }},
            "trace": {"sweeps": 200, "burn_in": 20, "quantity": "ELBO"},
            "fit": {"method": "vi", "unit": "read", "restarts": restarts, "best": 0,
                "converged": converged, "stopped_early": false, "epochs": 30, "asked": 180,
                "spread_per_unit": 0.01, "tol": 1e-5},
        },
    });
    if let Some(t) = traced {
        write::trace(&t, &dir.join("m6a_trace.parquet")).unwrap();
    }
    std::fs::write(dir.join(record_name("m6a")), record.to_string()).unwrap();
}

/// The modules of [`fitted`].
pub(super) fn load(converged: bool, traced: bool) -> Modules {
    let tmp = tempfile::tempdir().unwrap();
    fitted(tmp.path(), converged, traced);
    Modules::load(tmp.path(), "m6a", None).unwrap()
}

#[test]
fn the_trace_is_read_back_as_written() {
    let m = load(true, true);
    let t = m.fit.trace.as_ref().unwrap();
    assert_eq!(t.restarts, restarts([0.0, 0.5]));
    assert_eq!((m.fit.steps, t.warmup), (200, 20));
    assert!(m.fit.verdict().is_some_and(|v| v.converged()));
}

#[test]
fn the_header_says_whether_the_fit_converged() {
    for (converged, badge) in [(true, "ELBO ok"), (false, "ELBO check")] {
        let m = load(converged, true);
        let mut v = ModulesView::new("test", &m);
        let s = screen(&mut v, 120, 30);
        assert!(s.lines().next().unwrap().contains(badge), "{s}");
    }
    // No record, no badge.
    let m = super::loaded();
    let mut v = ModulesView::new("test", &m);
    assert!(!screen(&mut v, 120, 30).contains("ELBO"));
}

/// `t` opens the trace full pane, a line per restart under the verdict's
/// words and advice; `t` or Esc closes it, and the footer says so.
#[test]
fn t_opens_the_trace_with_the_verdict_and_esc_or_t_closes_it() {
    let m = load(false, true);
    let mut v = ModulesView::new("test", &m);
    assert!(
        screen(&mut v, 120, 30).contains("t trace"),
        "listed in the footer"
    );
    for close in [KeyCode::Esc, KeyCode::Char('t')] {
        press(&mut v, KeyCode::Char('t'));
        let s = screen(&mut v, 120, 30);
        assert!(s.contains("ELBO still rising after 30 epochs"), "{s}");
        assert!(s.contains("run again with --epochs 360"), "{s}");
        assert!(s.contains("close") && !s.contains("gene modules ·"), "{s}");
        let buf = drawn(&mut v, 120, 30);
        let braille = (0..30)
            .flat_map(|y| (0..120).map(move |x| (x, y)))
            .filter(|&at| {
                ('\u{2801}'..='\u{28ff}').contains(&buf[at].symbol().chars().next().unwrap())
            })
            .count();
        assert!(braille > 20, "the restarts drawn: {braille} cells");
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

/// Tables without a trace: `t` says so in one line and opens nothing; the
/// badge still comes from the record.
#[test]
fn tables_without_a_trace_say_why() {
    let m = load(true, false);
    let mut v = ModulesView::new("test", &m);
    press(&mut v, KeyCode::Char('t'));
    let s = screen(&mut v, 120, 30);
    assert!(s.contains("gene modules ·"), "nothing opened: {s}");
    assert!(s.contains("no m6a_trace.parquet"), "{s}");
    assert!(s.lines().next().unwrap().contains("ELBO ok"), "{s}");
}

/// Tables without their run record have neither a trace to draw (its
/// k-means steps are the record's) nor a verdict.
#[test]
fn without_a_record_there_is_no_trace_or_verdict() {
    let tmp = tempfile::tempdir().unwrap();
    fitted(tmp.path(), true, true);
    std::fs::remove_file(tmp.path().join(record_name("m6a"))).unwrap();
    let m = Modules::load(tmp.path(), "m6a", None).unwrap();
    assert!(m.fit.verdict().is_none() && m.fit.trace.is_err());
}

/// A trace whose restarts skip one (no rows for it) reads it back empty;
/// the page leaves it out, the others keeping their names.
#[test]
fn a_restart_with_no_points_is_left_off_the_page() {
    let tmp = tempfile::tempdir().unwrap();
    let [a, b] = <[_; 2]>::try_from(restarts([0.0, 0.5])).unwrap();
    fitted_with(tmp.path(), true, Some(vec![a, Vec::new(), b]));
    let m = Modules::load(tmp.path(), "m6a", None).unwrap();
    assert_eq!(m.fit.trace.as_ref().unwrap().restarts[1], Vec::new());
    let data = serde_json::to_value(crate::cluster::browse::report::data(
        &m,
        &Default::default(),
    ))
    .unwrap();
    let drawn = data["trace"]["restarts"].as_array().unwrap();
    let names: Vec<_> = drawn.iter().map(|c| c["name"].as_str().unwrap()).collect();
    assert_eq!(names, ["restart 0", "restart 2"]);
    assert!(drawn
        .iter()
        .all(|c| !c["points"].as_array().unwrap().is_empty()));
}

/// The record of an MCMC fit of faba before 0.19 (no fit by variational
/// inference in it, a split R-hat instead) is refused, saying to fit again.
#[test]
fn an_mcmc_fit_is_refused() {
    let tmp = tempfile::tempdir().unwrap();
    tables(tmp.path());
    let record = serde_json::json!({
        "job": "cluster", "status": "ok",
        "results": {"model": {}, "blockmodel": {"rhat_log_joint": 1.02},
            "trace": {"sweeps": 200, "burn_in": 20}},
    });
    let at = tmp.path().join(record_name("m6a"));
    std::fs::write(&at, record.to_string()).unwrap();
    let err = Modules::load(tmp.path(), "m6a", None).unwrap_err();
    assert!(
        format!("{err:#}").contains("MCMC fit of faba before 0.19"),
        "{err:#}"
    );
    let err = crate::cluster::browse::fit_record(&[at.to_string_lossy().into()]).unwrap_err();
    assert!(
        format!("{err:#}").contains("fit again with `faba cluster`"),
        "{err:#}"
    );
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
    write::trace(&restarts([-2.0, -2.1]), &dir.join("m6a_trace.parquet")).unwrap();
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
