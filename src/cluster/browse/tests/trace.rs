//! The fit's convergence in the browser: the badge in the header from the
//! run record, `t` opening the saved trace full pane with the verdict and
//! its advice under it, and older tables, which kept no trace, saying so.

use super::{drawn, press, screen, tables};
use crate::cluster::browse::{Modules, ModulesView};
use crate::cluster::{args, record_name, write};
use crate::run_record::RunRecord;
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

/// The fixture's tables in `dir` beside a blocks run's record, split
/// R-hat `rhat`, burn-in 20 of 200 sweeps, and with `traced` its trace.
pub(super) fn fitted(dir: &Path, rhat: f64, traced: bool) {
    fitted_with(dir, rhat, traced.then(|| chains([0.0, 0.5])), "2");
}

/// [`fitted`], the trace `traced` (none: no trace) and the run's
/// `--chains`.
fn fitted_with(dir: &Path, rhat: f64, traced: Option<Vec<Vec<(usize, f64)>>>, chains: &str) {
    tables(dir);
    let words = ["faba", "/x.zarr", "--model", "blocks", "--chains", chains];
    let a = crate::cluster::tests::parse_args(&words);
    let mut record = RunRecord::start("cluster", &dir.to_string_lossy())
        .file_name(&record_name("m6a", a.model))
        .options(&a);
    record.set_result("model", &args::recorded(&a));
    record.set_result("blockmodel", &serde_json::json!({"rhat_log_joint": rhat}));
    record.set_result("trace", &serde_json::json!({"sweeps": 200, "burn_in": 20}));
    if let Some(t) = traced {
        write::trace(&t, &dir.join("m6a_trace.parquet")).unwrap();
    }
    record.finish(&Ok(()));
}

/// The modules of [`fitted`].
pub(super) fn load(rhat: f64, traced: bool) -> Modules {
    let tmp = tempfile::tempdir().unwrap();
    fitted(tmp.path(), rhat, traced);
    Modules::load(tmp.path(), "m6a").unwrap()
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
        assert!(
            s.contains("run longer (--sweeps 400) and check the traces"),
            "{s}"
        );
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
    fitted_with(tmp.path(), 1.01, None, "1");
    let m = Modules::load(tmp.path(), "m6a").unwrap();
    let v = m.fit.verdict().unwrap();
    assert_eq!(v.said(), "one chain (R-hat 1.01)");
    assert!(v.advice().contains("--chains"), "{}", v.advice());
}

/// A trace whose chains skip one (no rows for it) reads it back empty;
/// the page leaves it out, the others keeping their names.
#[test]
fn a_chain_with_no_points_is_left_off_the_page() {
    let tmp = tempfile::tempdir().unwrap();
    let [a, b] = <[_; 2]>::try_from(chains([0.0, 0.5])).unwrap();
    fitted_with(tmp.path(), 1.01, Some(vec![a, Vec::new(), b]), "3");
    let m = Modules::load(tmp.path(), "m6a").unwrap();
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
