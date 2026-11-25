//! The view's state as JSON: what keys change it says, and set on a view
//! fresh from the same tables it arranges that view the same, `s`
//! drawing the same figure.

use super::*;
use crate::qc::site_plot::PanelArgs;
use crate::tui::live::{sheet_of, Live, Look};
use serde_json::{json, Value};

fn viewer(t: &SiteTable) -> SitePicker<'_> {
    let mut p = picker(t, Some(vec![3; t.len()]), SiteFilterArgs::default_values());
    p.view_only(PanelArgs::default().plan().unwrap(), Sheet::OWN);
    p
}

/// The figure `s` saves now, at the view's own size.
fn saved(p: &SitePicker) -> String {
    p.saved_figure(sheet_of(p, &Look::default()).unwrap())
        .unwrap()
}

#[test]
fn the_state_says_what_the_keys_arranged_and_sets_it_again() {
    let t = table(M6A, 300);
    let mut p = viewer(&t);
    let first = p.view_state();
    assert_eq!(first["focus"], "site_max_pv");
    assert_eq!(first["panel"], "thresholds");
    assert_eq!(first["two_d"], false);
    assert_eq!(first["saves"], json!({ "hist": "pv" }));

    // A knob up (round to the last), a threshold up, its scale stepped, the arms' look
    // turned, then the 2D panels and the one in focus.
    press(&mut p, KeyCode::Up);
    press(&mut p, KeyCode::Char('+'));
    press(&mut p, KeyCode::Char('x'));
    press(&mut p, KeyCode::Char('f'));
    press(&mut p, KeyCode::Char('u'));
    let hist = p.view_state();
    let knob = hist["focus"].as_str().unwrap().to_string();
    assert_ne!(knob, "site_max_pv");
    assert_ne!(hist["thresholds"][&knob], first["thresholds"][&knob]);
    assert_ne!(hist["x_scale"][&knob], first["x_scale"][&knob]);
    assert_eq!(
        (&hist["front"], &hist["p_axis"]),
        (&json!("mut"), &json!("neglog10"))
    );
    assert_eq!(hist["flags"], json!(qc_flags(&p.filter)));
    let hist_svg = saved(&p);

    press(&mut p, KeyCode::Char('v'));
    press(&mut p, KeyCode::Tab);
    press(&mut p, KeyCode::Tab);
    press(&mut p, KeyCode::Right);
    let two_d = p.view_state();
    assert_eq!(two_d["two_d"], true);
    assert_eq!(two_d["panel"], "pairs");
    assert_eq!(two_d["pair_at"], 1);
    assert!(two_d["saves"]["pair"].is_array(), "{}", two_d["saves"]);
    let two_d_svg = saved(&p);

    // Each state set on a view fresh from the table: the same state, and
    // `s` the same figure; the fresh view, before, draws another.
    for (state, svg) in [(&hist, &hist_svg), (&two_d, &two_d_svg)] {
        let mut q = viewer(&t);
        assert_ne!(&saved(&q), svg, "a fresh view draws another figure");
        q.set_view_state(state).unwrap();
        assert_eq!(&q.view_state(), state);
        assert_eq!(&saved(&q), svg);
    }
}

/// A key left out stays; the axis titles, the weight and the 2D scales
/// set as given; a key it cannot set is refused, saying why; a state
/// another view wrote is refused.
#[test]
fn a_state_sets_what_it_names_and_refuses_what_it_cannot() {
    let t = table(M6A, 300);
    let mut p = viewer(&t);
    let before = p.view_state();
    p.set_view_state(&json!({
        "titles": { "x": { "pv": "TITLE1" }, "sites": "TITLE2" },
        "weight": "sites",
        "scales": { "coverage": "linear" },
        "bins": 12,
        "thresholds": { "site_min_coverage": 5, "site_max_pv": null },
    }))
    .unwrap();
    let after = p.view_state();
    assert_eq!(after["titles"]["x"]["pv"], "TITLE1");
    assert_eq!(after["titles"]["sites"], "TITLE2");
    assert_eq!(after["weight"], "sites");
    assert_eq!(after["scales"]["coverage"], "linear");
    assert_eq!(after["bins"], 12);
    assert_eq!(after["thresholds"]["site_min_coverage"], 5.0);
    assert_eq!(after["thresholds"]["site_max_pv"], Value::Null);
    assert_eq!(after["focus"], before["focus"]);
    assert_eq!(after["front"], before["front"]);
    assert!(saved(&p).contains("TITLE1"), "the title drawn");

    for (bad, says) in [
        (
            json!({ "focus": "site_min_nothing" }),
            "no knob site_min_nothing",
        ),
        (json!({ "modality": "atoi" }), "no atoi sites"),
        (
            json!({ "x_scale": { "site_min_log_odds": "log" } }),
            "has no log scale",
        ),
        (json!({ "panel": "pairs" }), "only with two_d"),
        (json!({ "front": "both" }), "front"),
        (json!({ "scales": { "log_odds": "log" } }), "signed"),
    ] {
        let e = format!("{:#}", p.set_view_state(&bad).unwrap_err());
        assert!(e.contains(says), "{bad}: {e}");
    }
    let other = json!({ "view": "cluster_view" });
    let e = super::super::state::set_from_file(&mut p, &other)
        .unwrap_err()
        .to_string();
    assert!(e.contains("cluster_view's, not qc_view's"), "{e}");
}
