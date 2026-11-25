//! The view read only, as `faba qc-view` opens it: nothing to apply, no
//! flags to print, `q` leaving at once and Esc never, the keys it lists,
//! and `s` saving the panel on screen as a paper panel.

use super::*;
use crate::figure::style::{Palette, StyleArgs};
use crate::qc::site_plot::PanelArgs;

/// A read-only picker over one table, its panels as `panels` asks, laid
/// out for `sheet`.
fn viewer<'a>(
    t: &'a SiteTable,
    n_cells: Option<Vec<usize>>,
    panels: &PanelArgs,
    sheet: Sheet,
) -> SitePicker<'a> {
    let mut p = picker(t, n_cells, SiteFilterArgs::default_values());
    p.viewing = Some((panels.plan().unwrap(), sheet));
    p
}

fn sheet(w: f64, h: f64) -> Sheet {
    Sheet {
        page: (Some(w), Some(h)),
        look: StyleArgs {
            font_size: Some(6.0),
            panel_only: true,
            palette: Palette::Gray,
        },
    }
}

/// Ctrl+R is not taken, so no confirmation opens; `p` prints nothing;
/// Esc never leaves; `q` leaves at once, even with a threshold moved.
#[test]
fn nothing_is_applied_and_q_leaves_at_once() {
    let t = table(M6A, 300);
    let mut p = viewer(&t, None, &PanelArgs::default(), Sheet::OWN);
    assert!(!p.takes_apply());
    apply_key(&mut p);
    assert!(matches!(p.mode, Mode::Browse), "Ctrl+R opens nothing");
    press(&mut p, KeyCode::Char('+'));
    assert_ne!(qc_flags(&p.filter), qc_flags(&p.initial), "thresholds move");
    press(&mut p, KeyCode::Char('p'));
    assert!(
        matches!(p.mode, Mode::Browse) && !p.done(),
        "no flags to print"
    );
    press(&mut p, KeyCode::Esc);
    assert!(!p.done(), "Esc never leaves a top-level view");
    press(&mut p, KeyCode::Enter);
    assert!(!p.enter_hint, "Enter names no apply key");
    press(&mut p, KeyCode::Char('q'));
    assert!(
        matches!(p.decision, Some(Picked::Cancelled)),
        "q leaves at once"
    );
}

/// The footer, the `?` pop-up and the tabs name no apply key, button or
/// flags to print; `q` is quit, and the rest is as in `faba qc`.
#[test]
fn the_keys_listed_are_the_ones_bound() {
    let t = table(M6A, 300);
    let mut p = viewer(&t, None, &PanelArgs::default(), Sheet::OWN);
    let footer = p.footer_line().to_string();
    assert!(footer.contains("q quit"), "{footer}");
    for gone in ["Ctrl+R", "print flags", "cancel"] {
        assert!(!footer.contains(gone), "{gone}: {footer}");
    }
    let keys = crate::tui::help::shown(&mut p, 120, 50).expect("? opens the keys");
    for key in ["Tab/Shift+Tab", "x/y", "v", "s", "Ctrl+E", "/", "?", "quit"] {
        assert!(keys.contains(key), "{key}:\n{keys}");
    }
    for gone in ["Ctrl+R", "print the flags", "asks when changed"] {
        assert!(!keys.contains(gone), "{gone}:\n{keys}");
    }
    let screen = crate::tui::screen_text(&mut p, 150, 44);
    assert!(
        screen.contains("qc-view") && !screen.contains("apply"),
        "{screen}"
    );
}

/// Ctrl-C leaves at once, as `q` does.
#[test]
fn an_interrupt_leaves_at_once() {
    let t = table(M6A, 300);
    let mut p = viewer(&t, None, &PanelArgs::default(), Sheet::OWN);
    press(&mut p, KeyCode::Char('+'));
    p.interrupt();
    assert!(matches!(p.decision, Some(Picked::Cancelled)));
}

/// `s` saves the panel on screen: the focused knob's histogram, or in 2D
/// the focused pair with the marginals asked, each the size asked.
#[test]
fn s_saves_the_panel_on_screen_as_a_paper_panel() {
    let t = table(M6A, 300);
    let panels = PanelArgs {
        marginals: vec!["top".into()],
        ..PanelArgs::default()
    };
    let mut p = viewer(&t, Some(cells(300)), &panels, sheet(144.0, 108.0));
    let (plan, sheet) = p.viewing.clone().unwrap();
    focus_on(&mut p, Criterion::MinCoverage);
    let hist = p.screen_panel(&plan, sheet).unwrap();
    let head = hist.lines().next().unwrap();
    assert!(head.contains(r#"width="144" height="108""#), "{head}");
    let labels: Vec<String> = (crate::figure::labels(&hist).into_iter())
        .map(|l| l.1)
        .collect();
    assert!(labels.iter().any(|l| l.contains("coverage")), "{labels:?}");
    press(&mut p, KeyCode::Char('v'));
    let pair = p.screen_panel(&plan, sheet).unwrap();
    let (x, y, _) = p.view().shown_pairs()[p.pair_at].clone();
    let labels: Vec<String> = (crate::figure::labels(&pair).into_iter())
        .map(|l| l.1)
        .collect();
    for q in [x, y] {
        assert!(labels.iter().any(|l| l == q.title()), "{q:?}: {labels:?}");
    }
    assert_ne!(pair, hist);
}

/// `s` saves the panel on screen with the axis titles asked (`--xlab`,
/// `--ylab`), as `faba site-plot` draws them.
#[test]
fn s_saves_the_panel_on_screen_with_the_titles_asked() {
    let t = table(M6A, 300);
    let panels = PanelArgs {
        xlab: vec!["coverage=TITLE1 (x)".into()],
        ylab: vec!["sites=TITLE2".into()],
        ..PanelArgs::default()
    };
    let mut p = viewer(&t, Some(cells(300)), &panels, sheet(144.0, 108.0));
    let (plan, sheet) = p.viewing.clone().unwrap();
    focus_on(&mut p, Criterion::MinCoverage);
    let hist = p.screen_panel(&plan, sheet).unwrap();
    let labels: Vec<String> = (crate::figure::labels(&hist).into_iter())
        .map(|l| l.1)
        .collect();
    for want in ["TITLE1 (x)", "TITLE2"] {
        assert!(labels.iter().any(|l| l == want), "{want}: {labels:?}");
    }
}
