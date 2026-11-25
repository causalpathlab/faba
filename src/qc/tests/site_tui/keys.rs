//! The keys as the shared scheme has them: `w` the weight, `v` the 2D
//! look, Esc and `n` back from the confirm, the footer of each panel, and
//! the `?` pop-up.

use super::*;

fn footer(p: &SitePicker) -> String {
    p.footer_line().to_string()
}

#[test]
fn w_switches_the_weight_and_c_does_not() {
    let t = table(M6A, 300);
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    let start = p.weight;
    press(&mut p, KeyCode::Char('c'));
    assert_eq!(p.weight, start, "c is no longer the weight");
    press(&mut p, KeyCode::Char('w'));
    assert_ne!(p.weight, start);
    assert!(footer(&p).contains("w weight"), "{}", footer(&p));
}

#[test]
fn v_switches_the_2d_look_and_d_still_does_unshown() {
    let t = table(M6A, 300);
    let mut p = picker(&t, Some(cells(300)), SiteFilterArgs::default_values());
    press(&mut p, KeyCode::Char('v'));
    assert!(p.two_d);
    press(&mut p, KeyCode::Char('v'));
    assert!(!p.two_d);
    press(&mut p, KeyCode::Char('d'));
    assert!(p.two_d, "d is kept as an alias");
    let said = footer(&p);
    assert!(said.contains("v histograms"), "{said}");
    assert!(
        !said.contains("d histograms") && !said.contains("d 2D"),
        "{said}"
    );
}

#[test]
fn the_confirm_goes_back_on_esc_or_n_but_not_q() {
    let t = table(M6A, 300);
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    apply_key(&mut p);
    assert!(matches!(p.mode, Mode::Confirm));
    assert!(footer(&p).contains("Esc/n back"), "{}", footer(&p));
    press(&mut p, KeyCode::Char('q'));
    assert!(matches!(p.mode, Mode::Confirm), "q is quit, never back");
    press(&mut p, KeyCode::Char('n'));
    assert!(matches!(p.mode, Mode::Browse));
    apply_key(&mut p);
    press(&mut p, KeyCode::Esc);
    assert!(matches!(p.mode, Mode::Browse));
}

/// Each panel's footer names what works there: Shift+Tab back, PgUp/PgDn
/// in the genes, `A` only in 2D, and no scales on the panels, where `x`
/// and `y` are the axes.
#[test]
fn each_panels_footer_names_what_works_there() {
    let t = table(M6A, 300);
    let mut p = picker(&t, Some(cells(300)), SiteFilterArgs::default_values());
    let said = footer(&p);
    assert!(said.contains("Tab/Shift+Tab panel"), "{said}");
    assert!(said.contains("x/y scale"), "{said}");
    assert!(!said.contains("auto pairs"), "no A out of 2D: {said}");
    press(&mut p, KeyCode::Tab);
    assert!(footer(&p).contains("PgUp/PgDn page"), "{}", footer(&p));
    press(&mut p, KeyCode::Char('v'));
    assert!(footer(&p).contains("A auto pairs"), "{}", footer(&p));
    press(&mut p, KeyCode::Tab);
    assert_eq!(p.panel, Panel::Pairs);
    let said = footer(&p);
    assert!(!said.contains("x/y scale"), "{said}");
    assert!(said.contains("x/y X/Y axis"), "{said}");
}

#[test]
fn question_mark_lists_the_keys_but_types_in_the_find() {
    let t = table(M6A, 300);
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    let keys = crate::tui::help::shown(&mut p, 120, 50).expect("? opens the keys");
    for key in [
        "Tab/Shift+Tab",
        "↑/↓ k/j",
        "←/→ h/l",
        "w",
        "v",
        "Ctrl+E",
        "Ctrl+R",
        "PgUp/PgDn",
        "q",
    ] {
        assert!(keys.contains(key), "{key}:\n{keys}");
    }
    press(&mut p, KeyCode::Tab);
    press(&mut p, KeyCode::Char('/'));
    assert!(matches!(p.mode, Mode::Find));
    assert!(crate::tui::help::shown(&mut p, 120, 50).is_none());
    assert_eq!(p.list.find, "?", "? is typed into the find");
}
