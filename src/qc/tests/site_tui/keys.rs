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
        "q",
    ] {
        assert!(keys.contains(key), "{key}:\n{keys}");
    }
    assert!(!keys.contains("PgUp/PgDn"), "only in the genes:\n{keys}");
    press(&mut p, KeyCode::Tab);
    let keys = crate::tui::help::shown(&mut p, 120, 50).expect("? opens the keys");
    assert!(keys.contains("PgUp/PgDn"), "{keys}");
    press(&mut p, KeyCode::Char('/'));
    assert!(matches!(p.mode, Mode::Find));
    assert!(crate::tui::help::shown(&mut p, 120, 50).is_none());
    assert_eq!(p.list.find, "?", "? is typed into the find");
}

/// Two modalities with genes and cells, on the second knob nudged up a
/// step and the second gene, then `keys`: so that every move, step and
/// reset moves.
fn mid_way<'a>(tables: &'a [SiteTable; 2], keys: &[KeyCode]) -> SitePicker<'a> {
    let start = SiteFilterArgs::default_values();
    let views = (tables.iter())
        .map(|t| SiteView::new(t, Some(cells(t.len())), &start))
        .collect();
    let mut p = SitePicker::new("x", views, start, Meta::Unavailable("none".into()));
    for code in [KeyCode::Down, KeyCode::Char('+'), KeyCode::Char(']')]
        .iter()
        .chain(keys)
    {
        press(&mut p, *code);
    }
    p
}

#[test]
fn the_keys_listed_are_the_keys_bound_in_each_state() {
    use crate::tui::help::assert_lists_its_keys;
    let tables = [with_genes(table(M6A, 300)), with_genes(table(ATOI, 300))];
    let (tab, v, x) = (KeyCode::Tab, KeyCode::Char('v'), KeyCode::Char('x'));
    for (what, keys) in [
        ("thresholds", &[][..]),
        ("genes", &[tab][..]),
        // A panel's axis chosen by hand, so `A` has pairs to put back.
        ("the panels", &[v, tab, tab, x][..]),
        ("thresholds in 2D", &[v, tab, tab, x, tab][..]),
        ("genes in 2D", &[v, tab, tab, x, tab, tab][..]),
        ("the confirm", &[][..]),
        ("leaving", &[KeyCode::Char('q')][..]),
        ("printing", &[KeyCode::Char('p')][..]),
    ] {
        assert_lists_its_keys(what, || {
            let mut p = mid_way(&tables, keys);
            if what == "the confirm" {
                apply_key(&mut p);
            }
            p
        });
    }
    assert_lists_its_keys("read only", || {
        let mut p = mid_way(&tables, &[]);
        p.viewing = Some((
            crate::qc::site_plot::PanelArgs::default().plan().unwrap(),
            crate::figure::style::Sheet::OWN,
        ));
        p
    });
    // No other modality for `m` and `M` to switch to.
    let one = with_genes(table(M6A, 300));
    assert_lists_its_keys("one modality", || {
        let mut p = picker(&one, None, SiteFilterArgs::default_values());
        for code in [KeyCode::Down, KeyCode::Char('+'), KeyCode::Char(']')] {
            press(&mut p, code);
        }
        p
    });
}

/// On a knob of one scale the footer names `y` alone: `x` would step
/// nothing.
#[test]
fn the_footer_names_only_the_scale_key_that_acts() {
    let t = with_genes(table(M6A, 300));
    let make = || {
        let mut p = picker(&t, None, SiteFilterArgs::default_values());
        focus_on(&mut p, Criterion::MinLogOdds);
        // Off its start and off the first gene, so `r` and `[` act.
        press(&mut p, KeyCode::Char('+'));
        press(&mut p, KeyCode::Char(']'));
        p
    };
    let said = footer(&make());
    assert!(said.contains("y scale") && !said.contains("x/y"), "{said}");
    crate::tui::help::assert_lists_its_keys("a knob of one scale", make);
    let mut p = make();
    focus_on(&mut p, Criterion::MinCoverage);
    assert!(footer(&p).contains("x/y scale"), "{}", footer(&p));
}
