//! The form `faba modules` asks the fit's settings in: prefilled from the
//! flags and the preset, refilled on a preset switch, checked as typed, and
//! what Ctrl+R and Esc leave.

use super::model::modules;
use super::*;
use crate::cluster::browse::ask::{self, Asked, FitForm, ROWS};
use crate::tui::{apply_key, feed, screen_text};

/// The form over `words`' settings, as `faba modules` would ask them.
fn form(words: &[&str]) -> FitForm {
    let tmp = tempfile::tempdir().unwrap();
    let args = modules(tmp.path(), words).blockmodel().unwrap();
    FitForm::new(&args, "fitting", "no block tables", "open the old tables")
}

fn press(view: &mut FitForm, code: KeyCode) {
    feed(view, KeyEvent::new(code, KeyModifiers::NONE));
}

/// Move to the row of `long`.
fn go(view: &mut FitForm, long: &str) {
    press(view, KeyCode::Home);
    let at = ROWS.iter().position(|r| *r == long).unwrap();
    (0..at).for_each(|_| press(view, KeyCode::Down));
}

/// Type `text` over the row of `long`.
fn type_into(view: &mut FitForm, long: &str, text: &str) {
    go(view, long);
    press(view, KeyCode::Enter);
    (0..32).for_each(|_| press(view, KeyCode::Backspace));
    text.chars().for_each(|c| press(view, KeyCode::Char(c)));
    press(view, KeyCode::Enter);
}

#[test]
fn the_form_is_prefilled_from_the_flags_then_the_preset() {
    let mut view = form(&["--chains", "2", "--gene-roles", "4"]);
    for (long, value) in [
        ("model", "blocks"),
        ("chains", "2"),
        ("gene-roles", "4"),
        ("shift", "on"),
        ("profile-strength", "inf"),
        ("global", "symmetric"),
        ("sweeps", "150"),
        ("cell-roles", "8"),
    ] {
        assert_eq!(view.value(long), value, "--{long}");
    }
    let shown = screen_text(&mut view, 100, 30);
    for long in ROWS {
        assert!(shown.contains(long), "--{long} not shown:\n{shown}");
    }
    assert!(shown.contains("Ctrl+R") && shown.contains("Esc"), "{shown}");
}

#[test]
fn a_preset_switch_refills_only_the_untouched_rows() {
    let mut view = form(&["--gene-roles", "4"]);
    type_into(&mut view, "cell-sd", "0.5");
    go(&mut view, "model");
    press(&mut view, KeyCode::Right);
    assert_eq!(view.value("model"), "link");
    assert_eq!(view.value("shift"), "off");
    assert_eq!(view.value("profile-strength"), "5");
    assert_eq!(view.value("global"), "hdp");
    // Given on the line, and typed: kept.
    assert_eq!(view.value("gene-roles"), "4");
    assert_eq!(view.value("cell-sd"), "0.5");
}

#[test]
fn an_invalid_combination_is_said_and_blocks_the_fit() {
    let mut view = form(&[]);
    go(&mut view, "shift");
    press(&mut view, KeyCode::Right);
    assert_eq!(view.value("shift"), "off");
    let shown = screen_text(&mut view, 120, 30);
    assert!(shown.contains("need the block shift"), "{shown}");
    feed(&mut view, apply_key());
    assert!(!view.done());
    // The blocked note says why in full, not only beside the row.
    let shown = screen_text(&mut view, 120, 30);
    let note = shown
        .lines()
        .find(|l| l.contains("is blocked"))
        .unwrap_or("");
    assert!(note.contains("need the block shift"), "{shown}");
    // Back on: the fit can start.
    press(&mut view, KeyCode::Right);
    feed(&mut view, apply_key());
    assert!(view.done());
}

#[test]
fn ctrl_r_fits_with_the_values_edited() {
    let mut view = form(&[]);
    type_into(&mut view, "chains", "2");
    type_into(&mut view, "profile-strength", "10");
    // Enter never applies: it only types.
    assert!(!view.done());
    feed(&mut view, apply_key());
    assert!(view.done());
    let args = view.chosen().unwrap();
    assert_eq!(args.chains, 2);
    assert_eq!(args.model, ModelKind::Blocks);
    assert_eq!(args.profile_strength.map(|s| s.0), Some(10.0));
    assert!(args.given("chains"));
    let saved = crate::cluster::command::command(&args);
    let at = saved.iter().position(|w| w == "--chains").unwrap();
    assert_eq!(saved[at + 1], "2", "{saved:?}");
    assert!(saved.contains(&"blocks".to_string()), "{saved:?}");
}

#[test]
fn esc_leaves_without_settings() {
    let mut view = form(&[]);
    press(&mut view, KeyCode::Esc);
    assert!(view.done());
    assert!(view.chosen().is_none());
    assert!(matches!(view.outcome(), Some(Asked::Back)));
}

#[test]
fn ctrl_c_quits_rather_than_going_back() {
    let mut view = form(&[]);
    view.interrupt();
    assert!(view.done());
    assert!(matches!(view.outcome(), Some(Asked::Quit)));
}

#[test]
fn the_footer_says_what_esc_leaves_to() {
    let mut view = form(&[]);
    let shown = screen_text(&mut view, 160, 30);
    assert!(shown.contains("open the old tables"), "{shown}");
    assert!(shown.contains("Ctrl+C"), "{shown}");
}

#[test]
fn ctrl_r_while_typing_fits_with_the_text() {
    let mut view = form(&[]);
    go(&mut view, "chains");
    press(&mut view, KeyCode::Enter);
    (0..32).for_each(|_| press(&mut view, KeyCode::Backspace));
    press(&mut view, KeyCode::Char('2'));
    feed(&mut view, apply_key());
    assert!(view.done());
    assert_eq!(view.chosen().unwrap().chains, 2);
}

#[test]
fn a_flag_given_at_its_default_is_not_marked_changed() {
    let mut view = form(&[
        "--cell-sd",
        "1",
        "--model",
        "link",
        "--profile-strength",
        "5",
    ]);
    assert_eq!(view.value("cell-sd"), "1");
    assert_eq!(view.value("profile-strength"), "5");
    let shown = screen_text(&mut view, 120, 30);
    assert!(shown.contains("1 changed"), "{shown}");
}

#[test]
fn a_held_enter_says_what_plain_enter_does() {
    let mut view = form(&[]);
    feed(
        &mut view,
        KeyEvent::new(KeyCode::Enter, KeyModifiers::SHIFT),
    );
    let shown = screen_text(&mut view, 120, 30);
    assert!(!shown.contains("Enter does nothing"), "{shown}");
    assert!(shown.contains("Enter types"), "{shown}");
}

#[test]
fn batch_mode_skips_the_form() {
    let tmp = tempfile::tempdir().unwrap();
    let args = modules(tmp.path(), &["--chains", "3"])
        .blockmodel()
        .unwrap();
    let kept = ask::chosen(args, false, |_| panic!("the form was shown")).unwrap();
    let Asked::Fit(kept) = kept else {
        panic!("not fitted")
    };
    assert_eq!(kept.chains, 3);
}
