use super::*;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};

#[test]
fn go_is_shift_enter_or_capital_g() {
    let k = |c, m| KeyEvent::new(c, m);
    assert!(is_go(&k(KeyCode::Enter, KeyModifiers::SHIFT)));
    assert!(is_go(&k(KeyCode::Char('G'), KeyModifiers::NONE)));
    assert!(!is_go(&k(KeyCode::Enter, KeyModifiers::NONE)));
    assert!(!is_go(&k(KeyCode::Char('g'), KeyModifiers::NONE)));
}

#[test]
fn centered_fits_inside_small_areas() {
    let r = centered(Rect::new(0, 0, 10, 4), 40, 12);
    assert_eq!((r.width, r.height), (10, 4));
    assert_eq!(first_visible(9, 10, 4), 6);
}

#[test]
fn shift_enter_is_asked_once_and_released() {
    let mut s = ShiftEnter::default();
    s.arm();
    assert!(!s.on, "not wanted, never asked");
}
