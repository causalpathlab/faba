use super::*;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};

#[test]
fn the_apply_key_is_ctrl_r_and_nothing_else() {
    let k = |c, m| KeyEvent::new(c, m);
    let ctrl_shift = KeyModifiers::CONTROL | KeyModifiers::SHIFT;
    assert!(is_apply(&k(KeyCode::Char('r'), KeyModifiers::CONTROL)));
    assert!(is_apply(&k(KeyCode::Char('R'), ctrl_shift)));
    assert!(is_apply(&apply_key()));
    for m in [
        KeyModifiers::NONE,
        KeyModifiers::SHIFT,
        KeyModifiers::CONTROL,
    ] {
        assert!(!is_apply(&k(KeyCode::Enter, m)), "{m:?}+Enter");
    }
    assert!(!is_apply(&k(KeyCode::Char('j'), KeyModifiers::CONTROL)));
    for c in ['G', 'g', 'y', 'r', 'R'] {
        assert!(!is_apply(&k(KeyCode::Char(c), KeyModifiers::NONE)), "{c}");
    }
    // Typed into a find, it types nothing.
    let mut f = legume_tui::browse::Find::default();
    assert!(!f.key(&apply_key()) && f.text.is_empty());
}

#[test]
fn modified_keys_are_stray_but_the_apply_key_and_capitals_are_not() {
    let k = |c, m| KeyEvent::new(c, m);
    for m in [
        KeyModifiers::SHIFT,
        KeyModifiers::ALT,
        KeyModifiers::SUPER,
        KeyModifiers::CONTROL,
    ] {
        assert!(is_stray(&k(KeyCode::Enter, m)), "{m:?}+Enter");
    }
    assert!(is_stray(&k(KeyCode::Char('s'), KeyModifiers::CONTROL)));
    assert!(is_stray(&k(KeyCode::Char('j'), KeyModifiers::CONTROL)));
    assert!(is_stray(&k(KeyCode::Char('q'), KeyModifiers::ALT)));
    assert!(!is_stray(&k(KeyCode::Char('r'), KeyModifiers::CONTROL)));
    assert!(!is_stray(&k(KeyCode::Enter, KeyModifiers::NONE)));
    assert!(!is_stray(&k(KeyCode::Char('R'), KeyModifiers::SHIFT)));
    assert!(!is_stray(&k(KeyCode::BackTab, KeyModifiers::SHIFT)));
}

/// A view that keeps what it is handed.
#[derive(Default)]
struct Keeps {
    apply: bool,
    /// The chords it takes.
    chords: Vec<KeyEvent>,
    keys: Vec<KeyEvent>,
    stray: usize,
}

impl Screen for Keeps {
    fn render(&mut self, _: &mut Frame) {}
    fn handle_key(&mut self, key: KeyEvent) {
        self.keys.push(key);
    }
    fn interrupt(&mut self) {}
    fn done(&self) -> bool {
        false
    }
}

impl View for Keeps {
    fn stray_enter(&mut self) {
        self.stray += 1;
    }
    fn takes_apply(&self) -> bool {
        self.apply
    }
    fn takes_chord(&self, key: &KeyEvent) -> bool {
        self.chords.contains(key)
    }
}

#[test]
fn feed_hands_on_only_what_a_view_takes() {
    let k = |c, m| KeyEvent::new(c, m);
    for apply in [false, true] {
        let mut v = Keeps {
            apply,
            ..Keeps::default()
        };
        feed(&mut v, apply_key());
        feed(&mut v, k(KeyCode::Char('s'), KeyModifiers::CONTROL));
        feed(&mut v, k(KeyCode::Enter, KeyModifiers::SHIFT));
        feed(&mut v, k(KeyCode::Char('r'), KeyModifiers::NONE));
        let want = if apply {
            vec![apply_key(), k(KeyCode::Char('r'), KeyModifiers::NONE)]
        } else {
            vec![k(KeyCode::Char('r'), KeyModifiers::NONE)]
        };
        assert_eq!(v.keys, want, "apply taken: {apply}");
        assert_eq!(v.stray, 1);
    }
}

#[test]
fn a_view_gets_only_the_chords_it_names() {
    let k = |c, m| KeyEvent::new(c, m);
    let ctrl_e = k(KeyCode::Char('e'), KeyModifiers::CONTROL);
    let ctrl_s = k(KeyCode::Char('s'), KeyModifiers::CONTROL);
    let ctrl_enter = k(KeyCode::Enter, KeyModifiers::CONTROL);
    // Taking none, as every view does by default: all dropped.
    let mut none = Keeps::default();
    for key in [ctrl_e, ctrl_s, ctrl_enter] {
        feed(&mut none, key);
    }
    assert!(none.keys.is_empty());
    assert_eq!(none.stray, 1, "Ctrl+Enter still says what applies");
    // Naming Ctrl+E: that one comes, the others are dropped as before.
    let mut named = Keeps {
        chords: vec![ctrl_e],
        ..Keeps::default()
    };
    for key in [ctrl_e, ctrl_s, ctrl_enter] {
        feed(&mut named, key);
    }
    assert_eq!(named.keys, vec![ctrl_e]);
    assert_eq!(named.stray, 1);
}

/// A `faba kill` is known through the context faba's errors gather, and
/// through the I/O error [`run_view`] carries it in.
#[test]
fn a_stop_is_known_through_the_context_it_gathers() {
    let stopped = || legume_tui::view::Stopped::new(&TOOL);
    let e = anyhow::Error::from(stopped()).context("picking the inputs");
    assert!(is_stopped(&e));
    let io = std::io::Error::other(stopped());
    assert!(is_stopped(&anyhow::Error::from(io).context("the view")));
    assert!(!is_stopped(&anyhow::anyhow!("no such file")));
}
