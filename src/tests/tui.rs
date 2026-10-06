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
    let mut f = Find::default();
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

#[test]
fn centered_fits_inside_small_areas() {
    let r = centered(Rect::new(0, 0, 10, 4), 40, 12);
    assert_eq!((r.width, r.height), (10, 4));
    assert_eq!(first_visible(9, 10, 4), 6);
}

#[test]
fn output_folders_are_numbered_and_checked() {
    let tmp = tempfile::tempdir().unwrap();
    let base = tmp.path();
    assert_eq!(next_free(base, "out"), base.join("out"));
    std::fs::create_dir(base.join("out")).unwrap();
    assert_eq!(next_free(base, "out"), base.join("out2"));
    let out = base.join("out").to_string_lossy().into_owned();
    assert_eq!(output_problem(&out), None, "an empty folder is fine");
    std::fs::write(base.join("out/x"), b"").unwrap();
    assert!(output_problem(&out)
        .unwrap()
        .contains("already contains files"));
    let file = base.join("out/x").to_string_lossy().into_owned();
    assert!(output_problem(&file).unwrap().contains("is a file"));
    assert!(output_problem("").is_some());
}

/// A view that keeps what it is handed.
#[derive(Default)]
struct Keeps {
    apply: bool,
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
