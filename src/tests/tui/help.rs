use super::*;
use crate::tui::{feed, screen_text};
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};

/// A view that keeps what it is handed and lists two keys, or none while
/// typing.
#[derive(Default)]
struct Typer {
    typing: bool,
    keys: Vec<KeyEvent>,
}

impl Screen for Typer {
    fn render(&mut self, _: &mut Frame) {}
    fn handle_key(&mut self, key: KeyEvent) {
        self.keys.push(key);
    }
    fn interrupt(&mut self) {}
    fn done(&self) -> bool {
        false
    }
}

impl View for Typer {
    fn help(&self) -> Option<KeyList> {
        (!self.typing).then(|| keys(&[("q", "quit"), ("PgUp/PgDn", "page")]))
    }
}

fn key(c: KeyCode) -> KeyEvent {
    KeyEvent::new(c, KeyModifiers::NONE)
}

#[test]
fn question_mark_opens_the_keys_and_question_mark_or_esc_closes_them() {
    let mut view = Typer::default();
    let mut h = Helped::new(&mut view);
    assert!(!screen_text(&mut h, 60, 12).contains("PgUp/PgDn"));
    feed(
        &mut h,
        KeyEvent::new(KeyCode::Char('?'), KeyModifiers::SHIFT),
    );
    assert!(h.is_open());
    let text = screen_text(&mut h, 60, 12);
    assert!(text.contains("PgUp/PgDn  page"), "{text}");
    assert!(text.contains("q          quit"), "{text}");
    // Open, other keys never reach the view.
    feed(&mut h, key(KeyCode::Char('q')));
    assert!(h.is_open());
    feed(&mut h, key(KeyCode::Esc));
    assert!(!h.is_open());
    feed(&mut h, key(KeyCode::Char('?')));
    feed(&mut h, key(KeyCode::Char('?')));
    assert!(!h.is_open());
    assert!(view.keys.is_empty(), "{:?}", view.keys);
}

#[test]
fn while_typing_question_mark_is_typed() {
    let mut view = Typer {
        typing: true,
        ..Typer::default()
    };
    let mut h = Helped::new(&mut view);
    feed(&mut h, key(KeyCode::Char('?')));
    assert!(!h.is_open());
    assert_eq!(view.keys, vec![key(KeyCode::Char('?'))]);
}

#[test]
fn a_long_list_takes_columns_on_a_short_screen() {
    struct Many;
    impl Screen for Many {
        fn render(&mut self, _: &mut Frame) {}
        fn handle_key(&mut self, _: KeyEvent) {}
        fn interrupt(&mut self) {}
        fn done(&self) -> bool {
            false
        }
    }
    impl View for Many {
        fn help(&self) -> Option<KeyList> {
            Some(
                (0..20)
                    .map(|i| (format!("k{i}"), format!("does {i}")))
                    .collect(),
            )
        }
    }
    let mut view = Many;
    let mut h = Helped::new(&mut view);
    feed(&mut h, key(KeyCode::Char('?')));
    let text = screen_text(&mut h, 80, 14);
    for i in 0..20 {
        assert!(text.contains(&format!("does {i}")), "k{i}: {text}");
    }
}
