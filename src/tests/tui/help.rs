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

#[test]
fn a_label_names_its_keys() {
    let ctrl = |c| KeyEvent::new(KeyCode::Char(c), KeyModifiers::CONTROL);
    let shift_tab = KeyEvent::new(KeyCode::BackTab, KeyModifiers::SHIFT);
    let chars = |s: &str| s.chars().map(|c| key(KeyCode::Char(c))).collect::<Vec<_>>();
    let mut arrows = vec![key(KeyCode::Up), key(KeyCode::Down)];
    arrows.extend(chars("kj"));
    assert_eq!(named("↑/↓ k/j"), arrows);
    assert_eq!(named("Tab/Shift+Tab"), vec![key(KeyCode::Tab), shift_tab]);
    assert_eq!(named("0-9 e"), chars("0123456789e"));
    assert_eq!(named("-/+ , . ="), chars("-+,.="));
    assert_eq!(named("[ ] /"), chars("[]/"));
    assert_eq!(named("Ctrl+E"), vec![ctrl('e')]);
    assert_eq!(
        named("g/Esc"),
        vec![key(KeyCode::Char('g')), key(KeyCode::Esc)]
    );
    assert!(named("click").is_empty());
}

/// Takes `q`, and `x` only while `x_on`; lists `x` when `lists_x`.
struct Two {
    x_on: bool,
    lists_x: bool,
    done: bool,
    x: usize,
}

impl Screen for Two {
    fn render(&mut self, f: &mut Frame) {
        f.render_widget(Line::from(format!("x {}", self.x)), f.area());
    }
    fn handle_key(&mut self, key: KeyEvent) {
        match key.code {
            KeyCode::Char('q') => self.done = true,
            KeyCode::Char('x') if self.x_on => self.x += 1,
            _ => {}
        }
    }
    fn interrupt(&mut self) {}
    fn done(&self) -> bool {
        self.done
    }
}

impl View for Two {
    fn help(&self) -> Option<KeyList> {
        let mut k = vec![("q", "quit")];
        if self.lists_x {
            k.push(("x", "more"));
        }
        Some(keys(&k))
    }
}

#[test]
fn the_check_finds_a_key_listed_but_inert_and_one_bound_but_unlisted() {
    let two = |x_on, lists_x| {
        move || Two {
            x_on,
            lists_x,
            done: false,
            x: 0,
        }
    };
    assert_lists_its_keys("both", two(true, true));
    assert_lists_its_keys("neither", two(false, false));
    let x = || vec!["Char('x')".to_string()];
    assert_eq!(disagree(two(false, true)), (x(), vec![]));
    assert_eq!(disagree(two(true, false)), (vec![], x()));
}
