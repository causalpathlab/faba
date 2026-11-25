//! The `?` pop-up every view has: each key the view takes as it is now,
//! the ones its footer leaves out among them, centred over the view until
//! `?` or Esc closes it. A view lists its keys ([`View::help`]); while it
//! is typing (a find line, a name, a path) it lists none and `?` is typed.

use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseEvent};
use ratatui::text::{Line, Span};
use ratatui::Frame;

use super::ui::{DIM, HIGHLIGHT, PLAIN};
use super::{is_click_or_wheel, Screen, View};

/// A view's keys as the pop-up lists them: the key, then what it does.
pub type KeyList = Vec<(String, String)>;

/// `pairs` (a footer's, say) as a [`KeyList`].
pub fn keys(pairs: &[(&str, &str)]) -> KeyList {
    pairs
        .iter()
        .map(|&(k, what)| (k.to_string(), what.to_string()))
        .collect()
}

/// Whether `key` is `?` as typed, with Shift or nothing else held.
pub fn is_help(key: &KeyEvent) -> bool {
    key.code == KeyCode::Char('?') && (key.modifiers - KeyModifiers::SHIFT).is_empty()
}

/// A view with its `?` pop-up: keys go to the view while the pop-up is
/// closed; open, `?` or Esc closes it and nothing else reaches the view.
pub struct Helped<'v, V> {
    view: &'v mut V,
    open: bool,
}

impl<'v, V: View> Helped<'v, V> {
    pub fn new(view: &'v mut V) -> Self {
        Self { view, open: false }
    }

    /// Whether the pop-up is showing.
    #[cfg(test)]
    pub fn is_open(&self) -> bool {
        self.open
    }
}

impl<V: View> View for Helped<'_, V> {
    fn stray_enter(&mut self) {
        if !self.open {
            self.view.stray_enter();
        }
    }

    fn takes_apply(&self) -> bool {
        !self.open && self.view.takes_apply()
    }

    fn takes_chord(&self, key: &KeyEvent) -> bool {
        !self.open && self.view.takes_chord(key)
    }
}

impl<V: View> Screen for Helped<'_, V> {
    fn render(&mut self, frame: &mut Frame) {
        self.view.render(frame);
        if !self.open {
            return;
        }
        // Listed as the view is now; a view that started typing (it cannot
        // while the pop-up is open) would list nothing.
        let list = self.view.help().unwrap_or_default();
        draw(frame, &list);
    }

    fn handle_key(&mut self, key: KeyEvent) {
        if self.open {
            if is_help(&key) || key.code == KeyCode::Esc {
                self.open = false;
            }
        } else if is_help(&key) && self.view.help().is_some() {
            self.open = true;
        } else {
            self.view.handle_key(key);
        }
    }

    fn interrupt(&mut self) {
        self.view.interrupt();
    }

    fn done(&self) -> bool {
        self.view.done()
    }

    fn quit(&mut self) {
        self.view.quit();
    }

    fn pending_work(&self) -> Option<String> {
        self.view.pending_work()
    }

    fn do_work(&mut self) {
        self.view.do_work();
    }

    fn tick(&mut self) -> bool {
        self.view.tick()
    }

    fn takes_mouse(&self) -> bool {
        self.view.takes_mouse()
    }

    /// Open, a click closes the pop-up; the wheel does nothing.
    fn mouse(&mut self, event: MouseEvent) -> bool {
        if !self.open {
            return self.view.mouse(event);
        }
        if is_click_or_wheel(&event) && super::wheel_key(&event).is_none() {
            self.open = false;
            return true;
        }
        false
    }

    fn reports_chords(&self) -> bool {
        self.view.reports_chords()
    }

    fn take_notice(&mut self) -> Option<String> {
        self.view.take_notice()
    }

    fn status(&self) -> Option<String> {
        self.view.status()
    }
}

/// The pop-up listing `list`, in as many columns as the screen's height
/// asks for.
fn draw(frame: &mut Frame, list: &KeyList) {
    let area = frame.area();
    // Two rows of border, a blank and the closing note.
    let rows = usize::from(area.height.saturating_sub(4)).max(1);
    let columns = list.len().div_ceil(rows).max(1);
    let per = list.len().div_ceil(columns).max(1);
    let key_w = list
        .iter()
        .map(|(k, _)| k.chars().count())
        .max()
        .unwrap_or(0);
    let entry = |(k, what): &(String, String)| {
        vec![
            Span::styled(format!(" {k:<key_w$}  "), HIGHLIGHT),
            Span::styled(what.clone(), PLAIN),
        ]
    };
    let width = |i: usize| -> usize {
        list.iter()
            .skip(i * per)
            .take(per)
            .map(|(_, what)| what.chars().count())
            .max()
            .unwrap_or(0)
    };
    let mut lines: Vec<Line<'static>> = (0..per)
        .map(|r| {
            let mut spans = Vec::new();
            for c in 0..columns {
                let Some(pair) = list.get(c * per + r) else {
                    break;
                };
                let mut e = entry(pair);
                if c + 1 < columns {
                    let pad = width(c) - pair.1.chars().count() + 2;
                    e.push(Span::raw(" ".repeat(pad)));
                }
                spans.extend(e);
            }
            Line::from(spans)
        })
        .collect();
    lines.push(Line::from(""));
    lines.push(Line::from(Span::styled(" ? or Esc closes", DIM)));
    super::popup(frame, area, " keys ", lines);
}

/// The screen after `?` in `view` as it is, `w` x `h`: the pop-up over
/// it; none when `?` opened nothing (the view was typing, and took it).
#[cfg(test)]
pub fn shown<V: View>(view: &mut V, w: u16, h: u16) -> Option<String> {
    let mut helped = Helped::new(view);
    super::feed(
        &mut helped,
        KeyEvent::new(KeyCode::Char('?'), KeyModifiers::NONE),
    );
    helped.open.then(|| super::screen_text(&mut helped, w, h))
}

#[cfg(test)]
#[path = "../tests/tui/help.rs"]
mod tests;
