//! The `?` pop-up every view has: each key the view takes as it is now,
//! the ones its footer leaves out among them, centred over the view until
//! `?` or Esc closes it. A view lists its keys ([`View::help`]); while it
//! is typing (a find line, a name, a path) it lists none and `?` is typed.

use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseEvent};
use ratatui::text::{Line, Span};
use ratatui::Frame;

use super::live::Control;
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

/// `keys` with `s` (a figure view's save key) saying `what` it saves.
pub fn say_save(keys: &mut [(&str, &'static str)], what: &'static str) {
    for (k, does) in keys.iter_mut() {
        if *k == "s" {
            *does = what;
        }
    }
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
    /// Asked of the view at each tick, from outside it (see
    /// [`super::live`]): whether that changed it.
    served: Option<Served<'v, V>>,
}

/// What a view is asked at each tick, and how it answers: whether it
/// changed.
type Served<'v, V> = (&'v Control, fn(&Control, &mut V) -> bool);

impl<'v, V: View> Helped<'v, V> {
    pub fn new(view: &'v mut V) -> Self {
        Self {
            view,
            open: false,
            served: None,
        }
    }

    /// The view asked what `control` is asked at each tick too, answered
    /// by `serve` (see [`super::live`]).
    pub fn serving(mut self, control: &'v Control, serve: fn(&Control, &mut V) -> bool) -> Self {
        self.served = Some((control, serve));
        self
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
        let ticked = self.view.tick();
        let served = (self.served).is_some_and(|(control, serve)| serve(control, self.view));
        ticked || served
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

/// The keys `label` names as a pop-up writes them: words apart, keys in a
/// word apart by `/` (`↑/↓ k/j`, `Tab/Shift+Tab`, `[ ]`), a range of
/// digits (`0-9`), a chord (`Ctrl+E`). `click`, `Ctrl-drag` and `type`
/// name no key; `arrows` the four.
#[cfg(test)]
pub fn named(label: &str) -> Vec<KeyEvent> {
    let plain = |c| KeyEvent::new(c, KeyModifiers::NONE);
    let mut out = Vec::new();
    for word in label.split_whitespace() {
        // A `/` alone is the key; within a word it parts two.
        let parts: Vec<&str> = match word {
            "/" => vec!["/"],
            _ => word.split('/').collect(),
        };
        for part in parts {
            let code = match part {
                "click" | "type" | "Ctrl-drag" => continue,
                "arrows" => {
                    let arrows = [KeyCode::Up, KeyCode::Down, KeyCode::Left, KeyCode::Right];
                    out.extend(arrows.map(plain));
                    continue;
                }
                "Shift+Tab" => {
                    out.push(KeyEvent::new(KeyCode::BackTab, KeyModifiers::SHIFT));
                    continue;
                }
                _ if part.starts_with("Ctrl+") && part.chars().count() == 6 => {
                    let c = part.chars().last().unwrap().to_ascii_lowercase();
                    out.push(KeyEvent::new(KeyCode::Char(c), KeyModifiers::CONTROL));
                    continue;
                }
                _ if part.len() == 3 && part.as_bytes()[1] == b'-' && part != "-" => {
                    let b = part.as_bytes();
                    out.extend((b[0]..=b[2]).map(|c| plain(KeyCode::Char(char::from(c)))));
                    continue;
                }
                "↑" => KeyCode::Up,
                "↓" => KeyCode::Down,
                "←" => KeyCode::Left,
                "→" => KeyCode::Right,
                "Tab" => KeyCode::Tab,
                "Esc" => KeyCode::Esc,
                "Enter" => KeyCode::Enter,
                "Space" => KeyCode::Char(' '),
                "Backspace" => KeyCode::Backspace,
                "Del" => KeyCode::Delete,
                "Home" => KeyCode::Home,
                "End" => KeyCode::End,
                "PgUp" => KeyCode::PageUp,
                "PgDn" => KeyCode::PageDown,
                _ if part.chars().count() == 1 => KeyCode::Char(part.chars().next().unwrap()),
                _ => panic!("{part:?} in {label:?} names no key"),
            };
            out.push(plain(code));
        }
    }
    out
}

/// Every key a test tries on a view: the printable characters but `?`
/// (the pop-up's), the named keys, and Ctrl with each letter but `c`.
#[cfg(test)]
fn every_key() -> Vec<KeyEvent> {
    let plain = |c| KeyEvent::new(c, KeyModifiers::NONE);
    let mut keys: Vec<KeyEvent> = (' '..='~')
        .filter(|&c| c != '?')
        .map(|c| plain(KeyCode::Char(c)))
        .collect();
    keys.extend(
        [
            KeyCode::Esc,
            KeyCode::Enter,
            KeyCode::Tab,
            KeyCode::Backspace,
            KeyCode::Delete,
            KeyCode::Up,
            KeyCode::Down,
            KeyCode::Left,
            KeyCode::Right,
            KeyCode::Home,
            KeyCode::End,
            KeyCode::PageUp,
            KeyCode::PageDown,
        ]
        .map(plain),
    );
    keys.push(KeyEvent::new(KeyCode::BackTab, KeyModifiers::SHIFT));
    // Ctrl+C is the session's: it interrupts every view.
    keys.extend(
        ('a'..='z')
            .filter(|&c| c != 'c')
            .map(|c| KeyEvent::new(KeyCode::Char(c), KeyModifiers::CONTROL)),
    );
    keys
}

/// How [`super::enter_hint`] starts.
#[cfg(test)]
const HINT: &str = "Enter does nothing here";

/// What a key can be seen to have done to a view: what it draws (each
/// cell's look too: a focus may be a colour only), says and is about to
/// do.
#[cfg(test)]
pub fn seen<V: View>(view: &mut V) -> String {
    let backend = ratatui::backend::TestBackend::new(160, 50);
    let mut t = ratatui::Terminal::new(backend).expect("a test backend draws");
    t.draw(|f| view.render(f)).expect("a test backend draws");
    let notice = view.take_notice();
    format!(
        "{:?} {} {:?} {:?} {notice:?}",
        t.backend().buffer(),
        view.done(),
        view.status(),
        view.pending_work()
    )
}

/// The keys of the view `make` builds, as it is then, that its `?` list
/// or its footer ([`View::footer_keys`]) and its key handler disagree on:
/// listed (in either) but doing nothing, then bound but not in the `?`
/// list. A key is bound when it changes what `look` sees of the view,
/// each tried on a view of its own (the apply key too, which is no chord
/// to a view); a chord when the view takes it. The keys `untried` names
/// (as a pop-up writes them: one that downloads, say) are taken as bound,
/// unpressed. While the view lists nothing (it is typing) only its footer
/// is checked; panics when it names no key at all there, or `look` sees
/// it change by itself.
#[cfg(test)]
pub fn disagree_by<V: View>(
    mut make: impl FnMut() -> V,
    look: impl Fn(&mut V) -> String,
    untried: &str,
) -> (Vec<String>, Vec<String>) {
    let untried = named(untried);
    let first = make();
    let names = |list: Option<KeyList>| -> Option<Vec<KeyEvent>> {
        list.map(|l| {
            l.iter()
                .flat_map(|(k, _)| named(k))
                .filter(|k| !is_help(k))
                .collect()
        })
    };
    let help = names(first.help());
    let footer = names(first.footer_keys()).unwrap_or_default();
    assert!(
        help.is_some() || !footer.is_empty(),
        "the view lists its keys here"
    );
    let mut steady = make();
    assert_eq!(look(&mut steady), look(&mut steady), "draws the same twice");
    let (mut inert, mut unlisted) = (Vec::new(), Vec::new());
    for key in every_key() {
        let in_help = help.as_ref().is_some_and(|h| h.contains(&key));
        let listed = in_help || footer.contains(&key);
        // Typing, every key types: only the footer's are asked after.
        if help.is_none() && !listed {
            continue;
        }
        let mut view = make();
        let bound = if untried.contains(&key) {
            true
        } else if super::is_stray(&key) {
            view.takes_chord(&key)
        } else {
            let before = look(&mut view);
            super::feed(&mut view, key);
            let after = look(&mut view);
            // Saying Enter does nothing is not doing something.
            after != before && !after.contains(HINT)
        };
        if listed && !bound {
            inert.push(said(&key));
        } else if bound && !in_help && help.is_some() {
            unlisted.push(said(&key));
        }
    }
    (inert, unlisted)
}

/// [`disagree_by`] what the view draws, says and is about to do.
#[cfg(test)]
pub fn disagree<V: View>(make: impl FnMut() -> V) -> (Vec<String>, Vec<String>) {
    disagree_by(make, seen, "")
}

/// A key as a test says it: its code, and `+Ctrl` for a chord.
#[cfg(test)]
fn said(key: &KeyEvent) -> String {
    let ctrl = key.modifiers.contains(KeyModifiers::CONTROL);
    format!("{:?}{}", key.code, if ctrl { "+Ctrl" } else { "" })
}

/// Assert the view `make` builds lists every key it takes as it is then,
/// and no key it does not (see [`disagree`]); `what` names the state.
#[cfg(test)]
#[track_caller]
pub fn assert_lists_its_keys<V: View>(what: &str, make: impl FnMut() -> V) {
    agree(what, "", disagree(make));
}

/// [`assert_lists_its_keys`] on a view at an end of a list too short to
/// be anywhere else: the keys `end` names (as a pop-up writes them) may
/// do nothing there.
#[cfg(test)]
#[track_caller]
pub fn assert_lists_its_keys_at_an_end<V: View>(what: &str, end: &str, make: impl FnMut() -> V) {
    agree(what, end, disagree(make));
}

/// Fail on what [`disagree`] found in state `what`, but the keys `end`
/// names doing nothing.
#[cfg(test)]
#[track_caller]
pub fn agree(what: &str, end: &str, (mut inert, unlisted): (Vec<String>, Vec<String>)) {
    let ends: Vec<String> = named(end).iter().map(said).collect();
    inert.retain(|k| !ends.contains(k));
    assert!(
        inert.is_empty() && unlisted.is_empty(),
        "{what}: listed but doing nothing {inert:?}; bound but not listed {unlisted:?}"
    );
}

#[cfg(test)]
#[path = "../tests/tui/help.rs"]
mod tests;
