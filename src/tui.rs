//! Pieces the terminal views share: centred pop-ups, lists that keep their
//! selection in view and narrow as you type, the key rules every view is fed
//! by (the apply key above all), typed paths, and the rules for an output
//! folder.

pub mod browser;

use std::path::{Path, PathBuf};

use data_beans::interactive::ui::{panel, Screen, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseEvent};
use ratatui::layout::Rect;
use ratatui::text::{Line, Span};
use ratatui::widgets::{Clear, Paragraph};
use ratatui::Frame;

/// A `w` x `h` rectangle centred in `area`, shrunk to fit it.
pub fn centered(area: Rect, w: u16, h: u16) -> Rect {
    let (w, h) = (w.min(area.width), h.min(area.height));
    Rect::new(
        area.x + (area.width - w) / 2,
        area.y + (area.height - h) / 2,
        w,
        h,
    )
}

/// Clear a `w` x `h` pop-up centred in `area`, frame it in a panel titled
/// `title`, and return the inside.
pub fn popup_frame(frame: &mut Frame, area: Rect, w: u16, h: u16, title: String) -> Rect {
    let rect = centered(area, w, h);
    frame.render_widget(Clear, rect);
    let block = panel(title, true);
    let inner = block.inner(rect);
    frame.render_widget(block, rect);
    inner
}

/// A row of buttons to click, `[ label ]` each, two apart: the line, and
/// each button's column (from the line's start) and width.
pub fn buttons(labels: &[&str]) -> (Line<'static>, Vec<(u16, u16)>) {
    let style = HIGHLIGHT.add_modifier(ratatui::style::Modifier::REVERSED);
    let mut spans = vec![Span::raw(" ")];
    let mut spots = Vec::new();
    let mut x = 1;
    for (k, label) in labels.iter().enumerate() {
        if k > 0 {
            spans.push(Span::raw("  "));
            x += 2;
        }
        let text = format!("[ {label} ]");
        let w = text.chars().count() as u16;
        spans.push(Span::styled(text, style));
        spots.push((x, w));
        x += w;
    }
    (Line::from(spans), spots)
}

/// `lines` in a pop-up sized to them, at least 40 wide, centred in `area`;
/// the inside, where they went.
pub fn popup(frame: &mut Frame, area: Rect, title: &str, lines: Vec<Line<'static>>) -> Rect {
    let width = lines.iter().map(Line::width).max().unwrap_or(0) as u16 + 4;
    let inner = popup_frame(
        frame,
        area,
        width.max(40),
        lines.len() as u16 + 2,
        title.into(),
    );
    frame.render_widget(Paragraph::new(lines), inner);
    inner
}

/// A pop-up with `buttons` on its first line (where a long body cannot push
/// them out of sight), then `body`; each button's spot goes into `hits`.
pub fn button_popup<T: Clone>(
    frame: &mut Frame,
    area: Rect,
    title: &str,
    buttons_: &[(&str, T)],
    body: Vec<Line<'static>>,
    hits: &Hits<T>,
) -> Rect {
    let labels: Vec<&str> = buttons_.iter().map(|(l, _)| *l).collect();
    let (line, spots) = buttons(&labels);
    let mut lines = vec![line, Line::raw("")];
    lines.extend(body);
    let inner = popup(frame, area, title, lines);
    for (&(dx, w), (_, what)) in spots.iter().zip(buttons_) {
        hits.add(Rect::new(inner.x + dx, inner.y, w, 1), what.clone());
    }
    inner
}

/// A tab line drawn in `area`: the names, `active` lit, then `button` if
/// any; each name's and the button's spot goes into `hits`.
pub fn tab_bar<T: Clone>(
    area: Rect,
    names: &[&str],
    active: usize,
    button: Option<(&str, T)>,
    hits: &Hits<T>,
    tab: impl Fn(usize) -> T,
) -> Line<'static> {
    let mut spans = vec![Span::raw(" ")];
    let mut x = area.x + 1;
    for (i, name) in names.iter().enumerate() {
        let style = if i == active {
            HIGHLIGHT
        } else {
            data_beans::interactive::ui::DIM
        };
        let w = name.chars().count() as u16;
        hits.add(Rect::new(x, area.y, w, 1), tab(i));
        spans.push(Span::styled(format!("{name}  "), style));
        x += w + 2;
    }
    if let Some((label, what)) = button {
        let (line, spots) = buttons(&[label]);
        if let Some(&(dx, w)) = spots.first() {
            hits.add(Rect::new(x + dx, area.y, w, 1), what);
        }
        spans.extend(line.spans);
    }
    Line::from(spans)
}

/// How many of `width` cells a bar shows filled for `done` of `total`.
pub fn filled(done: u64, total: u64, width: usize) -> usize {
    match total {
        0 => 0,
        _ => (done.min(total) as f64 / total as f64 * width as f64).round() as usize,
    }
}

/// The x axis of a text plot in `area` without its tick marks: the axis
/// line runs through where the plot put a `┴` under each label. For labels
/// that name a span (a metagene's regions), where a tick would read as a
/// position.
pub fn plain_axis(buf: &mut ratatui::buffer::Buffer, area: Rect) {
    for y in area.top()..area.bottom() {
        for x in area.left()..area.right() {
            if buf[(x, y)].symbol() == "┴" {
                buf[(x, y)].set_symbol("─");
            }
        }
    }
}

/// The marker of a list's highlighted row.
pub fn marker(on: bool) -> Span<'static> {
    Span::styled(if on { "▸ " } else { "  " }, HIGHLIGHT)
}

/// The `rows` lines of `lines` around line `at`.
pub fn window(lines: Vec<Line<'static>>, at: usize, rows: usize) -> Vec<Line<'static>> {
    let first = first_visible(at, lines.len(), rows);
    lines.into_iter().skip(first).take(rows).collect()
}

/// The first of `rows` visible items of a `len`-long list that keeps item
/// `at` near the middle.
pub fn first_visible(at: usize, len: usize, rows: usize) -> usize {
    at.saturating_sub(rows.saturating_sub(1) / 2)
        .min(len.saturating_sub(rows))
}

/// `at` moved by `d`, kept within `0..=last`.
pub fn moved(at: usize, d: isize, last: usize) -> usize {
    at.saturating_add_signed(d).min(last)
}

/// The home folder, if known.
pub fn home() -> Option<PathBuf> {
    std::env::var_os("HOME").map(PathBuf::from)
}

/// `p` with the home folder as `~`.
pub fn tilde(p: &Path) -> String {
    if let Some(rest) = home().and_then(|h| p.strip_prefix(h).ok().map(Path::to_path_buf)) {
        return if rest.as_os_str().is_empty() {
            "~".into()
        } else {
            format!("~/{}", rest.display())
        };
    }
    p.display().to_string()
}

/// A path as typed into a line: `~` is home, and the rest is made absolute
/// with `.` and `..` resolved (see [`browser::normalize`]).
pub fn typed_path(text: &str) -> PathBuf {
    let path = match (text.strip_prefix('~'), home()) {
        (Some(rest), Some(h)) if rest.is_empty() || rest.starts_with('/') => {
            h.join(rest.trim_start_matches('/'))
        }
        _ => PathBuf::from(text),
    };
    browser::normalize(&path)
}

/// `base/stem`, or `base/stem2`, `base/stem3`, ... : the first that does
/// not exist yet.
pub(crate) fn next_free(base: &Path, stem: &str) -> PathBuf {
    let mut out = base.join(stem);
    let mut n = 2;
    while out.exists() {
        out = base.join(format!("{stem}{n}"));
        n += 1;
    }
    out
}

/// Why `out` cannot take a run's outputs, if it cannot: it is unnamed, a
/// file, or a folder with files in it.
pub(crate) fn output_problem(out: &str) -> Option<String> {
    let path = Path::new(out);
    if out.is_empty() {
        return Some("name a directory".into());
    }
    if path.is_file() {
        return Some(format!("{out} is a file"));
    }
    let non_empty = std::fs::read_dir(path).is_ok_and(|mut d| d.next().is_some());
    non_empty.then(|| format!("{out} already contains files; choose an empty one"))
}

/// The apply key, as footers name it: Ctrl+R, a control byte every
/// terminal sends as it is, so it is never taken for Enter.
pub const APPLY_KEYS: &str = "Ctrl+R";

/// What a view says when Enter (plain or held with Shift, Alt or Ctrl) did
/// nothing: the key that `does` it.
pub fn enter_hint(does: &str) -> String {
    format!("Enter does nothing here; {APPLY_KEYS} {does}")
}

/// A full-screen view fed by the shared key rules (see [`feed`]); the mouse,
/// notices and which keys the terminal reports are data_beans' `Screen`'s.
pub trait View: Screen {
    /// Enter held with Shift, Alt or Ctrl came, which must not act as
    /// Enter: say what the apply key does here.
    fn stray_enter(&mut self) {}

    /// Whether the view has something for the apply key to do; elsewhere it
    /// is dropped, so it never acts as a letter key.
    fn takes_apply(&self) -> bool {
        false
    }
}

/// The ↑/↓ a turn of the wheel stands for.
pub fn wheel_key(m: &MouseEvent) -> Option<KeyEvent> {
    use ratatui::crossterm::event::MouseEventKind;
    let code = match m.kind {
        MouseEventKind::ScrollUp => KeyCode::Up,
        MouseEventKind::ScrollDown => KeyCode::Down,
        _ => return None,
    };
    Some(KeyEvent::new(code, KeyModifiers::NONE))
}

/// Where a view drew what, for the mouse: rectangles, each with what a
/// click there means. Kept behind a `RefCell` so drawing, which only
/// reads the view, can write it.
pub struct Hits<T> {
    spots: std::cell::RefCell<Vec<(Rect, T)>>,
}

impl<T> Default for Hits<T> {
    fn default() -> Self {
        Self {
            spots: std::cell::RefCell::default(),
        }
    }
}

impl<T: Clone> Hits<T> {
    /// Forget the last draw's; call at the top of a draw.
    pub fn clear(&self) {
        self.spots.borrow_mut().clear();
    }

    pub fn add(&self, area: Rect, what: T) {
        self.spots.borrow_mut().push((area, what));
    }

    /// A list drawn in `area` from item `first`: `what(i)` for item `i` on
    /// each of its rows, up to `len` items.
    pub fn rows(&self, area: Rect, first: usize, len: usize, what: impl Fn(usize) -> T) {
        for (k, i) in (first..len).take(area.height as usize).enumerate() {
            let row = Rect::new(area.x, area.y + k as u16, area.width, 1);
            self.add(row, what(i));
        }
    }

    /// A cell where `what` was drawn, to click in tests.
    #[cfg(test)]
    pub fn spot(&self, what: &T) -> Option<(u16, u16)>
    where
        T: PartialEq,
    {
        let spots = self.spots.borrow();
        let (r, _) = spots.iter().rev().find(|(_, t)| t == what)?;
        Some((r.x, r.y))
    }

    /// What is at a cell: the last drawn there, which is on top.
    pub fn at(&self, column: u16, row: u16) -> Option<T> {
        let pos = ratatui::layout::Position::new(column, row);
        self.spots
            .borrow()
            .iter()
            .rev()
            .find(|(r, _)| r.contains(pos))
            .map(|(_, t)| t.clone())
    }
}

/// Hand `key` to `view` by the rules every view shares: Enter held with
/// Shift, Alt or Ctrl goes to [`View::stray_enter`], other Ctrl and Alt
/// chords are dropped, so Ctrl+S is never `s`, and so is the apply key in a
/// view with nothing for it to do. Views and the widgets in them can take
/// what comes as meant.
pub fn feed(view: &mut impl View, key: KeyEvent) {
    if is_stray(&key) {
        if key.code == KeyCode::Enter {
            view.stray_enter();
        }
    } else if !is_apply(&key) || view.takes_apply() {
        view.handle_key(key);
    }
}

/// Run `view` full screen until it is done (data_beans' `run_screen`), its
/// keys fed by [`feed`].
pub fn run_view<V: View>(view: &mut V) -> anyhow::Result<()> {
    data_beans::interactive::ui::run_screen_with(view, |v, key| feed(v, key))
}

/// Whether a mouse event is one the views act on: a left click or a turn of
/// the wheel.
pub fn is_click_or_wheel(m: &MouseEvent) -> bool {
    use ratatui::crossterm::event::{MouseButton, MouseEventKind};
    matches!(
        m.kind,
        MouseEventKind::Down(MouseButton::Left)
            | MouseEventKind::ScrollUp
            | MouseEventKind::ScrollDown
    )
}

/// Whether `key` is the apply key: Ctrl+R (with Shift or not, which
/// terminals do not tell apart), and nothing else, so a stray key cannot
/// write anything.
pub fn is_apply(key: &KeyEvent) -> bool {
    key.modifiers.contains(KeyModifiers::CONTROL) && matches!(key.code, KeyCode::Char('r' | 'R'))
}

/// The apply key, as a button sends it.
pub fn apply_key() -> KeyEvent {
    KeyEvent::new(KeyCode::Char('r'), KeyModifiers::CONTROL)
}

/// Whether `key`, as the terminal sent it, is held with a modifier that
/// makes it another key (Ctrl, Alt, or Shift on Enter) and is not the apply
/// key. Views never see such keys, so Shift+Enter, Alt+Enter or Ctrl+Enter
/// never acts as a plain Enter, nor Ctrl+S as `s`.
fn is_stray(key: &KeyEvent) -> bool {
    let chord = KeyModifiers::CONTROL
        | KeyModifiers::ALT
        | KeyModifiers::SUPER
        | KeyModifiers::META
        | KeyModifiers::HYPER;
    let held = key.modifiers.intersects(chord)
        || (key.code == KeyCode::Enter && key.modifiers.contains(KeyModifiers::SHIFT));
    held && !is_apply(key)
}

/// What has been typed into a file list to narrow it: letters, digits and
/// punctuation type, so a file pane takes no letter keys of its own.
#[derive(Default, Debug, Clone)]
pub struct Find {
    pub text: String,
}

impl Find {
    /// Take `key` if it edits the text: a typed character (not Space, which
    /// panes keep for choosing), or Backspace and Esc while there is text.
    pub fn key(&mut self, key: &KeyEvent) -> bool {
        match key.code {
            KeyCode::Char(c) if c != ' ' && !is_apply(key) => self.text.push(c),
            KeyCode::Backspace if !self.text.is_empty() => {
                self.text.pop();
            }
            KeyCode::Esc if !self.text.is_empty() => self.text.clear(),
            _ => return false,
        }
        true
    }

    /// Where `key` jumps when typing a path: `/` with nothing typed to the
    /// root, `~` to home, and `/` after a name to that folder in `cwd` (the
    /// name as typed, else `dir`, the highlighted folder), so typing
    /// `/data/runs/` walks there.
    pub fn jump(&self, key: &KeyEvent, cwd: &Path, dir: Option<PathBuf>) -> Option<PathBuf> {
        let KeyCode::Char(c @ ('/' | '~')) = key.code else {
            return None;
        };
        match (c, self.text.as_str()) {
            ('/', "") => cwd.ancestors().last().map(Path::to_path_buf),
            ('~', "") => home(),
            ('/', "..") => cwd.parent().map(Path::to_path_buf),
            ('/', name) => Some(cwd.join(name)).filter(|p| p.is_dir()).or(dir),
            _ => None,
        }
    }

    /// A pane title's tail naming the text, if any.
    pub fn tag(&self) -> String {
        if self.text.is_empty() {
            String::new()
        } else {
            format!("· find: {} ", self.text)
        }
    }
}

/// What a narrowed list says when nothing matches the typed text.
pub const NOTHING_MATCHES: &str = "nothing matches; Backspace or Esc";

/// A list of named items narrowed by what has been typed, with a cursor
/// over those shown: a folder's listing, a site's releases or species.
pub struct FindList<T> {
    all: Vec<T>,
    /// `all`'s names in lower case, to match against.
    lower: Vec<String>,
    /// Indices into `all` of the items shown.
    shown: Vec<usize>,
    /// The cursor, among the items shown.
    pub at: usize,
    pub find: Find,
}

impl<T> Default for FindList<T> {
    fn default() -> Self {
        Self {
            all: Vec::new(),
            lower: Vec::new(),
            shown: Vec::new(),
            at: 0,
            find: Find::default(),
        }
    }
}

impl<T: AsRef<str>> FindList<T> {
    /// Show all of `all`, with nothing typed and the cursor at `at`.
    pub fn set(&mut self, all: Vec<T>, at: usize) {
        self.lower = all.iter().map(|x| x.as_ref().to_lowercase()).collect();
        self.shown = (0..all.len()).collect();
        self.all = all;
        self.at = at.min(self.shown.len().saturating_sub(1));
        self.find = Find::default();
    }

    /// The items shown, in order.
    pub fn shown(&self) -> impl Iterator<Item = &T> {
        self.shown.iter().map(|&i| &self.all[i])
    }

    pub fn len(&self) -> usize {
        self.shown.len()
    }

    pub fn is_empty(&self) -> bool {
        self.shown.is_empty()
    }

    /// The item under the cursor.
    pub fn highlighted(&self) -> Option<&T> {
        self.shown.get(self.at).map(|&i| &self.all[i])
    }

    pub fn step(&mut self, d: isize) {
        self.at = moved(self.at, d, self.len().saturating_sub(1));
    }

    /// Put the cursor on the item shown named `name`; whether there is one.
    pub fn select(&mut self, name: &str) -> bool {
        let at = self.shown().position(|x| x.as_ref() == name);
        self.at = at.unwrap_or(self.at);
        at.is_some()
    }

    /// Take `key` if it edits the typed text (see [`Find::key`]) or moves
    /// the cursor (↑/↓, PgUp/PgDn; letters type, so not `j`/`k`).
    pub fn key(&mut self, key: &KeyEvent) -> bool {
        if self.find.key(key) {
            self.narrow();
            return true;
        }
        let d = match key.code {
            KeyCode::Up => -1,
            KeyCode::Down => 1,
            KeyCode::PageUp => -10,
            KeyCode::PageDown => 10,
            _ => return false,
        };
        self.step(d);
        true
    }

    /// Show the items holding the typed text, any case, never `..`; the
    /// cursor on the first that starts with it, else the first.
    fn narrow(&mut self) {
        let text = self.find.text.to_lowercase();
        let lower = &self.lower;
        self.shown = (0..lower.len())
            .filter(|&i| text.is_empty() || (lower[i] != ".." && lower[i].contains(&text)))
            .collect();
        self.at = self
            .shown
            .iter()
            .position(|&i| lower[i].starts_with(&text))
            .unwrap_or(0);
    }

    /// What an empty list says: `empty`, or that nothing matches.
    pub fn empty_note<'a>(&self, empty: &'a str) -> &'a str {
        if self.find.text.is_empty() {
            empty
        } else {
            NOTHING_MATCHES
        }
    }
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
