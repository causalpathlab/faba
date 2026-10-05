//! Pieces the terminal views share: centred pop-ups, the scroll window of a
//! list that keeps its selection in view, the go keys, and the rules for an
//! output folder.

use std::path::{Path, PathBuf};

use data_beans::interactive::ui::panel;
use ratatui::crossterm::event::{
    KeyCode, KeyEvent, KeyModifiers, KeyboardEnhancementFlags, PopKeyboardEnhancementFlags,
    PushKeyboardEnhancementFlags,
};
use ratatui::layout::Rect;
use ratatui::text::Line;
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

/// `lines` in a pop-up sized to them, at least 40 wide, centred in `area`.
pub fn popup(frame: &mut Frame, area: Rect, title: &str, lines: Vec<Line<'static>>) {
    let width = lines.iter().map(Line::width).max().unwrap_or(0) as u16 + 4;
    let inner = popup_frame(
        frame,
        area,
        width.max(40),
        lines.len() as u16 + 2,
        title.into(),
    );
    frame.render_widget(Paragraph::new(lines), inner);
}

/// The first of `rows` visible items of a `len`-long list that keeps item
/// `at` near the middle.
pub fn first_visible(at: usize, len: usize, rows: usize) -> usize {
    at.saturating_sub(rows.saturating_sub(1) / 2)
        .min(len.saturating_sub(rows))
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

/// The go keys, as footers name them.
pub const GO_KEYS: &str = "⇧Enter/G";

/// Whether `key` asks to go ahead: Shift+Enter, or `G` where the terminal
/// cannot tell Shift+Enter from Enter.
pub fn is_go(key: &KeyEvent) -> bool {
    match key.code {
        KeyCode::Enter => key.modifiers.contains(KeyModifiers::SHIFT),
        KeyCode::Char('G') => true,
        _ => false,
    }
}

/// Asks the terminal to report Shift+Enter (the kitty keyboard protocol;
/// others ignore the request) on the screen the view draws on, and takes the
/// request back when the view ends.
#[derive(Default)]
pub struct ShiftEnter {
    /// Ask at the next draw.
    pub want: bool,
    /// The request is in force.
    pub on: bool,
}

impl ShiftEnter {
    pub fn wanted() -> Self {
        ShiftEnter {
            want: true,
            on: false,
        }
    }

    /// Call at the top of `render`: terminals keep the main and alternate
    /// screens' keyboard modes apart, so ask once the view's screen is up.
    pub fn arm(&mut self) {
        if std::mem::take(&mut self.want) {
            let flags = KeyboardEnhancementFlags::DISAMBIGUATE_ESCAPE_CODES;
            self.on = ratatui::crossterm::execute!(
                std::io::stdout(),
                PushKeyboardEnhancementFlags(flags)
            )
            .is_ok();
        }
    }

    /// Call when the view ends.
    pub fn release(&mut self) {
        if std::mem::take(&mut self.on) {
            let _ = ratatui::crossterm::execute!(std::io::stdout(), PopKeyboardEnhancementFlags);
        }
    }
}

#[cfg(test)]
#[path = "tests/tui.rs"]
mod tests;
