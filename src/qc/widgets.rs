//! Pieces the `faba qc` views share: centred pop-ups and the scroll window of
//! a list that keeps its selection in view.

use data_beans::interactive::ui::panel;
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
