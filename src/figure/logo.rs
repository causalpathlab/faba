//! A two-row faba mark for the top-right corner of the full-screen views.

use ratatui::buffer::Buffer;
use ratatui::layout::Rect;
use ratatui::style::{Color, Modifier, Style};

/// One bean in quadrant blocks, five columns by two rows; the notches in
/// its top row are its eyes.
pub(crate) const BEAN: [&str; 2] = ["▟▛█▜▙", "▝███▘"];

/// Four beans, as in the `--help` logo, in two dried-fava tones.
pub(crate) const TONES: [Color; 4] = [
    Color::Rgb(196, 152, 84),
    Color::Rgb(150, 105, 55),
    Color::Rgb(196, 152, 84),
    Color::Rgb(150, 105, 55),
];

/// Columns the mark takes: four beans a column apart, a space, the name.
const WIDTH: u16 = 4 * 5 + 3 + 1 + 4;

/// Draw the mark at the right end of the first two rows of `area`, left of
/// the session keys when they end the title row (see
/// [`crate::tui::title_keys_width`]), when it leaves `reserve` columns free
/// on its left for the rows' own text.
pub fn draw_mini_logo(buf: &mut Buffer, area: Rect, reserve: u16) {
    draw_beside(buf, area, reserve, crate::tui::title_keys_width());
}

/// [`draw_mini_logo`], `keys` columns at the right end of the top row taken.
fn draw_beside(buf: &mut Buffer, area: Rect, reserve: u16, keys: u16) {
    let right = WIDTH + 1 + keys;
    if area.height < 2 || area.width < right.saturating_add(reserve) {
        return;
    }
    let x = area.right() - right;
    for (b, &tone) in TONES.iter().enumerate() {
        let style = Style::default().fg(tone);
        for (dy, row) in BEAN.iter().enumerate() {
            buf.set_string(x + 6 * b as u16, area.y + dy as u16, row, style);
        }
    }
    let name = Style::default().fg(TONES[0]).add_modifier(Modifier::BOLD);
    buf.set_string(x + WIDTH - 4, area.y + 1, "faba", name);
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rows(width: u16) -> Vec<String> {
        let area = Rect::new(0, 0, width, 3);
        let mut buf = Buffer::empty(area);
        draw_mini_logo(&mut buf, area, 40);
        (0..3)
            .map(|y| {
                (0..width)
                    .map(|x| buf[(x, y)].symbol().to_string())
                    .collect()
            })
            .collect()
    }

    #[test]
    fn sits_in_the_top_right_corner_of_two_rows() {
        let r = rows(80);
        assert!(r[0].trim_end().ends_with("▟▛█▜▙ ▟▛█▜▙ ▟▛█▜▙ ▟▛█▜▙"));
        assert!(r[1].trim_end().ends_with("▝███▘ ▝███▘ ▝███▘ ▝███▘ faba"));
        assert!(r[2].trim().is_empty());
    }

    #[test]
    fn stays_out_of_a_narrow_screen() {
        assert!(rows(50).iter().all(|r| r.trim().is_empty()));
    }
}
