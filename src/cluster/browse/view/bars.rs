//! Mirror bars in eighths of a cell, as the pileup draws reads: converted
//! up in the accent, unconverted hanging below.

use super::super::Count;
use data_beans::interactive::ui::{ACCENTED, PLAIN};
use ratatui::buffer::Buffer;
use ratatui::style::{Modifier, Style};

/// Blocks from the bottom of a cell, in eighths.
const RISING: [&str; 9] = [" ", "▁", "▂", "▃", "▄", "▅", "▆", "▇", "█"];

/// The taller side of the tallest of `reads`, the top of their scale.
pub(super) fn peak<'r>(reads: impl IntoIterator<Item = &'r Count>) -> f64 {
    reads
        .into_iter()
        .map(|r| r.on.max(r.off()))
        .fold(0.0, f64::max)
}

/// A mirror bar in column `x` from row `y`: `lines` rows of converted reads
/// rising to the zero line, then `lines` of unconverted hanging from it, on
/// a scale where `peak` reads fill a side but its last eighth, so that rows
/// stacked one under another never touch.
pub(super) fn mirror(buf: &mut Buffer, x: u16, y: u16, lines: u16, reads: Count, peak: f64) {
    let (up, down) = (
        eighths(reads.on, peak, lines),
        eighths(reads.off(), peak, lines),
    );
    for k in 0..lines {
        let e = up.saturating_sub(8 * k as usize).min(8);
        if e > 0 {
            buf.set_string(x, y + lines - 1 - k, RISING[e], ACCENTED);
        }
        let e = down.saturating_sub(8 * k as usize).min(8);
        if e > 0 {
            let (glyph, style) = hanging(e);
            buf.set_string(x, y + lines + k, glyph, style);
        }
    }
}

/// Eighths of a side of `lines` rows that `value` fills, out of `peak`;
/// at least one for any reads.
pub fn eighths(value: f64, peak: f64, lines: u16) -> usize {
    if value <= 0.0 || peak <= 0.0 {
        return 0;
    }
    let most = 8 * lines as usize - 1;
    ((value / peak * most as f64).round() as usize).clamp(1, most)
}

/// A bar `e` eighths deep hanging from the top of a cell: a rising block
/// reversed, so its blank lower part shows the cell's ink above it.
pub fn hanging(e: usize) -> (&'static str, Style) {
    match e {
        8 => (RISING[8], PLAIN),
        e => (RISING[8 - e], PLAIN.add_modifier(Modifier::REVERSED)),
    }
}
