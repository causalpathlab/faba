//! Pictures of a drawn view, as the README shows them: a terminal buffer
//! drawn cell by cell as SVG (dark page, each cell's colours and weight,
//! lines and blocks drawn as shapes so they join), then as PNG.

use ratatui::buffer::Buffer;
use ratatui::style::{Color, Modifier};
use std::fmt::Write;

/// Where the pictures go: the folder `FABA_SCREENSHOTS_DIR` names, if any.
pub const DIR_VAR: &str = "FABA_SCREENSHOTS_DIR";

/// A cell's width and height, and the font's size, in points: cells a
/// whole number of pixels, so filled neighbours leave no seam.
const CELL_W: f64 = 4.0;
const CELL_H: f64 = 8.5;
const FONT: f64 = 6.6;
/// The margin around the cells, in points.
const PAD: f64 = 2.0;
/// Pixels per point in the PNG.
const SCALE: f32 = 2.0;

const PAGE_BG: (u8, u8, u8) = (0x1d, 0x20, 0x26);
const PAGE_FG: (u8, u8, u8) = (0xdc, 0xdf, 0xe4);
const FONTS: &str = "Menlo, 'DejaVu Sans Mono', monospace";

/// Write `buf` as `{name}.png` into the folder [`DIR_VAR`] names; nothing
/// when it names none.
pub fn save(name: &str, buf: &Buffer) {
    let Some(dir) = std::env::var_os(DIR_VAR) else {
        eprintln!("{DIR_VAR} is not set: {name}.png not written");
        return;
    };
    let dir = std::path::PathBuf::from(dir);
    std::fs::create_dir_all(&dir).unwrap();
    let png = crate::figure::png(&svg(buf), SCALE).unwrap();
    std::fs::write(dir.join(format!("{name}.png")), png).unwrap();
}

/// `c` as RGB, a terminal's default colours standing in for `Reset`.
fn rgb(c: Color, default: (u8, u8, u8)) -> (u8, u8, u8) {
    const ANSI: [(u8, u8, u8); 16] = [
        (0x1d, 0x20, 0x26),
        (0xe0, 0x6c, 0x75),
        (0x98, 0xc3, 0x79),
        (0xe5, 0xc0, 0x7b),
        (0x61, 0xaf, 0xef),
        (0xc6, 0x78, 0xdd),
        (0x56, 0xb6, 0xc2),
        (0xab, 0xb2, 0xbf),
        (0x5c, 0x63, 0x70),
        (0xf0, 0x8a, 0x92),
        (0xb5, 0xe0, 0x96),
        (0xf5, 0xd4, 0x95),
        (0x8c, 0xc8, 0xff),
        (0xdc, 0x9c, 0xf0),
        (0x7c, 0xd4, 0xde),
        (0xff, 0xff, 0xff),
    ];
    let cube = |v: u8| if v == 0 { 0 } else { 55 + 40 * v };
    match c {
        Color::Reset => default,
        Color::Rgb(r, g, b) => (r, g, b),
        Color::Indexed(i) if i < 16 => ANSI[usize::from(i)],
        Color::Indexed(i) if i < 232 => {
            let i = i - 16;
            (cube(i / 36), cube(i / 6 % 6), cube(i % 6))
        }
        Color::Indexed(i) => {
            let v = 8 + 10 * (i - 232);
            (v, v, v)
        }
        named => {
            let i = [
                Color::Black,
                Color::Red,
                Color::Green,
                Color::Yellow,
                Color::Blue,
                Color::Magenta,
                Color::Cyan,
                Color::Gray,
                Color::DarkGray,
                Color::LightRed,
                Color::LightGreen,
                Color::LightYellow,
                Color::LightBlue,
                Color::LightMagenta,
                Color::LightCyan,
                Color::White,
            ]
            .iter()
            .position(|&n| n == named)
            .unwrap_or(7);
            ANSI[i]
        }
    }
}

/// `a` moved `t` of the way to `b`.
fn mix(a: (u8, u8, u8), b: (u8, u8, u8), t: f64) -> (u8, u8, u8) {
    let m = |x: u8, y: u8| (f64::from(x) + (f64::from(y) - f64::from(x)) * t).round() as u8;
    (m(a.0, b.0), m(a.1, b.1), m(a.2, b.2))
}

fn hex((r, g, b): (u8, u8, u8)) -> String {
    format!("#{r:02x}{g:02x}{b:02x}")
}

/// The arms of a box-drawing line from a cell's middle: left, right, up,
/// down, each 0 (none), 1 (light) or 2 (heavy); and whether the corner is
/// rounded.
fn arms(c: char) -> Option<([u8; 4], bool)> {
    let a = match c {
        '─' => [1, 1, 0, 0],
        '━' => [2, 2, 0, 0],
        '│' => [0, 0, 1, 1],
        '┃' => [0, 0, 2, 2],
        '┌' => [0, 1, 0, 1],
        '┐' => [1, 0, 0, 1],
        '└' => [0, 1, 1, 0],
        '┘' => [1, 0, 1, 0],
        '┏' => [0, 2, 0, 2],
        '┓' => [2, 0, 0, 2],
        '┗' => [0, 2, 2, 0],
        '┛' => [2, 0, 2, 0],
        '├' => [0, 1, 1, 1],
        '┤' => [1, 0, 1, 1],
        '┬' => [1, 1, 0, 1],
        '┴' => [1, 1, 1, 0],
        '┼' => [1, 1, 1, 1],
        '╴' => [1, 0, 0, 0],
        '╶' => [0, 1, 0, 0],
        '╵' => [0, 0, 1, 0],
        '╷' => [0, 0, 0, 1],
        '╭' => return Some(([0, 1, 0, 1], true)),
        '╮' => return Some(([1, 0, 0, 1], true)),
        '╰' => return Some(([0, 1, 1, 0], true)),
        '╯' => return Some(([1, 0, 1, 0], true)),
        _ => return None,
    };
    Some((a, false))
}

/// A share of a cell: x, y, width and height, each of the cell's own.
type Part = (f64, f64, f64, f64);

/// The parts of a cell a block element fills, as (x, y, w, h) shares of
/// it, and the opacity it fills them at.
fn block(c: char) -> Option<(Vec<Part>, f64)> {
    let eighth = |n: u32| f64::from(n) / 8.0;
    let solid = |r: Vec<Part>| Some((r, 1.0));
    let q = |tl: bool, tr: bool, bl: bool, br: bool| {
        let mut r = Vec::new();
        for (on, x, y) in [
            (tl, 0.0, 0.0),
            (tr, 0.5, 0.0),
            (bl, 0.0, 0.5),
            (br, 0.5, 0.5),
        ] {
            if on {
                r.push((x, y, 0.5, 0.5));
            }
        }
        r
    };
    match c {
        '█' => solid(vec![(0.0, 0.0, 1.0, 1.0)]),
        '▀' => solid(vec![(0.0, 0.0, 1.0, 0.5)]),
        '▔' => solid(vec![(0.0, 0.0, 1.0, eighth(1))]),
        '▕' => solid(vec![(eighth(7), 0.0, eighth(1), 1.0)]),
        '▁'..='▇' => {
            let n = c as u32 - '▁' as u32 + 1;
            solid(vec![(0.0, 1.0 - eighth(n), 1.0, eighth(n))])
        }
        '▉'..='▏' => {
            let n = 8 - (c as u32 - '▉' as u32 + 1);
            solid(vec![(0.0, 0.0, eighth(n), 1.0)])
        }
        '▐' => solid(vec![(0.5, 0.0, 0.5, 1.0)]),
        '░' => Some((vec![(0.0, 0.0, 1.0, 1.0)], 0.25)),
        '▒' => Some((vec![(0.0, 0.0, 1.0, 1.0)], 0.5)),
        '▓' => Some((vec![(0.0, 0.0, 1.0, 1.0)], 0.75)),
        '▖' => solid(q(false, false, true, false)),
        '▗' => solid(q(false, false, false, true)),
        '▘' => solid(q(true, false, false, false)),
        '▙' => solid(q(true, false, true, true)),
        '▚' => solid(q(true, false, false, true)),
        '▛' => solid(q(true, true, true, false)),
        '▜' => solid(q(true, true, false, true)),
        '▝' => solid(q(false, true, false, false)),
        '▞' => solid(q(false, true, true, false)),
        '▟' => solid(q(false, true, true, true)),
        _ => None,
    }
}

/// `v` points to the nearest whole pixel.
fn snap(v: f64) -> f64 {
    (v * f64::from(SCALE)).round() / f64::from(SCALE)
}

/// `buf` as an SVG picture of a terminal showing it.
pub fn svg(buf: &Buffer) -> String {
    let area = buf.area;
    let (w, h) = (
        f64::from(area.width) * CELL_W + 2.0 * PAD,
        f64::from(area.height) * CELL_H + 2.0 * PAD,
    );
    let mut back = String::new();
    let mut fore = String::new();
    for y in 0..area.height {
        for x in 0..area.width {
            let cell = &buf[(area.x + x, area.y + y)];
            let (px, py) = (PAD + f64::from(x) * CELL_W, PAD + f64::from(y) * CELL_H);
            let m = cell.modifier;
            let (mut fg, mut bg) = (rgb(cell.fg, PAGE_FG), rgb(cell.bg, PAGE_BG));
            if m.contains(Modifier::REVERSED) {
                std::mem::swap(&mut fg, &mut bg);
            }
            if m.contains(Modifier::DIM) {
                fg = mix(fg, bg, 0.45);
            }
            if bg != PAGE_BG {
                let _ = write!(
                    back,
                    r#"<rect x="{px:.2}" y="{py:.2}" width="{:.2}" height="{:.2}" fill="{}"/>"#,
                    CELL_W,
                    CELL_H,
                    hex(bg)
                );
            }
            let sym = cell.symbol();
            let mut chars = sym.chars();
            let (Some(c), None) = (chars.next(), chars.next()) else {
                if !sym.is_empty() {
                    text(&mut fore, px, py, sym, fg, m);
                }
                continue;
            };
            if c == ' ' {
                continue;
            }
            let colour = hex(fg);
            if let Some((a, round)) = arms(c) {
                lines(&mut fore, px, py, a, round, &colour);
            } else if let Some((parts, alpha)) = block(c) {
                for (bx, by, bw, bh) in parts {
                    let _ = write!(
                        fore,
                        r#"<rect x="{:.2}" y="{:.2}" width="{:.2}" height="{:.2}" fill="{colour}" fill-opacity="{alpha}"/>"#,
                        px + bx * CELL_W,
                        py + snap(by * CELL_H),
                        bw * CELL_W,
                        snap(bh * CELL_H),
                    );
                }
            } else if ('\u{2801}'..='\u{28ff}').contains(&c) {
                braille(&mut fore, px, py, c as u32 - 0x2800, &colour);
            } else if c != '\u{2800}' {
                text(&mut fore, px, py, sym, fg, m);
            }
        }
    }
    format!(
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{w}" height="{h}" viewBox="0 0 {w} {h}"><rect width="{w}" height="{h}" fill="{}"/>{back}<g font-family="{FONTS}" font-size="{FONT}">{fore}</g></svg>"#,
        hex(PAGE_BG)
    )
}

/// One cell's text at its place.
fn text(out: &mut String, px: f64, py: f64, sym: &str, fg: (u8, u8, u8), m: Modifier) {
    let esc = sym
        .replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;");
    let weight = if m.contains(Modifier::BOLD) {
        r#" font-weight="bold""#
    } else {
        ""
    };
    let style = if m.contains(Modifier::ITALIC) {
        r#" font-style="italic""#
    } else {
        ""
    };
    let base = py + CELL_H / 2.0 + FONT * 0.36;
    let _ = write!(
        out,
        r#"<text x="{:.2}" y="{base:.2}" fill="{}"{weight}{style} xml:space="preserve">{esc}</text>"#,
        px,
        hex(fg)
    );
    if m.contains(Modifier::UNDERLINED) {
        let _ = write!(
            out,
            r#"<rect x="{px:.2}" y="{:.2}" width="{CELL_W}" height="0.4" fill="{}"/>"#,
            base + 1.0,
            hex(fg)
        );
    }
}

/// A box-drawing cell as strokes from its middle.
fn lines(out: &mut String, px: f64, py: f64, a: [u8; 4], round: bool, colour: &str) {
    // The middle, set so a light stroke covers whole pixels.
    let (cx, cy) = (px + CELL_W / 2.0 + 0.25, py + CELL_H / 2.0);
    let width = |n: u8| if n == 2 { 1.0 } else { 0.5 };
    if round {
        let hx = if a[0] > 0 { px } else { px + CELL_W };
        let vy = if a[2] > 0 { py } else { py + CELL_H };
        let r = CELL_W / 2.0;
        let (sx, sy) = (cx + (hx - cx).signum() * r, cy + (vy - cy).signum() * r);
        let _ = write!(
            out,
            r#"<path d="M{hx:.2},{cy:.2} L{sx:.2},{cy:.2} Q{cx:.2},{cy:.2} {cx:.2},{sy:.2} L{cx:.2},{vy:.2}" fill="none" stroke="{colour}" stroke-width="0.5"/>"#
        );
        return;
    }
    // Square caps carry each arm over the middle, so arms meet without a notch.
    let ends = [(px, cy), (px + CELL_W, cy), (cx, py), (cx, py + CELL_H)];
    for (n, (ex, ey)) in a.iter().zip(ends) {
        if *n > 0 {
            let _ = write!(
                out,
                r#"<line x1="{cx:.2}" y1="{cy:.2}" x2="{ex:.2}" y2="{ey:.2}" stroke="{colour}" stroke-width="{}" stroke-linecap="square"/>"#,
                width(*n)
            );
        }
    }
}

/// A braille cell's dots: two columns of four.
fn braille(out: &mut String, px: f64, py: f64, bits: u32, colour: &str) {
    const DOTS: [(u32, u32); 8] = [
        (0, 0),
        (0, 1),
        (0, 2),
        (1, 0),
        (1, 1),
        (1, 2),
        (0, 3),
        (1, 3),
    ];
    for (i, (col, row)) in DOTS.iter().enumerate() {
        if bits & (1 << i) != 0 {
            let _ = write!(
                out,
                r#"<circle cx="{:.2}" cy="{:.2}" r="0.55" fill="{colour}"/>"#,
                px + (f64::from(*col) + 0.5) * CELL_W / 2.0,
                py + (f64::from(*row) + 0.5) * CELL_H / 4.0,
            );
        }
    }
}
