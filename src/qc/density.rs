//! What faba qc's 2D figures share: things binned two ways into squares,
//! kept apart from the rest, a square's colour deeper as it holds more; and
//! the frame they are drawn in, its ticks and titles.

use crate::figure::style::{Palette, OWN_FONT};
use crate::figure::{mixed, rgb, Anchor, Canvas, FAINT, INK, MUTED};

/// A kept site or barcode, drawn in this.
pub const KEPT: &str = "#2a78d6";
/// One the cut removes.
pub const REMOVED: &str = "#eb6834";
/// Squares a side.
pub const BINS: usize = 44;

/// Things binned `bins` a side, kept and the rest per square, row-major
/// from the bottom left; and the busiest square's count.
#[derive(Clone, Debug, PartialEq, serde::Serialize)]
pub struct Squares {
    pub bins: usize,
    pub kept: Vec<u32>,
    pub other: Vec<u32>,
    pub top: u32,
}

impl Squares {
    /// From kept and the rest per square, `bins` a side.
    pub fn new(bins: usize, kept: Vec<u32>, other: Vec<u32>) -> Self {
        let top = (kept.iter().zip(&other))
            .map(|(k, o)| k + o)
            .max()
            .unwrap_or(0)
            .max(1);
        Self {
            bins,
            kept,
            other,
            top,
        }
    }

    /// Square `at`'s colour over white in `palette`'s kept and removed
    /// colours, by which most of it is, deeper as it holds more of the
    /// busiest; none when empty.
    pub fn tone(&self, at: usize, palette: Palette) -> Option<(u8, u8, u8)> {
        let (k, o) = (self.kept[at], self.other[at]);
        let n = k + o;
        (n > 0).then(|| {
            let depth = 0.18 + 0.82 * (n as f64).ln_1p() / (self.top as f64).ln_1p();
            let roles = palette.roles();
            let base = if k >= o { roles.kept } else { roles.removed };
            mixed((255.0, 255.0, 255.0), rgb(base), depth)
        })
    }

    /// Drawn over `plot`'s box, in the canvas' palette.
    pub fn draw(&self, c: &mut Canvas, plot: &Plot) {
        let palette = c.style().palette;
        let (cw, ch) = (plot.pw / self.bins as f64, plot.ph / self.bins as f64);
        for by in 0..self.bins {
            for bx in 0..self.bins {
                if let Some((r, g, b)) = self.tone(by * self.bins + bx, palette) {
                    let (x, y) = (
                        plot.px + bx as f64 * cw,
                        plot.py + plot.ph - (by + 1) as f64 * ch,
                    );
                    c.rect(
                        x,
                        y,
                        cw + 0.05,
                        ch + 0.05,
                        &format!("#{r:02x}{g:02x}{b:02x}"),
                    );
                }
            }
        }
    }
}

/// A plot's box in a figure and the ranges it spans.
pub struct Plot {
    pub px: f64,
    pub py: f64,
    pub pw: f64,
    pub ph: f64,
    pub xr: (f64, f64),
    pub yr: (f64, f64),
    /// The box `(x, y, w, h)` it sits in, its titles under and beside it.
    outer: (f64, f64, f64, f64),
}

/// Ticks' and titles' size, in points of the 8-point layout: as a bar
/// panel's, which the canvas scales to its font.
pub(super) const SIZE: f64 = OWN_FONT;

impl Plot {
    /// The plot in the box `(x, y, w, h)` over `xr` by `yr`, with room
    /// for ticks and titles at the layout's 8 points.
    pub fn on((x, y, w, h): (f64, f64, f64, f64), xr: (f64, f64), yr: (f64, f64)) -> Self {
        let (l, r, t, b) = (SIZE * 4.6, SIZE, SIZE * 0.6, SIZE * 3.6);
        Self {
            px: x + l,
            py: y + t,
            pw: (w - l - r).max(1.0),
            ph: (h - t - b).max(1.0),
            xr,
            yr,
            outer: (x, y, w, h),
        }
    }

    /// Where `(vx, vy)` falls.
    pub fn at(&self, vx: f64, vy: f64) -> (f64, f64) {
        let span = |(lo, hi): (f64, f64)| (hi - lo).max(1e-12);
        (
            self.px + (vx - self.xr.0) / span(self.xr) * self.pw,
            self.py + self.ph - (vy - self.yr.0) / span(self.yr) * self.ph,
        )
    }

    /// Its axes, ticks (drawn positions and labels) and titles. Every
    /// tick is marked; where their labels would run into each other, every
    /// second one is labelled (or third, ...), from the first.
    pub fn axes(&self, c: &mut Canvas, ticks: [&[(f64, String)]; 2], titles: [&str; 2]) {
        let (px, py, pw, ph, size) = (self.px, self.py, self.pw, self.ph, SIZE);
        c.line(px, py + ph, px + pw, py + ph, FAINT, 0.8);
        c.line(px, py, px, py + ph, FAINT, 0.8);
        let written = c.label_size(size);
        let xs: Vec<f64> = ticks[0]
            .iter()
            .map(|(v, _)| self.at(*v, self.yr.0).0)
            .collect();
        let wide: Vec<f64> = (ticks[0].iter())
            .map(|(_, l)| crate::figure::text_width(l, written))
            .collect();
        let every = thinned(&xs, &wide, 0.4 * size);
        for (k, ((_, label), &x)) in ticks[0].iter().zip(&xs).enumerate() {
            c.line(x, py + ph, x, py + ph + 3.0, FAINT, 0.8);
            if k % every == 0 {
                c.text(x, py + ph + size * 1.4, label, size, Anchor::Middle, MUTED);
            }
        }
        let ys: Vec<f64> = ticks[1]
            .iter()
            .map(|(v, _)| self.at(self.xr.0, *v).1)
            .collect();
        let every = thinned(&ys, &vec![written; ys.len()], 0.15 * size);
        for (k, ((_, label), &y)) in ticks[1].iter().zip(&ys).enumerate() {
            c.line(px - 3.0, y, px, y, FAINT, 0.8);
            if k % every != 0 {
                continue;
            }
            c.text(
                px - size * 0.5,
                y + size * 0.35,
                label,
                size,
                Anchor::End,
                MUTED,
            );
        }
        let (x, y, _, h) = self.outer;
        c.text(
            px + pw / 2.0,
            y + h - size * 0.6,
            titles[0],
            size,
            Anchor::Middle,
            INK,
        );
        c.vtext(x + size * 0.9, py + ph / 2.0, titles[1], size, INK);
    }
}

/// How many ticks apart labels go so none runs into the next: labels
/// centred at `at`, each `wide` along the axis, `gap` clear at least; 1
/// when every one fits, and as many as there are (the first alone) when
/// no two do.
pub(super) fn thinned(at: &[f64], wide: &[f64], gap: f64) -> usize {
    let clear = |every: usize| {
        let kept = || (0..at.len()).step_by(every);
        (kept().zip(kept().skip(1)))
            .all(|(a, b)| (at[b] - at[a]).abs() >= (wide[a] + wide[b]) / 2.0 + gap)
    };
    (1..=at.len().max(1)).find(|&e| clear(e)).unwrap_or(1)
}

/// Ticks at each decade from 10^ceil(lo) to 10^floor(hi) on a log10
/// axis, labelled as counts read.
pub fn decades((lo, hi): (f64, f64)) -> Vec<(f64, String)> {
    ((lo.ceil() as i32)..=(hi.floor() as i32))
        .map(|k| (k as f64, crate::tui::ui::compact(10f64.powi(k))))
        .collect()
}
