//! The 2D view (`v`): the thresholds two at a time, in four panels in
//! place of the histogram, gene and metagene (see [`crate::qc::pairs`]).
//! Tab reaches the panels; there ←/→ pick one, `x`/`X` and `y`/`Y` step its
//! axes, and `A` lets Auto choose the pairs again. The thresholds still
//! move live from their table.

use super::*;
use crate::figure::style::Sheet;
use ratatui::style::Color;
use std::sync::Arc;

/// The 2D figure's own size, in points, when the page sets neither side.
const OWN_SIZE: (f64, f64) = (720.0, 640.0);

impl SitePicker<'_> {
    /// `v`: the 2D view in place of the histograms, or back.
    pub(super) fn toggle_two_d(&mut self) {
        self.two_d = !self.two_d;
        if !self.two_d && self.panel == Panel::Pairs {
            self.panel = Panel::Thresholds;
        }
        self.ensure_pairs();
    }

    /// The modality on screen's panels, Auto's when none were picked yet:
    /// picked only once the 2D view shows them.
    pub(super) fn ensure_pairs(&mut self) {
        if self.two_d && self.view().pairs.is_empty() {
            self.views[self.modality].auto_pairs();
        }
    }

    /// `A`: Auto's pairs for the modality on screen, at the thresholds now.
    pub(super) fn auto_pairs(&mut self) {
        self.views[self.modality].auto_pairs();
        self.pair_at = self.pair_at.min(self.view().pairs.len().saturating_sub(1));
    }

    /// Focus the panel `delta` on, around.
    pub(super) fn step_pair(&mut self, delta: isize) {
        let n = self.view().pairs.len().max(1) as isize;
        self.pair_at = (self.pair_at as isize + delta).rem_euclid(n) as usize;
    }

    /// Step the focused panel's x (`along_x`) or y axis by `delta` among the
    /// modality's axes, past the other axis; chosen by hand from then on.
    pub(super) fn step_axis(&mut self, along_x: bool, delta: isize) {
        let (m, k) = (self.modality, self.pair_at);
        let axes = pairs::axes(&self.views[m].sites());
        let Some(pair) = self.views[m].pairs.get_mut(k) else {
            return;
        };
        let (moved, other) = if along_x {
            (pair.0, pair.1)
        } else {
            (pair.1, pair.0)
        };
        let n = axes.len() as isize;
        let Some(at) = axes.iter().position(|&a| a == moved) else {
            return;
        };
        let mut next = at as isize;
        for _ in 0..n {
            next = (next + delta).rem_euclid(n);
            if axes[next as usize] != other {
                break;
            }
        }
        let to = axes[next as usize];
        if along_x {
            pair.0 = to;
        } else {
            pair.1 = to;
        }
        pair.2 = "chosen by hand".into();
    }

    /// The 2D view in `area`: the modality's panels two by two.
    pub(super) fn render_pairs(&mut self, frame: &mut Frame, area: Rect) {
        self.ensure_pairs();
        let [top, bottom] =
            Layout::vertical([Constraint::Fill(1), Constraint::Fill(1)]).areas(area);
        let cells: Vec<Rect> = [top, bottom]
            .into_iter()
            .flat_map(|row| {
                let [l, r] =
                    Layout::horizontal([Constraint::Fill(1), Constraint::Fill(1)]).areas(row);
                [l, r]
            })
            .collect();
        let n = self.view().pairs.len();
        for (k, &cell) in cells.iter().enumerate().take(n) {
            self.render_pair(frame, cell, k);
        }
    }

    /// Panel `k` in `area`: its squares, then the box's counts and why it
    /// is shown.
    fn render_pair(&mut self, frame: &mut Frame, area: Rect, k: usize) {
        let (x, y, why) = self.view().pairs[k].clone();
        let focused = self.panel == Panel::Pairs && self.pair_at == k;
        let block = panel(format!(" {} by {} ", y.short(), x.short()), focused);
        let inner = block.inner(area);
        frame.render_widget(block, area);
        let [plot, said] =
            Layout::vertical([Constraint::Min(3), Constraint::Length(2)]).areas(inner);
        // Binned once, at the size it is drawn: an image's squares, or one a
        // cell of the terminal's own.
        let images = self.controls.images().is_some();
        let bins = match images {
            true => pairs::BINS,
            false => (plot.width.max(plot.height) as usize).max(1),
        };
        let (key, g) = self.grid_of((x, y), bins);
        frame.render_widget(
            Paragraph::new(vec![
                Line::styled(g.in_box(true), DIM),
                Line::styled(why, DIM),
            ])
            .wrap(ratatui::widgets::Wrap { trim: true }),
            said,
        );
        if !images {
            draw_cells(frame, plot, &g, self.arm_look.front);
            return;
        }
        // Its image is drawn again only when what it shows changed.
        if self.pair_drawn[k] != Some(key) {
            self.pair_plots[k].invalidate();
            self.pair_drawn[k] = Some(key);
        }
        let mut image = std::mem::take(&mut self.pair_plots[k]);
        let how = pairs::PairLook {
            arms: self.arm_look,
            titles: Some(&self.titles),
        };
        let none = (Default::default(), 0.0);
        if let Some(picker) = self.controls.images() {
            image.render(frame, plot, picker, |w, h| {
                figure::svg(w, h, |c| g.draw_with(c, (0.0, 0.0, w, h), "", none, &how))
            });
        }
        self.pair_plots[k] = image;
    }

    /// The panel of `axes` binned `bins` a side, at the thresholds now:
    /// binned again only when the modality, its axes, the fail bits or the
    /// size changed since (see [`GridKey`]).
    fn grid_of(&mut self, (x, y): (Axis, Axis), bins: usize) -> (GridKey, Arc<pairs::Grid>) {
        let key = (self.modality, x, y, self.view().generation, bins);
        if let Some((_, g)) = self.pair_grids.iter().find(|(k, _)| *k == key) {
            return (key, g.clone());
        }
        let g = Arc::new(pairs::grid(
            &self.view().sites(),
            (x, y),
            &self.filter,
            bins,
        ));
        // Four panels, and room to step through their axes.
        if self.pair_grids.len() >= 16 {
            self.pair_grids.remove(0);
        }
        self.pair_grids.push((key, g.clone()));
        (key, g)
    }

    /// The 2D figure of modality `m` at the thresholds now, laid out for
    /// `sheet`: its own size grown to the font, a side asked for checked
    /// as a knob figure's is ([`super::knob_least`]).
    pub(super) fn pairs_figure(&self, m: usize, sheet: Sheet) -> String {
        let v = &self.views[m];
        let title = format!("{} sites, two thresholds at a time", v.table.modality);
        let size = sheet.fill_own(OWN_SIZE);
        let style = sheet.style(size);
        pairs::figure(
            &v.sites(),
            &v.shown_pairs(),
            &self.filter,
            size,
            (&title, self.arm_look),
            &style,
        )
    }
}

/// The panel drawn in the terminal's own cells, without images: a square
/// a cell, the arm in front first (see [`pairs::Grid::front_tone`]): the
/// signal arm's kept solid in its colour, its removed shaded, the control
/// arm's solid gray, each deeper as it holds more.
fn draw_cells(frame: &mut Frame, area: Rect, g: &pairs::Grid, front: crate::qc::arms::Front) {
    let (w, h) = (area.width as usize, area.height as usize);
    if w < 4 || h < 2 {
        return;
    }
    let buf = frame.buffer_mut();
    for row in 0..h {
        for col in 0..w {
            let bx = col * g.squares.bins / w;
            let by = (h - 1 - row) * g.squares.bins / h;
            let at = by * g.squares.bins + bx;
            let palette = crate::figure::style::Palette::Colour;
            let (glyph, (r, gr, b)) = match g.front_tone(at, palette, front) {
                Some(pairs::Tone::Fill(rgb)) => ('█', rgb),
                Some(pairs::Tone::Outlined(rgb)) => ('░', invert(rgb)),
                None => continue,
            };
            let cell = &mut buf[(area.x + col as u16, area.y + row as u16)];
            cell.set_char(glyph).set_fg(Color::Rgb(r, gr, b));
        }
    }
}

/// A dark outline's colour as a dark terminal shows it: as light.
fn invert((r, g, b): (u8, u8, u8)) -> (u8, u8, u8) {
    (255 - r, 255 - g, 255 - b)
}

/// The view's figures add each modality's 2D panels.
pub(super) fn pairs_stem(modality: &str) -> String {
    format!("{modality}_pairs")
}
