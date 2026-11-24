//! The 2D view (`d`): the thresholds two at a time, in four panels in
//! place of the histogram, gene and metagene (see [`crate::qc::pairs`]).
//! Tab reaches the panels; there ←/→ pick one, `x`/`X` and `y`/`Y` step its
//! axes, and `A` lets Auto choose the pairs again. The thresholds still
//! move live from their table.

use super::*;
use crate::figure::size::Page;
use ratatui::style::Color;

/// The 2D figure's own size, in points, when the page sets neither side.
const OWN_SIZE: (f64, f64) = (720.0, 640.0);

impl SitePicker<'_> {
    /// `d`: the 2D view in place of the histograms, or back.
    pub(super) fn toggle_two_d(&mut self) {
        self.two_d = !self.two_d;
        if !self.two_d && self.panel == Panel::Pairs {
            self.panel = Panel::Thresholds;
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
        let g = pairs::grid(&self.view().sites(), (x, y), &self.filter, pairs::BINS);
        frame.render_widget(
            Paragraph::new(vec![
                Line::styled(g.in_box_short(), DIM),
                Line::styled(why, DIM),
            ])
            .wrap(ratatui::widgets::Wrap { trim: true }),
            said,
        );
        let mut image = std::mem::take(&mut self.pair_plots[k]);
        let drawn = self.controls.images().is_some_and(|picker| {
            image.render(frame, plot, picker, |w, h| {
                figure::svg(w, h, |c| g.draw(c, (0.0, 0.0, w, h), ""))
            })
        });
        self.pair_plots[k] = image;
        if !drawn {
            draw_cells(frame, plot, &self.view().sites(), (x, y), &self.filter);
        }
    }

    /// The 2D figure of modality `m` at the thresholds now, laid out for
    /// `page`.
    pub(super) fn pairs_figure(&self, m: usize, page: Page) -> String {
        let v = &self.views[m];
        let title = format!("{} sites, two thresholds at a time", v.table.modality);
        let size = crate::figure::size::fill(page, || OWN_SIZE);
        pairs::figure(&v.sites(), &v.pairs, &self.filter, size, &title)
    }
}

/// The panel drawn in the terminal's own cells, without images: a square
/// a cell, kept or removed by its colour and deeper as it holds more.
fn draw_cells(frame: &mut Frame, area: Rect, s: &Sites, axes: (Axis, Axis), f: &SiteFilterArgs) {
    let (w, h) = (area.width as usize, area.height as usize);
    if w < 4 || h < 2 {
        return;
    }
    // A grid of the area's shape: as many squares as cells across.
    let g = pairs::grid(s, axes, f, w.max(h));
    let buf = frame.buffer_mut();
    for row in 0..h {
        for col in 0..w {
            let bx = col * g.bins / w;
            let by = (h - 1 - row) * g.bins / h;
            if let Some((r, gr, b)) = g.colour(by * g.bins + bx) {
                let cell = &mut buf[(area.x + col as u16, area.y + row as u16)];
                cell.set_char('█').set_fg(Color::Rgb(r, gr, b));
            }
        }
    }
}

/// The view's figures add each modality's 2D panels.
pub(super) fn pairs_stem(modality: &str) -> String {
    format!("{modality}_pairs")
}
