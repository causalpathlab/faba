//! The genes strip: the pinned genes side by side, in rows as many as fit
//! across, each on its own axis (its gene model's, else its sites' span)
//! under its model, an "all cells" track on a scale of its own above one
//! track per cell module on a shared one.

use super::bars::{mirror, peak};
use super::{ModulesView, Pane, Spot};
use crate::cluster::browse::data::bounds;
use crate::cluster::browse::Count;
use crate::site_analysis::miami::genemodel::GeneModel;
use crate::site_analysis::pileup::query_symbol;
use crate::site_analysis::pileup::tui::draw_gene_line;
use crate::tui::hidden_marks;
use data_beans::interactive::ui::{compact, ACCENTED, DIM, HIGHLIGHT, PLAIN};
use ratatui::buffer::Buffer;
use ratatui::layout::Rect;

/// Width of the strip's track labels, `C12  12.3k`.
const TRACK_LABEL: u16 = 12;

/// What the strip says before anything is pinned.
const EMPTY: &str =
    "Space on a site pins its gene here, or Alt-drag its row down; pins can come from any module";

/// How the strip's rows fall: a model row (when models are drawn), the
/// tracks, `lines` rows a side each, then the pointer and the names.
struct Rows {
    model: u16,
    lines: u16,
    tracks: usize,
}

impl Rows {
    fn top(&self, area: Rect, t: usize) -> u16 {
        area.y + self.model + t as u16 * 2 * self.lines
    }
}

impl ModulesView<'_> {
    /// The pinned genes in rows, as many as fit across each: the rows that
    /// fit the pane, scrolled to keep the focused gene's row in view.
    pub(super) fn draw_genes(&mut self, buf: &mut Buffer, area: Rect) {
        if self.pins.is_empty() {
            let len = EMPTY.chars().count() as u16;
            let x = area.x + area.width.saturating_sub(len) / 2;
            let y = area.y + area.height / 2;
            buf.set_stringn(x, y, EMPTY, area.width as usize, DIM);
            return;
        }
        let m = self.m;
        let widths: Vec<u16> = self
            .pins
            .iter()
            .map(|p| segment_width(m.genes(p.role)[p.gene].len()))
            .collect();
        let packed = pack(&widths, area.width.saturating_sub(TRACK_LABEL));
        let model = u16::from(!self.models.is_empty() || self.pending.is_some());
        let tracks = 1 + m.cell_roles - self.col_from;
        // A line a side per track, the pointer, the names and the rule.
        let least = model + 2 * tracks as u16 + 3;
        let fit = ((area.height / least) as usize).clamp(1, packed.len());
        let at = packed
            .iter()
            .position(|r| r.contains(&self.pin))
            .unwrap_or(0);
        self.row_from = self.row_from.min(at).max((at + 1).saturating_sub(fit));
        let height = area.height / fit as u16;
        let shown = packed.iter().skip(self.row_from).take(fit);
        for (k, pins) in shown.enumerate() {
            let y = area.y + k as u16 * height;
            // Each row but the last ends in a rule above the next.
            let last = k + 1 == fit;
            let h = if last { height } else { height - 1 };
            if !last {
                let rule = "─".repeat(area.width as usize);
                buf.set_stringn(area.x, y + h, rule, area.width as usize, DIM);
            }
            let row = Rect::new(area.x, y, area.width, h);
            self.draw_row(buf, row, pins.clone(), &widths, model);
        }
    }

    /// One row of the strip: pinned genes `pins` side by side in `area`.
    fn draw_row(
        &self,
        buf: &mut Buffer,
        area: Rect,
        pins: std::ops::Range<usize>,
        widths: &[u16],
        model: u16,
    ) {
        let rows = self.strip_rows(area, model);
        self.track_labels(buf, area, &rows);
        let mut x = area.x + TRACK_LABEL;
        for k in pins.clone() {
            if k > pins.start {
                for y in area.y..rows.top(area, rows.tracks) {
                    buf.set_string(x, y, "│", DIM);
                }
                x += 1;
            }
            self.draw_segment(buf, area, &rows, x, k, widths[k]);
            x += widths[k];
        }
    }

    /// Track `t`'s place among a column's reads: all cells first, then
    /// the cell modules from the first shown.
    fn track(&self, t: usize) -> usize {
        if t == 0 {
            0
        } else {
            t + self.col_from
        }
    }

    /// The rows one row of the strip has room for.
    fn strip_rows(&self, area: Rect, model: u16) -> Rows {
        let tracks = 1 + self.m.cell_roles - self.col_from;
        let avail = area.height.saturating_sub(model + 2);
        let lines = (avail / (2 * tracks as u16)).max(1);
        Rows {
            model,
            lines,
            tracks: ((avail / (2 * lines)) as usize).min(tracks),
        }
    }

    /// All cells, then each cell module shown, with its posterior size over
    /// the fitted cells; the cell modules hidden above and below counted
    /// under them.
    fn track_labels(&self, buf: &mut Buffer, area: Rect, rows: &Rows) {
        let m = self.m;
        let fitted: f64 = m.cell_sizes.iter().sum();
        for t in 0..rows.tracks {
            let y = rows.top(area, t) + rows.lines - 1;
            let (name, n, style) = match self.track(t) {
                0 => ("all".to_string(), fitted, ACCENTED),
                t => {
                    let c = m.columns[t - 1];
                    (format!("C{}", c + 1), m.cell_sizes[c], PLAIN)
                }
            };
            let label = format!("{name:<4}{:>6}", m.cells_label(n));
            buf.set_stringn(area.x, y, label, TRACK_LABEL as usize - 1, style);
        }
        let below = m.cell_roles - (self.col_from + rows.tracks.saturating_sub(1));
        let hidden = hidden_marks(self.col_from, below);
        let y = area.bottom().saturating_sub(2);
        buf.set_stringn(area.x, y, hidden, TRACK_LABEL as usize - 1, DIM);
    }

    /// Pin `k`'s gene in the `w` columns from `x`.
    fn draw_segment(&self, buf: &mut Buffer, area: Rect, rows: &Rows, x: u16, k: usize, w: u16) {
        let m = self.m;
        let pin = self.pins[k];
        let span = m.genes(pin.role)[pin.gene].clone();
        let gene = m.gene_of(pin.role, span.start);
        let members = &m.members[pin.role][span.clone()];
        let pos: Vec<i64> = members.iter().map(|s| m.sites[s.site].pos).collect();
        let model = self.models.get(gene);
        let extent = match model {
            Some(model) => (model.lo, model.hi - 1),
            None => bounds(pos.iter().copied()),
        };
        let focused = k == self.pin;
        let cursor = focused && self.pane == Pane::Genes;
        self.hits
            .add(Rect::new(x, area.y, w, area.height), Spot::Gene(k));
        // Per column: each track's reads (all cells first), and its first
        // member.
        let tracks = 1 + m.cell_roles;
        let mut cols = vec![vec![Count::default(); tracks]; w as usize];
        let mut first: Vec<Option<usize>> = vec![None; w as usize];
        let mut pointer = None;
        for (j, (member, col)) in members.iter().zip(columns(&pos, extent, w)).enumerate() {
            let col = col as usize;
            for (t, &r) in member.reads.iter().enumerate() {
                cols[col][0].add(r);
                cols[col][t + 1].add(r);
            }
            first[col].get_or_insert(span.start + j);
            if span.start + j == self.pin_site && focused {
                pointer = Some(col as u16);
            }
        }
        for (col, f) in first.iter().enumerate() {
            if let Some(i) = *f {
                let r = Rect::new(x + col as u16, area.y, 1, area.height);
                self.hits.add(r, Spot::Strip(k, i));
            }
        }
        if let (Some(model), true) = (model, rows.model > 0) {
            draw_model(buf, (x, area.y, w), model, extent, focused);
        }
        let all_peak = peak(cols.iter().map(|c| &c[0]));
        let cell_peak = peak(cols.iter().flat_map(|c| &c[1..]));
        for t in 0..rows.tracks {
            let top = if t == 0 { all_peak } else { cell_peak };
            for (j, col) in cols.iter().enumerate() {
                mirror(
                    buf,
                    x + j as u16,
                    rows.top(area, t),
                    rows.lines,
                    col[self.track(t)],
                    top,
                );
            }
        }
        if let (Some(p), true) = (pointer, cursor) {
            // The site's numbers beside its bars, over what is there.
            for t in 0..rows.tracks {
                let r = cols[p as usize][self.track(t)];
                let Some(f) = r.fraction() else { continue };
                let text = format!(" {:.0}% {} ", 100.0 * f, compact(r.total));
                let len = text.chars().count() as u16;
                let tx = if x + p + 1 + len <= area.right() {
                    x + p + 1
                } else {
                    (x + p).saturating_sub(len)
                };
                let style = if t == 0 { ACCENTED } else { PLAIN };
                buf.set_string(tx, rows.top(area, t) + rows.lines - 1, text, style);
            }
            buf.set_string(x + p, area.bottom() - 2, "▲", ACCENTED);
        }
        let style = if focused { HIGHLIGHT } else { DIM };
        let name = format!("{} M{}", query_symbol(gene), pin.role + 1);
        buf.set_stringn(x, area.bottom() - 1, name, w as usize - 1, style);
    }
}

/// A gene's model across a segment (`x`, `y`, `w`) over `extent`: a thin
/// arrowed intron line, heavy exons, in the accent for the gene in focus.
fn draw_model(
    buf: &mut Buffer,
    (x, y, w): (u16, u16, u16),
    g: &GeneModel,
    extent: (i64, i64),
    mine: bool,
) {
    let col = |p: i64| x + column(p, extent, w);
    let style = if mine { ACCENTED } else { PLAIN };
    draw_gene_line(buf, y, col, g, extent, style);
}

/// Columns a gene's segment takes in the strip: two per site and a margin,
/// within bounds that keep a name legible and a long gene from crowding out
/// the rest.
pub fn segment_width(sites: usize) -> u16 {
    (2 * sites + 2).clamp(10, 32) as u16
}

/// Each position's column within a segment `width` wide over `extent`.
pub fn columns(pos: &[i64], extent: (i64, i64), width: u16) -> Vec<u16> {
    pos.iter().map(|&p| column(p, extent, width)).collect()
}

/// Position `p`'s column within a segment `width` wide over `extent`: the
/// inner columns span it, positions outside it clipped to its ends.
fn column(p: i64, (lo, hi): (i64, i64), width: u16) -> u16 {
    if hi <= lo {
        return width / 2;
    }
    let last = width.saturating_sub(2).max(1);
    let t = (p - lo) as f64 / (hi - lo) as f64;
    1 + (t.clamp(0.0, 1.0) * (last - 1) as f64).round() as u16
}

/// Segments `widths` wide packed into rows `room` columns wide, in order,
/// each after the first in a row taking a divider: each row's span, at
/// least one segment a row.
pub fn pack(widths: &[u16], room: u16) -> Vec<std::ops::Range<usize>> {
    let mut rows = Vec::new();
    let (mut start, mut used) = (0, 0u16);
    for (k, &w) in widths.iter().enumerate() {
        let need = if k > start { w + 1 } else { w };
        if k > start && used + need > room {
            rows.push(start..k);
            (start, used) = (k, w);
        } else {
            used += need;
        }
    }
    if start < widths.len() {
        rows.push(start..widths.len());
    }
    rows
}
