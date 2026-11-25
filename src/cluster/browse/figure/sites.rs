//! One gene module's sites as the browser's sites pane draws them under
//! its bars look: a row a site in the module's order, its gene's symbol and
//! position on the left (marked when its gene is pinned), a mirror bar in
//! each cell module (converted rising, unconverted hanging, on the row's
//! own scale) with its % beside it where there is room, then the site's
//! raw reads, converted, unconverted and cells over every cell read; the
//! rows that fit, the rest counted under them.

use super::{width, Ink, BLOCK_UP, BLOCK_UP_GRAY, MIRROR_LEAST, MIRROR_OFF};
use crate::cluster::browse::view::peak;
use crate::cluster::browse::Modules;
use crate::figure::{esc, Anchor, Canvas, FAINT, INK, MUTED};
use crate::site_analysis::pileup::{channel_names, query_symbol};
use crate::tui::ui::compact;

/// A site row's depth, in font sizes: the browser's two lines.
const ROW: f64 = 1.9;
/// The heads' line, in font sizes.
const HEAD: f64 = 1.5;
/// A cell module's column, in font sizes at most, and at least.
const CELL_MOST: f64 = 4.2;
const CELL_LEAST: f64 = 1.2;
/// The mark before a pinned gene's sites, as the browser's.
const PIN_MARK: &str = "• ";

/// The sites of gene role `g` drawn in `(x, y, w, h)` on `c`, the genes at
/// places `pinned` among its genes marked.
pub(super) fn draw(
    c: &mut Canvas,
    m: &Modules,
    g: usize,
    pinned: &[usize],
    (x, y, w, h): (f64, f64, f64, f64),
    ink: Ink,
) {
    let f = ink.f;
    let members = &m.members[g];
    let spans = m.genes(g);
    let gene_of = |i: usize| spans.iter().position(|s| s.contains(&i));
    let labels: Vec<String> = (0..members.len())
        .map(|i| {
            let s = &m.sites[members[i].site];
            let mark = match gene_of(i).is_some_and(|k| pinned.contains(&k)) {
                true => PIN_MARK,
                false => "",
            };
            format!("{mark}{} {}", query_symbol(&s.gene), s.pos)
        })
        .collect();
    let label_w = labels.iter().map(|l| width(l, f)).fold(0.0, f64::max) + f;
    // The counts' four columns, each as wide as `000.0k` and a gap.
    let field = width("000.0k", f) + 0.4 * f;
    let counts_w = 4.0 * field;
    let cols = m.columns.len().max(1);
    let cell = ((w - label_w - counts_w) / cols as f64).clamp(CELL_LEAST * f, CELL_MOST * f);
    let numbered = cell >= 3.2 * f;
    let (on, off) = channel_names(&m.modality);
    let initial = |s: &str| s.chars().next().unwrap_or('?').to_string();
    let heads = [
        "reads".to_string(),
        initial(on),
        initial(off),
        "cells".into(),
    ];
    let gx = x + label_w;
    let counts_x = gx + cell * cols as f64;
    let head_y = y + HEAD * f * 0.75;
    for (j, &role) in m.columns.iter().enumerate() {
        let name = format!("C{}", role + 1);
        c.text(gx + cell * j as f64, head_y, &name, f, Anchor::Start, MUTED);
    }
    for (k, head) in heads.iter().enumerate() {
        let at = counts_x + field * (k as f64 + 1.0);
        c.text(at, head_y, head, f, Anchor::End, MUTED);
    }
    let up_fill = match ink.gray {
        true => BLOCK_UP_GRAY,
        false => BLOCK_UP,
    };
    let row = ROW * f;
    let body = h - HEAD * f;
    let mut fit = ((body / row).floor() as usize).min(members.len());
    // The rows left out are counted on a line of their own.
    if fit < members.len() {
        fit = ((body - 1.2 * f) / row).floor().max(0.0) as usize;
    }
    for (i, member) in members.iter().take(fit).enumerate() {
        let top = y + HEAD * f + row * i as f64;
        let mid = top + row / 2.0;
        c.text(x, mid + f / 2.5, &esc(&labels[i]), f, Anchor::Start, INK);
        // The row's own scale, as the browser's bars are.
        let peak = peak(&member.reads).max(1e-12);
        let reach = row / 2.0 - 1.0;
        for (j, r) in member.reads.iter().enumerate() {
            let cx = gx + cell * j as f64;
            let bw = (cell * 2.0 / 6.0).min(0.9 * f);
            c.line(cx, mid, cx + bw, mid, FAINT, 0.3);
            let side = |reads: f64| match reads > 0.0 {
                true => (reach * reads / peak).max(MIRROR_LEAST * f),
                false => 0.0,
            };
            let (up, down) = (side(r.up()), side(r.off()));
            c.rect(cx, mid - up, bw, up, up_fill);
            c.rect(cx, mid, bw, down, MIRROR_OFF);
            if let (true, Some(n)) = (numbered, r.number()) {
                let tx = cx + bw + 0.25 * f;
                c.text(tx, mid + f / 2.5, &n, f * 0.9, Anchor::Start, MUTED);
            }
        }
        for (k, said) in counts(m, member.site).iter().enumerate() {
            let at = counts_x + field * (k as f64 + 1.0);
            c.text(at, mid + f / 2.5, said, f, Anchor::End, MUTED);
        }
    }
    if fit < members.len() {
        let more = members.len() - fit;
        let ty = y + HEAD * f + row * fit as f64 + f;
        let said = format!("+{more} more sites");
        c.text(x, ty, &said, f, Anchor::Start, MUTED);
    }
}

/// A site's raw reads over every cell read, converted, unconverted and
/// the cells with any, as the browser writes them; dashes without them.
fn counts(m: &Modules, site: usize) -> [String; 4] {
    match &m.stats {
        Some(s) => [
            compact((s.on[site] + s.off[site]) as f64),
            compact(s.on[site] as f64),
            compact(s.off[site] as f64),
            compact(s.cells[site] as f64),
        ],
        None => std::array::from_fn(|_| "-".to_string()),
    }
}
