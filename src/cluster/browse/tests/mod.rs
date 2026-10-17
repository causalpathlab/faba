//! Tests of `faba modules`: a small set of block tables written to a
//! temporary folder, loaded and browsed. The fixtures are here; the tests
//! are by what they cover.

use super::*;
use crate::cluster::read::{Cell, Site, SiteStats};
use crate::cluster::write;
use data_beans::interactive::ui::Screen;
use ratatui::backend::TestBackend;
use ratatui::buffer::Buffer;
use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers, MouseEvent, MouseEventKind};
use ratatui::Terminal;

mod data;
mod gene;
mod report;
mod strip;
mod view;

fn site(gene: &str, pos: i64) -> Site {
    Site {
        gene: gene.into(),
        chr: "chr1".into(),
        pos,
    }
}

/// Two cell roles and three gene roles; gene role 3 has no sites. Gene
/// role 1's shifts differ by 0.2 across the cell roles, gene role 2's by 1.5.
///
/// Sites: GENE1 at 100 and 300, GENE2 at 50, GENE3 at 10. GENE1 is in gene
/// role 1 only; GENE2 has reads in roles 1 and 2 but role 1 is its main
/// one; GENE3 is in role 2.
///
/// Cells B0..B4: B0, B1 and B2 in cell role 1, B3 in role 2, B4 dropped;
/// B2 was placed, not fitted.
fn tables(dir: &Path) {
    tables_with(dir, false)
}

/// [`tables`], or with the two cell roles swapped, so the panes show them
/// in the other order.
fn tables_with(dir: &Path, swap: bool) {
    let cell = |c: usize| if swap { 1 - c } else { c };
    let sites = vec![
        site("GENE1", 100),
        site("GENE1", 300),
        site("GENE2", 50),
        site("GENE3", 10),
    ];
    // Blocks `c * 3 + g`; site-major `s * 6 + k`.
    let (blocks, roles) = (6, 3);
    let mut on = vec![0.0; sites.len() * blocks];
    let mut off = vec![0.0; sites.len() * blocks];
    let mut put = |s: usize, c: usize, g: usize, m: f64, u: f64| {
        let c = cell(c);
        on[s * blocks + c * roles + g] = m;
        off[s * blocks + c * roles + g] = u;
    };
    put(0, 0, 0, 8.0, 2.0);
    put(0, 1, 0, 1.0, 9.0);
    put(1, 0, 0, 5.0, 5.0);
    put(2, 0, 0, 1.0, 1.0);
    put(2, 1, 1, 3.0, 1.0);
    put(3, 0, 1, 30.0, 10.0);
    put(3, 1, 1, 10.0, 30.0);
    let lift = vec![0.0; on.len()];
    let path = |name: &str| dir.join(format!("m6a_{name}.parquet"));
    write::block_sites(
        &sites,
        &on,
        &off,
        &lift,
        blocks,
        roles,
        &path("block_sites"),
    )
    .unwrap();
    let mut shift = [0.1, -1.0, 0.0, 0.3, 0.5, 0.0];
    if swap {
        shift.rotate_left(roles);
    }
    write::blocks(&shift, &[0.1; 6], &[1.0; 6], roles, &path("blocks")).unwrap();
    let weights = [0.9, 0.1, 0.0, 0.6, 0.4, 0.0, 0.0, 1.0, 0.0];
    let genes = ["GENE1", "GENE2", "GENE3"];
    write::gene_roles(&genes, &weights, roles, &path("gene_roles")).unwrap();
    let cells: Vec<Cell> = (0..5)
        .map(|i| Cell {
            barcode: format!("B{i}").into(),
            batch: "b".into(),
            fit: true,
        })
        .collect();
    let mut member = [[0.8f32, 0.2], [0.7, 0.3], [0.6, 0.4], [0.1, 0.9]];
    if swap {
        member.iter_mut().for_each(|m| m.reverse());
    }
    let rows: Vec<write::CommunityRow> = cells
        .iter()
        .enumerate()
        .map(|(i, cell)| write::CommunityRow {
            cell,
            membership: member.get(i).map(|m| m.as_slice()),
            level: None,
        })
        .collect();
    write::communities(&rows, 2, "role", &path("cell_roles")).unwrap();
    // The fitted cells' links: B0, B1 and B3.
    let fitted = [0usize, 1, 3];
    let mut links = vec![0.0; fitted.len() * blocks];
    for (r, _) in fitted.iter().enumerate() {
        links[r * blocks] = 1.0;
    }
    let labels = [
        (
            "barcode",
            fitted.iter().map(|&i| &*cells[i].barcode).collect(),
        ),
        ("batch", fitted.iter().map(|&i| &*cells[i].batch).collect()),
    ];
    let enrichment = vec![1.0; links.len()];
    write::block_members(
        &labels,
        &links,
        &enrichment,
        blocks,
        roles,
        &path("block_cells"),
    )
    .unwrap();
    let stats = SiteStats {
        on: vec![20, 10, 4, 40],
        off: vec![20, 10, 2, 40],
        cells: vec![9, 4, 3, 12],
    };
    write::site_stats(&sites, &stats, &path("site_stats")).unwrap();
}

fn loaded() -> Modules {
    let tmp = tempfile::tempdir().unwrap();
    tables(tmp.path());
    Modules::load(tmp.path(), "m6a").unwrap()
}

/// The modules of [`tables_with`] the cell roles swapped.
fn swapped() -> Modules {
    let tmp = tempfile::tempdir().unwrap();
    tables_with(tmp.path(), true);
    Modules::load(tmp.path(), "m6a").unwrap()
}

fn count(on: f64, total: f64) -> Count {
    Count { on, total }
}

fn press(v: &mut ModulesView, code: KeyCode) {
    v.handle_key(KeyEvent::new(code, KeyModifiers::NONE));
}

fn typed(v: &mut ModulesView, text: &str) {
    for ch in text.chars() {
        press(v, KeyCode::Char(ch));
    }
}

/// One frame of `v` on a `w` by `h` terminal, its styles too.
fn drawn(v: &mut ModulesView, w: u16, h: u16) -> Buffer {
    let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
    term.draw(|f| v.render(f)).unwrap();
    term.backend().buffer().clone()
}

fn screen(v: &mut ModulesView, w: u16, h: u16) -> String {
    let buf = drawn(v, w, h);
    (0..h)
        .map(|y| (0..w).map(|x| buf[(x, y)].symbol()).collect::<String>())
        .collect::<Vec<_>>()
        .join("\n")
}

fn mouse(v: &mut ModulesView, kind: MouseEventKind, at: (u16, u16), modifiers: KeyModifiers) {
    v.mouse(MouseEvent {
        kind,
        column: at.0,
        row: at.1,
        modifiers,
    });
}

/// The screen cell of the first line holding `text`, at its first character.
fn cell_of(s: &str, text: &str) -> (u16, u16) {
    let (y, line) = s
        .lines()
        .enumerate()
        .find(|(_, l)| l.contains(text))
        .unwrap();
    let x = line[..line.find(text).unwrap()].chars().count();
    (x as u16, y as u16)
}

/// The line of the modules grid's column heads: its corner and the names.
fn heads(s: &str) -> String {
    let below = s
        .lines()
        .position(|l| l.contains("gene modules ·"))
        .unwrap();
    s.lines().nth(below + 1).unwrap().to_string()
}
