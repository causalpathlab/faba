//! `--go-labels N` on the modules grid: each gene module's representative
//! terms on its row, a line each with its q after it, as wide as the page
//! leaves, and a row with none blank unless `--go-none`.

use super::super::figure::{draw_with, go_table, Look};
use super::go_column::{with_go, with_sets};
use super::margins::{grid, grid_of, paper, texts, Text, F};
use crate::figure::style::{Palette, Sheet, StyleArgs};

/// Where `t` ends on the right, measured as GO text is laid out.
fn end(t: &Text) -> f64 {
    let wide = super::super::figure::go_width(&t.s, F);
    match t.anchor.as_str() {
        "start" => t.x + wide,
        "end" => t.x,
        _ => t.x + wide / 2.0,
    }
}

/// The grid's row of gene module 1, the one with terms.
fn first(m: &super::Modules) -> usize {
    m.diagonal.iter().position(|&g| g == 0).unwrap()
}

/// `w` by `h` points, 7-point text, the plot alone, in grays.
fn page(w: f64, h: f64) -> Sheet {
    Sheet {
        page: (Some(w), Some(h)),
        look: StyleArgs {
            font_size: Some(F),
            panel_only: true,
            palette: Palette::Gray,
        },
    }
}

/// The texts right of the grid on the row centred at `y`, top down.
fn row_texts(svg: &str, right: f64, y: f64, cell: f64) -> Vec<Text> {
    let mut row: Vec<Text> = (texts(svg).into_iter())
        .filter(|t| t.x > right && (t.y - y).abs() < cell / 2.0)
        .collect();
    row.sort_by(|a, b| a.y.total_cmp(&b.y));
    row
}

/// Each row carries its module's representative terms, up to N, each its
/// name (shortened only with an ellipsis) and then its q as the browser
/// writes it; a row with none is blank, and says so with `--go-none`; GO
/// is read only with the flag.
#[test]
fn each_row_lists_its_representative_terms_with_their_q_after_them() {
    let tmp = tempfile::tempdir().unwrap();
    let (m, mut ask) = with_go(tmp.path());
    let table = go_table(&m, Some(3), &ask).unwrap().unwrap();
    for (n, none) in [(1, false), (2, false), (3, true)] {
        let g = grid(|g| {
            g.go_labels = Some(n);
            g.go_none = none;
        });
        let sheet = paper(Palette::Gray);
        let svg = draw_with(&m, Look::Hinton, &g, Some(&table), "t", sheet).unwrap();
        let at = grid_of(&svg, &m);
        for (&g, &y) in m.diagonal.iter().zip(&at.rows) {
            let row = row_texts(&svg, at.right(), y, at.cell);
            let terms = table.of(g).unwrap().representative(n);
            if terms.is_empty() {
                let said: Vec<&str> = row.iter().map(|t| t.s.as_str()).collect();
                match none {
                    true => assert_eq!(said, ["none at q <= 0.05"]),
                    false => assert!(said.is_empty(), "{said:?}"),
                }
                continue;
            }
            assert_eq!(row.len(), terms.len(), "M{}: {row:?}", g + 1);
            for (t, h) in row.iter().zip(&terms) {
                let q = format!(" (q {})", super::super::go::q_said(h.q));
                let name = t.s.strip_suffix(&q).unwrap_or_else(|| panic!("{t:?}"));
                let name = name.trim_end_matches('…');
                assert!(h.named().starts_with(name) && name.len() > 3, "{t:?}");
            }
        }
    }
    // The fixture's sets overlap: the second term shown is not the
    // second by p, which is mostly within the first.
    let m1 = table.of(0).unwrap();
    assert_ne!(m1.representative(2)[1].id, m1.hits[1].id);
    // Not read without the flag: a source that is not there is no error.
    ask.gmt = Some(tmp.path().join("none.gmt").to_string_lossy().into());
    assert!(go_table(&m, None, &ask).unwrap().is_none());
    assert!(go_table(&m, Some(2), &ask).is_err());
}

/// A row's N terms each on a line of its own within the row's band, the
/// lines no closer than the font, whenever the row holds N lines and half
/// a line between rows; fewer terms where it does not; every line's
/// capitals a third of the font inside the band, so a row's lines never
/// run into the next row's: over pages of many heights.
#[test]
fn n_terms_show_when_the_row_holds_n_lines() {
    let tmp = tempfile::tempdir().unwrap();
    let (m, ask) = with_go(tmp.path());
    for n in [2, 3] {
        let table = go_table(&m, Some(n), &ask).unwrap().unwrap();
        let has = table.of(0).unwrap().representative(n).len();
        assert_eq!(has, n, "the fixture's M1 has {n} terms");
        let g = grid(|g| g.go_labels = Some(n));
        let (mut full, mut fewer) = (0, 0);
        for h in (40..200).step_by(3) {
            let sheet = page(540.0, h as f64);
            let Ok(svg) = draw_with(&m, Look::Hinton, &g, Some(&table), "t", sheet) else {
                continue;
            };
            let at = grid_of(&svg, &m);
            let (y, cell) = (at.rows[first(&m)], at.cell);
            let row = row_texts(&svg, at.right(), y, cell);
            let holds = ((cell - 0.5 * F) / F).floor().max(1.0) as usize;
            assert_eq!(row.len(), n.min(holds), "cell {cell}: {row:?}");
            match row.len() == n {
                true => full += 1,
                false => fewer += 1,
            }
            for pair in row.windows(2) {
                assert!(pair[1].y - pair[0].y >= F - 1e-9, "cell {cell}: {row:?}");
            }
            for t in &row {
                let inside = F / 3.0 - 1e-9;
                assert!(
                    t.y - 0.75 * F >= y - cell / 2.0 + inside,
                    "cell {cell}: {t:?}"
                );
                assert!(t.y <= y + cell / 2.0 - inside, "cell {cell}: {t:?}");
            }
        }
        assert!(full >= 5 && fewer >= 1, "N {n}: {full} full, {fewer} fewer");
    }
}

/// Long names take the width the page leaves after the grid: the
/// rightmost text ends near the page's right edge, none runs off it, and
/// a name is shortened at a word's end, with an ellipsis; short names
/// take only what they need.
#[test]
fn long_names_take_the_width_the_page_leaves_and_no_more() {
    let tmp = tempfile::tempdir().unwrap();
    let long = |k: usize| {
        format!("SET{k}_A_VERY_LONG_DESCRIPTIVE_NAME_OF_A_GENE_SET_THAT_GOES_ON_AND_ON_FOR_A_WHILE")
    };
    let (m, ask) = with_sets(tmp.path(), vec![(long(1), 0..10), (long(2), 10..20)]);
    let table = go_table(&m, Some(2), &ask).unwrap().unwrap();
    let g = grid(|g| g.go_labels = Some(2));
    let (w, h) = (540.0, 230.4);
    let svg = draw_with(&m, Look::Hinton, &g, Some(&table), "t", page(w, h)).unwrap();
    for t in texts(&svg) {
        assert!(end(&t) <= w + 1e-9, "{t:?} ends at {}", end(&t));
    }
    let rightmost = texts(&svg).iter().map(end).fold(0.0, f64::max);
    assert!(rightmost >= w - 3.0 * F, "{rightmost}");
    let at = grid_of(&svg, &m);
    let row = row_texts(&svg, at.right(), at.rows[first(&m)], at.cell);
    assert_eq!(row.len(), 2, "{row:?}");
    for t in &row {
        let name = t.s.split(" (q ").next().unwrap();
        let kept = name.strip_suffix('…').expect("shortened");
        assert!(
            long(1).starts_with(kept) || long(2).starts_with(kept),
            "{t:?}"
        );
        let next = long(1).chars().nth(kept.chars().count());
        assert!(matches!(next, Some('_' | ' ' | '-')), "{t:?}");
    }
    let (m, ask) = with_sets(tmp.path(), vec![("SET1_A".into(), 0..10)]);
    let table = go_table(&m, Some(2), &ask).unwrap().unwrap();
    let svg = draw_with(&m, Look::Hinton, &g, Some(&table), "t", page(w, h)).unwrap();
    let rightmost = texts(&svg).iter().map(end).fold(0.0, f64::max);
    assert!(rightmost < w - 10.0 * F, "{rightmost}");
    // Either way the grid is as large as with no labels: the height
    // bounds it, and the labels take only the width left.
    let bare = draw_with(&m, Look::Hinton, &grid(|_| {}), None, "t", page(w, h)).unwrap();
    let cell = grid_of(&bare, &m).cell;
    assert!((grid_of(&svg, &m).cell - cell).abs() < 1e-9);
    assert!((at.cell - cell).abs() < 1e-9);
}
