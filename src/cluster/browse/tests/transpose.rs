//! The transposed grid (cell modules as rows, gene modules as columns),
//! the screen's genes beside it on its rows, and the genes' scale and
//! look: where the model goes, the "all" track and the tracks' labels.

use super::super::figure::screen::{self, GenesLook, ModelAt, Region, State, YScale};
use super::super::figure::{draw_with, GridArgs, Look};
use super::*;
use crate::figure::style::Sheet;
use rustc_hash::FxHashMap;

fn screen_of(m: &Modules, state: &State, grid: &GridArgs, genes: GenesLook) -> String {
    let placed = state.place(m).unwrap();
    let models = FxHashMap::default();
    let ask = screen::Ask {
        go: None,
        models: &models,
        grid,
        region: Region::All,
        genes,
        title: "a fit",
        sheet: Sheet::OWN,
    };
    screen::draw(m, &placed, &ask).unwrap()
}

/// The `y` attribute of the first `<text>` saying `said` exactly.
fn text_y(svg: &str, said: &str) -> Option<f64> {
    let end = svg.find(&format!(">{said}</text>"))?;
    let tag = &svg[svg[..end].rfind("<text")?..end];
    let at = tag.find(" y=\"")? + 4;
    tag[at..].split('"').next()?.parse().ok()
}

/// Transposed, the gene modules' names are turned over the columns and
/// the cell modules' are upright left of the rows, under every look.
#[test]
fn a_transposed_grid_turns_the_gene_modules_names() {
    let m = loaded();
    let grid = GridArgs {
        transpose: true,
        size_bars: vec![
            super::super::figure::Side::Left,
            super::super::figure::Side::Top,
        ],
        ..GridArgs::default()
    };
    for look in [
        Look::Shift,
        Look::Percent,
        Look::Z,
        Look::Hinton,
        Look::Bars,
    ] {
        let svg = draw_with(&m, look, &grid, None, "a fit", Sheet::OWN).unwrap();
        assert!(
            svg.contains("rotate(-90)") && svg.contains(">M1 ("),
            "{look:?}"
        );
        // C1 upright, left of the rows: a plain text, not turned.
        let c1 = svg.find(">C1</text>").expect("C1");
        let tag = &svg[svg[..c1].rfind("<text").unwrap()..c1];
        assert!(!tag.contains("rotate"), "{look:?}: {tag}");
        assert!(text_y(&svg, "C2").unwrap() > text_y(&svg, "C1").unwrap());
    }
    // Upright, nothing is turned.
    let svg = draw_with(
        &m,
        Look::Hinton,
        &GridArgs::default(),
        None,
        "a fit",
        Sheet::OWN,
    )
    .unwrap();
    assert!(!svg.contains("rotate(-90)"));
}

/// On a transposed screen the pinned genes' tracks are on the grid's
/// rows: each cell module's track label level with its row's name.
#[test]
fn a_transposed_screens_genes_share_the_grids_rows() {
    let m = loaded();
    let state = State {
        pins: vec!["GENE1".into(), "GENE3".into()],
        side: screen::Beside::None,
        ..State::default()
    };
    let grid = GridArgs {
        transpose: true,
        ..GridArgs::default()
    };
    let svg = screen_of(&m, &state, &grid, GenesLook::default());
    // Each of C1 and C2 written twice (the grid's row, the genes' track)
    // at one height on the page: the grid's in a figure placed over it.
    let attr = |tag: &str, name: &str| -> Option<f64> {
        let at = tag.find(&format!(" {name}=\""))? + name.len() + 3;
        tag[at..].split('"').next()?.parse().ok()
    };
    for c in ["C1", "C2"] {
        let ys: Vec<f64> = (svg.match_indices(&format!(">{c}</text>")))
            .filter_map(|(at, _)| {
                let tag = &svg[svg[..at].rfind("<text").unwrap()..at];
                let placed = (svg[..at].rfind("<svg x=")).and_then(|s| attr(&svg[s..at], "y"));
                Some(attr(tag, "y")? + placed.unwrap_or(0.0))
            })
            .collect();
        assert_eq!(ys.len(), 2, "{c}: {ys:?}");
        assert!((ys[0] - ys[1]).abs() < 0.6, "{c}: {ys:?}");
    }
}

/// The genes' look: "all" and the labels left out, the name over the
/// tracks when the model goes under them.
#[test]
fn the_genes_look_leaves_out_what_it_is_asked_to() {
    let m = loaded();
    let state = State {
        pins: vec!["GENE1".into()],
        side: screen::Beside::None,
        ..State::default()
    };
    let grid = GridArgs::default();
    let plain = screen_of(&m, &state, &grid, GenesLook::default());
    assert!(plain.contains(">all</text>"));
    let bare = GenesLook {
        no_all: true,
        no_track_labels: true,
        ..GenesLook::default()
    };
    let svg = screen_of(&m, &state, &grid, bare);
    assert!(!svg.contains(">all</text>"));
    let under = |look| text_y(&screen_of(&m, &state, &grid, look), "GENE1 M1").unwrap();
    let flipped = GenesLook {
        model_at: ModelAt::Bottom,
        ..GenesLook::default()
    };
    assert!(under(flipped) < under(GenesLook::default()));
}

/// The scales: each side's share of its peak, the peak filling it; the
/// smaller sides larger on sqrt and log than on linear.
#[test]
fn the_genes_scales_keep_the_peak_and_lift_the_rest() {
    for s in [YScale::Linear, YScale::Sqrt, YScale::Log10] {
        assert!((s.share(40.0, 40.0) - 1.0).abs() < 1e-12);
        assert_eq!(s.share(0.0, 40.0), 0.0);
        assert_eq!(s.share(3.0, 0.0), 0.0);
    }
    let lin = YScale::Linear.share(4.0, 400.0);
    assert!(YScale::Sqrt.share(4.0, 400.0) > lin && YScale::Log10.share(4.0, 400.0) > lin);
    assert_eq!(YScale::Linear.next().next().next(), YScale::Linear);
    let json = serde_json::to_string(&State {
        genes_yscale: YScale::Log10,
        ..State::default()
    })
    .unwrap();
    assert!(json.contains(r#""genes_yscale":"log10""#), "{json}");
    // As the pileup names it, and as `log1p` too.
    for said in ["log1p", "log10"] {
        let read: State = serde_json::from_str(&format!(r#"{{"genes_yscale":"{said}"}}"#)).unwrap();
        assert_eq!(read.genes_yscale, YScale::Log10, "{said}");
    }
}

/// Transposed, a gene module's GO terms are turned under its column (or
/// over the names with `--go-side top`), never right of the grid.
#[test]
fn a_transposed_grids_go_terms_are_turned_by_its_columns() {
    use super::go_column::with_go;
    let tmp = tempfile::tempdir().unwrap();
    let (m, ask) = with_go(tmp.path());
    let table = super::super::figure::go_table(&m, Some(2), &ask)
        .unwrap()
        .unwrap();
    let ys = |side| {
        let grid = GridArgs {
            transpose: true,
            go_labels: Some(2),
            go_side: side,
            ..GridArgs::default()
        };
        let svg = draw_with(&m, Look::Hinton, &grid, Some(&table), "a fit", Sheet::OWN).unwrap();
        if let Ok(dir) = std::env::var("FABA_TEST_SVG") {
            std::fs::write(format!("{dir}/go_{side:?}.svg"), &svg).unwrap();
        }
        // The turned texts with a q: the terms, and where they are.
        let turned: Vec<f64> = (svg.lines())
            .filter(|l| l.contains("rotate(-90)") && l.contains("(q "))
            .map(|l| {
                let at = l.find("translate(").unwrap() + 10;
                l[at..]
                    .split(')')
                    .next()
                    .unwrap()
                    .split(' ')
                    .nth(1)
                    .unwrap()
                    .parse()
                    .unwrap()
            })
            .collect();
        assert!(!turned.is_empty(), "{side:?}: no turned term");
        (turned, text_y(&svg, "C1").unwrap())
    };
    let (under, c1) = ys(None);
    assert!(under.iter().all(|&y| y > c1), "{under:?} under {c1}");
    let (over, c1) = ys(Some(super::super::figure::GoSide::Top));
    assert!(over.iter().all(|&y| y < c1), "{over:?} over {c1}");
}

/// `--grid-cap`: the colours topped at the quantile asked, the key's ends
/// saying so; uncapped, they are not marked.
#[test]
fn a_capped_grid_says_so_at_its_keys_ends() {
    let m = loaded();
    let capped = GridArgs {
        cap: Some(0.5),
        ..GridArgs::default()
    };
    let plain = draw_with(
        &m,
        Look::Shift,
        &GridArgs::default(),
        None,
        "a fit",
        Sheet::OWN,
    )
    .unwrap();
    let svg = draw_with(&m, Look::Shift, &capped, None, "a fit", Sheet::OWN).unwrap();
    assert!(!plain.contains('≥') && !plain.contains('≤'));
    assert!(svg.contains('≥') && svg.contains('≤'), "{svg}");
    assert_ne!(plain, svg);
}
