//! The screen as a figure: its state named and found in a fit, the cell
//! modules put in its order, each region drawn from it, the browser's `s`
//! saving the figure and its state, and `--screen` writing one screen or
//! many, each file named for what it shows.

use super::super::figure::screen::{self, Beside, Preset, Region, State};
use super::super::figure::GridArgs;
use super::super::view::Pin;
use super::*;
use crate::figure::style::Sheet;
use rustc_hash::FxHashMap;

/// The screen of `m` in `state`, the region asked for, at its own size.
fn drawn_screen(m: &Modules, state: &State, region: Region) -> Result<String, String> {
    let placed = state.place(m).map_err(|e| e.to_string())?;
    let models = FxHashMap::default();
    let ask = screen::Ask {
        go: None,
        models: &models,
        grid: &GridArgs::default(),
        region,
        genes: Default::default(),
        title: "a fit",
        sheet: Sheet::OWN,
    };
    screen::draw(m, &placed, &ask)
}

/// A state of the fixture: M1 in focus, its sites beside the grid, GENE1
/// pinned.
fn focused() -> State {
    State {
        module: Some("M1".into()),
        pins: vec!["GENE1".into()],
        ..State::default()
    }
}

/// A state says everything by the names the browser shows, and reads back
/// as it was written; a field it does not have is refused, one left out
/// is the browser's own.
#[test]
fn a_state_is_plain_json_by_the_browsers_names() {
    let state = State {
        preset: Preset::Blocks,
        cell_order: vec!["C2".into()],
        module: Some("M2".into()),
        side: Beside::Go,
        go_terms: Some(2),
        pins: vec!["GENE3".into(), "GENE1".into()],
        genes_yscale: screen::YScale::Sqrt,
    };
    let json = serde_json::to_string(&state).unwrap();
    assert!(json.contains(r#""preset":"blocks""#), "{json}");
    assert!(json.contains(r#""side":"go""#), "{json}");
    assert_eq!(serde_json::from_str::<State>(&json).unwrap(), state);
    assert_eq!(
        serde_json::from_str::<State>("{}").unwrap(),
        State::default()
    );
    assert!(serde_json::from_str::<State>(r#"{"look":"bars"}"#).is_err());
}

/// The names are found in the fit: a cell module as C2 or 2, the gene
/// module in focus, each pinned gene in its main module and in order; a
/// name the fit lacks is refused, saying so.
#[test]
fn a_state_is_placed_in_a_fit_by_its_names() {
    let m = loaded();
    let state = State {
        cell_order: vec!["2".into()],
        module: Some("M2".into()),
        pins: vec!["GENE3".into(), "GENE1".into(), "GENE3".into()],
        ..State::default()
    };
    let placed = state.place(&m).unwrap();
    assert_eq!(placed.columns[0], 1);
    assert_eq!(placed.module, Some(1));
    // GENE3 in M2, GENE1 in M1, each once, in the order named.
    assert_eq!(
        placed.pins,
        vec![Pin { role: 1, gene: 0 }, Pin { role: 0, gene: 0 }]
    );
    // Left out: the grid's first row in focus.
    let own = State::default().place(&m).unwrap();
    assert_eq!(own.module, m.diagonal.first().copied());
    for (bad, said) in [
        (
            State {
                cell_order: vec!["C9".into()],
                ..State::default()
            },
            "C9",
        ),
        (
            State {
                module: Some("M9".into()),
                ..State::default()
            },
            "M9",
        ),
        (
            State {
                pins: vec!["GENE9".into()],
                ..State::default()
            },
            "GENE9",
        ),
    ] {
        let e = bad.place(&m).unwrap_err().to_string();
        assert!(e.contains(said), "{e}");
    }
}

/// Ordered as a state names them, the cell modules are shown so and every
/// block's reads move with their column; a role named twice is refused.
#[test]
fn the_cell_modules_are_put_in_the_states_order() {
    let mut m = loaded();
    let was: Vec<(usize, Count)> = m.columns.iter().copied().zip(m.block(0).to_vec()).collect();
    let last = *m.columns.last().unwrap();
    let state = State {
        cell_order: vec![format!("C{}", last + 1)],
        ..State::default()
    };
    state.order(&mut m).unwrap();
    assert_eq!(m.columns[0], last);
    for (k, &c) in m.columns.iter().enumerate() {
        let before = was.iter().find(|(r, _)| *r == c).unwrap().1;
        assert_eq!(m.block(0)[k], before, "C{}", c + 1);
    }
    assert!(m.reorder_columns(&[0, 0]).is_err());
    // Drawn only once the fit is in the state's order.
    let fresh = loaded();
    let swapped = State {
        cell_order: vec![format!("C{}", fresh.columns[1] + 1)],
        ..State::default()
    };
    assert!(drawn_screen(&fresh, &swapped, Region::All).is_err());
}

/// The whole screen draws its three regions, each titled as the
/// browser's panes are: the grid with each row's genes and shift range,
/// the module in focus's sites (its pinned gene's marked), and the pinned
/// gene's tracks, all cells first, then each cell module with its size.
#[test]
fn the_screen_draws_the_grid_the_sites_and_the_pinned_genes() {
    let m = loaded();
    let svg = drawn_screen(&m, &focused(), Region::All).unwrap();
    for said in [
        "gene modules · diagonal",
        "sites · M1",
        "genes · 1 pinned",
        "Δ0.20",
        "• GENE1 100",
        "GENE1 M1",
        ">all<",
    ] {
        assert!(svg.contains(said), "{said}");
    }
    // The blocks preset: the rows by shift range, as bars.
    let blocks = State {
        preset: Preset::Blocks,
        ..focused()
    };
    let svg = drawn_screen(&m, &blocks, Region::All).unwrap();
    assert!(svg.contains("by shift range · bars"));
    // Nothing beside the grid: no sites.
    let alone = State {
        side: Beside::None,
        ..focused()
    };
    let svg = drawn_screen(&m, &alone, Region::All).unwrap();
    assert!(!svg.contains("sites ·"));
}

/// A region alone is that region: the grid has no sites or genes, the
/// genes no grid; one with nothing to draw is refused.
#[test]
fn a_region_is_drawn_alone_or_refused_when_empty() {
    let m = loaded();
    let grid = drawn_screen(&m, &focused(), Region::Grid).unwrap();
    assert!(grid.contains("gene modules") && !grid.contains("GENE1 M1"));
    let genes = drawn_screen(&m, &focused(), Region::Genes).unwrap();
    assert!(genes.contains("GENE1 M1") && !genes.contains("gene modules"));
    let side = drawn_screen(&m, &focused(), Region::Side).unwrap();
    assert!(side.contains("GENE1 300") && !side.contains("Δ"));
    let unpinned = State {
        pins: Vec::new(),
        ..focused()
    };
    let e = drawn_screen(&m, &unpinned, Region::Genes).unwrap_err();
    assert!(e.contains("pinned"), "{e}");
    let nothing = State {
        side: Beside::None,
        ..focused()
    };
    assert!(drawn_screen(&m, &nothing, Region::Side).is_err());
    // GO beside the grid, with no terms read: refused.
    let go = State {
        side: Beside::Go,
        ..focused()
    };
    assert!(drawn_screen(&m, &go, Region::Side).is_err());
}

/// The pinned genes go in the order pinned, and with the plot alone the
/// regions have no titles.
#[test]
fn the_pins_keep_their_order_and_the_plot_alone_has_no_titles() {
    let m = loaded();
    let state = State {
        pins: vec!["GENE3".into(), "GENE1".into()],
        ..focused()
    };
    let svg = drawn_screen(&m, &state, Region::Genes).unwrap();
    let at = |s: &str| svg.find(s).unwrap();
    assert!(at("GENE3 M2") < at("GENE1 M1"));
    let placed = state.place(&m).unwrap();
    let models = FxHashMap::default();
    let mut sheet = Sheet::OWN;
    sheet.look.panel_only = true;
    let ask = screen::Ask {
        go: None,
        models: &models,
        grid: &GridArgs::default(),
        region: Region::All,
        genes: Default::default(),
        title: "a fit",
        sheet,
    };
    let svg = screen::draw(&m, &placed, &ask).unwrap();
    assert!(!svg.contains("genes · ") && !svg.contains("a fit"));
}

/// The pinned genes' tracks are the browser's mirror: each column's
/// methylated reads in the accent rising from the track's midline, its
/// unmethylated in a light gray hanging from it, each bar a terminal
/// column wide at the page's scale and any side with reads at least an
/// eighth of a line deep, as the browser's least block; the sites pane's
/// bars hang in the same light gray.
#[test]
fn the_pinned_genes_are_the_browsers_mirror_bars() {
    use super::margins::rects;
    let m = loaded();
    let f = 8.0;
    let placed = focused().place(&m).unwrap();
    let models = FxHashMap::default();
    let mut sheet = Sheet::OWN;
    sheet.look.font_size = Some(f);
    let svg = |region| {
        let ask = screen::Ask {
            go: None,
            models: &models,
            grid: &GridArgs::default(),
            region,
            genes: Default::default(),
            title: "a fit",
            sheet,
        };
        screen::draw(&m, &placed, &ask).unwrap()
    };
    let light = |fill: &str| {
        let v = |k: usize| u8::from_str_radix(&fill[k..k + 2], 16).unwrap();
        fill.len() == 7 && v(1) == v(3) && v(3) == v(5) && v(1) >= 0xa0 && v(1) < 0xf0
    };
    let all = rects(&svg(Region::Genes));
    let up: Vec<_> = all
        .iter()
        .filter(|b| b.fill == crate::figure::ACCENT)
        .collect();
    let off: Vec<_> = all.iter().filter(|b| light(&b.fill)).collect();
    assert!(!up.is_empty() && !off.is_empty(), "{all:?}");
    for b in up.iter().chain(&off) {
        assert!(b.w >= 0.6 * f - 0.01, "{b:?} narrower than a column");
        assert!(
            b.h >= 1.2 * f / 8.0 - 0.01,
            "{b:?} under an eighth of a line"
        );
    }
    // Each methylated bar stands on the midline its unmethylated hangs from.
    for b in &up {
        let under = off
            .iter()
            .find(|o| (o.x - b.x).abs() < 0.01 && (o.y - (b.y + b.h)).abs() < 0.01);
        assert!(under.is_some(), "{b:?} has nothing hanging under it");
    }
    let sites = rects(&svg(Region::Side));
    assert!(sites.iter().any(|b| light(&b.fill)), "{sites:?}");
    assert!(!sites.iter().any(|b| b.fill == crate::figure::MUTED));
}

/// GO beside the grid: a gene module with no term is left blank, as the
/// grid figure's `--go-labels` leaves it, and says so with `--go-none`.
#[test]
fn a_go_row_with_no_term_is_blank_unless_go_none() {
    let tmp = tempfile::tempdir().unwrap();
    let (m, ask) = super::go_column::with_go(tmp.path());
    let table = super::super::figure::go_table(&m, Some(3), &ask)
        .unwrap()
        .unwrap();
    assert!(m
        .diagonal
        .iter()
        .any(|&g| table.of(g).is_none_or(|t| t.representative(3).is_empty())));
    let state = State {
        side: Beside::Go,
        ..focused()
    };
    let placed = state.place(&m).unwrap();
    let models = FxHashMap::default();
    let svg = |go_none| {
        let grid = GridArgs {
            go_none,
            ..GridArgs::default()
        };
        let ask = screen::Ask {
            go: Some(&table),
            models: &models,
            grid: &grid,
            region: Region::All,
            genes: Default::default(),
            title: "a fit",
            sheet: Sheet::OWN,
        };
        screen::draw(&m, &placed, &ask).unwrap()
    };
    let (plain, said) = (svg(false), svg(true));
    assert!(plain.contains("TERM1_SET") && !plain.contains("none at q"));
    assert!(said.contains("none at q &lt;= 0.05") || said.contains("none at q <= 0.05"));
}

/// The browser's state is what it shows: the look's preset, the cell
/// modules in the panes' order, the module in focus, the pins by key.
#[test]
fn the_browsers_state_is_what_it_shows() {
    let m = loaded();
    let mut v = ModulesView::new("t", &m);
    drawn(&mut v, 160, 40);
    let state = v.screen_state();
    assert_eq!(state.preset, Preset::Diagonal);
    assert_eq!(state.module, Some(format!("M{}", v.module + 1)));
    let cols: Vec<String> = m.columns.iter().map(|c| format!("C{}", c + 1)).collect();
    assert_eq!(state.cell_order, cols);
    v.select_module(0);
    press(&mut v, KeyCode::Tab);
    press(&mut v, KeyCode::Char(' '));
    bars(&mut v);
    let state = v.screen_state();
    assert_eq!(state.preset, Preset::Blocks);
    assert_eq!(state.pins, vec!["GENE1".to_string()]);
    assert_eq!(state.module.as_deref(), Some("M1"));
    // Drawn as `--screen` draws that state.
    let svg = v.screen_svg().unwrap();
    assert!(svg.contains("GENE1 M1") && svg.contains("by shift range"));
}

/// The screen the browser saves is titled as the grid figure is, by the
/// fit's folder name, its modality and how many modules, not by the
/// browser's header (the whole path and the model).
#[test]
fn the_saved_screen_is_titled_by_the_fits_name() {
    let m = loaded();
    let dir = Path::new("/a/b/FIT1");
    let mut v = super::super::opened_view(dir, &m);
    drawn(&mut v, 160, 40);
    let svg = v.screen_svg().unwrap();
    let title = format!(
        "FIT1 · {} · {} gene modules × {} cell modules",
        m.modality,
        m.order.len(),
        m.cell_roles
    );
    assert!(svg.contains(&title), "{title}");
    assert!(!svg.contains("/a/b"));
    assert_eq!(
        super::super::figure_title(dir, &m, Some("C2")),
        title.replacen("FIT1", "FIT1 · C2", 1)
    );
}

/// `s` in the browser saves the screen as a PDF under the name typed,
/// and its state beside it, which reads back as the browser's.
#[test]
fn s_saves_the_screen_and_its_state() {
    let m = loaded();
    let tmp = tempfile::tempdir().unwrap();
    let mut v = ModulesView::new("t", &m);
    drawn(&mut v, 160, 40);
    press(&mut v, KeyCode::Tab);
    press(&mut v, KeyCode::Char(' '));
    press(&mut v, KeyCode::Char('s'));
    // The prompt takes every key: clear the default name, type ours.
    for _ in 0..80 {
        press(&mut v, KeyCode::Backspace);
    }
    let prefix = tmp.path().join("shot");
    typed(&mut v, &prefix.to_string_lossy());
    press(&mut v, KeyCode::Enter);
    assert!(tmp.path().join("shot.pdf").is_file());
    let text = std::fs::read_to_string(tmp.path().join("shot.state.json")).unwrap();
    let state: State = serde_json::from_str(&text).unwrap();
    assert_eq!(state, v.screen_state());
    assert!(v
        .help_keys()
        .contains(&("s", "save the screen as a figure (PDF)")));
}

/// `faba cluster-view` with `args`, run.
fn cluster_view(args: &[&str]) -> anyhow::Result<()> {
    let line = std::iter::once("faba")
        .chain(std::iter::once("cluster-view"))
        .chain(args.iter().copied());
    let matches = crate::faba_command().try_get_matches_from(line)?;
    match crate::parse(&matches)?.commands {
        crate::Commands::ClusterView(ref args) => crate::run_modules(args),
        _ => unreachable!(),
    }
}

fn names(dir: &Path) -> Vec<String> {
    let mut names: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    names.sort();
    names
}

/// `--screen` writes the screen as a PDF and a PNG (or the kinds
/// `--format` names), never over a file already there; the options name
/// the state, and a name the fit lacks is refused before anything is
/// written.
#[test]
fn screen_writes_one_screen_and_never_replaces_a_file() {
    let tmp = tempfile::tempdir().unwrap();
    let fit = tmp.path().join("fit");
    std::fs::create_dir(&fit).unwrap();
    tables(&fit);
    let out = tmp.path().join("out");
    std::fs::create_dir(&out).unwrap();
    let prefix = out.join("shot");
    let (fit, prefix) = (fit.to_string_lossy(), prefix.to_string_lossy());
    let given = [
        &fit[..],
        "--screen",
        &prefix,
        "--pin",
        "GENE1,GENE3",
        "--module",
        "M1",
        "--cell-order",
        "C2",
        "--fig-width",
        "7",
    ];
    cluster_view(&given).unwrap();
    assert_eq!(names(&out), ["shot.pdf", "shot.png"]);
    let e = cluster_view(&given).unwrap_err().to_string();
    assert!(e.contains("exists"), "{e}");
    let svg = out.join("other");
    let other = svg.to_string_lossy();
    cluster_view(&[
        &fit, "--screen", &other, "--format", "svg", "--region", "grid",
    ])
    .unwrap();
    assert!(out.join("other.svg").is_file());
    let e = cluster_view(&[&fit, "--screen", &other, "--pin", "GENE9"]).unwrap_err();
    assert!(e.to_string().contains("GENE9"), "{e}");
    assert!(!out.join("other.pdf").exists());
}

/// A saved state draws the same screen again, each option given over it.
#[test]
fn screen_draws_a_saved_state_again() {
    let tmp = tempfile::tempdir().unwrap();
    tables(tmp.path());
    let state = tmp.path().join("view.state.json");
    std::fs::write(&state, serde_json::to_string(&focused()).unwrap()).unwrap();
    let out = tmp.path().join("again");
    let (dir, state, out) = (
        tmp.path().to_string_lossy(),
        state.to_string_lossy(),
        out.to_string_lossy(),
    );
    cluster_view(&[&dir, "--screen", &out, "--state", &state, "--format", "svg"]).unwrap();
    let svg = std::fs::read_to_string(format!("{out}.svg")).unwrap();
    assert!(svg.contains("GENE1 M1") && svg.contains("sites · M1"));
    let out2 = format!("{out}2");
    let given = [
        &dir[..],
        "--screen",
        &out2,
        "--state",
        &state,
        "--format",
        "svg",
        "--side",
        "none",
    ];
    cluster_view(&given).unwrap();
    let svg = std::fs::read_to_string(format!("{out2}.svg")).unwrap();
    assert!(!svg.contains("sites ·"));
}

/// Many screens in one call, each named for what it shows: a gene module
/// each, `{module}_screen`; every fit under a folder, `{fit}_..`, by its
/// folder there; `--pin-top` pins each module's first genes.
#[test]
fn screens_of_every_module_and_every_fit_are_named_for_them() {
    let tmp = tempfile::tempdir().unwrap();
    let sweep = tmp.path().join("sweep");
    for fit in ["c02", "c03"] {
        let dir = sweep.join(fit);
        std::fs::create_dir_all(&dir).unwrap();
        tables(&dir);
        let record = serde_json::json!({ "job": "cluster", "results": { "model": {}, "tables": { "modality": "m6a" } } });
        std::fs::write(dir.join("m6a_blockmodel.run.json"), record.to_string()).unwrap();
    }
    let out = tmp.path().join("modules");
    let (fit, out_s) = (
        sweep.join("c02").to_string_lossy().into_owned(),
        out.to_string_lossy().into_owned(),
    );
    cluster_view(&[
        &fit,
        "--screen",
        &out_s,
        "--each-module",
        "--pin-top",
        "1",
        "--format",
        "svg",
    ])
    .unwrap();
    assert_eq!(names(&out), ["M1_screen.svg", "M2_screen.svg"]);
    let m2 = std::fs::read_to_string(out.join("M2_screen.svg")).unwrap();
    assert!(m2.contains("sites · M2") && m2.contains("GENE3 M2") && !m2.contains("GENE1 M1"));
    let all = tmp.path().join("all");
    let (sweep_s, all_s) = (sweep.to_string_lossy(), all.to_string_lossy());
    let given = [
        &sweep_s[..],
        "--screen",
        &all_s,
        "--each-fit",
        "--each-module",
        "--format",
        "pdf",
    ];
    cluster_view(&given).unwrap();
    assert_eq!(
        names(&all),
        [
            "c02_M1_screen.pdf",
            "c02_M2_screen.pdf",
            "c03_M1_screen.pdf",
            "c03_M2_screen.pdf"
        ]
    );
}

/// `--each-group` draws a screen of every fitted group of a within
/// record, `{group}_screen`; of a single fit it is refused.
#[test]
fn screens_of_every_group_of_a_within_record() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    count_tables(dir);
    let fit = |modality: &str| {
        serde_json::json!({ "job": "cluster", "results": { "model": {}, "tables": { "modality": modality } } })
            .to_string()
    };
    std::fs::write(dir.join(crate::cluster::within::COUNT_RECORD), fit("count")).unwrap();
    let mut groups = Vec::new();
    for g in ["C1", "C2"] {
        let sub = dir.join("within/m6a").join(g);
        std::fs::create_dir_all(&sub).unwrap();
        tables(&sub);
        std::fs::write(sub.join("m6a_blockmodel.run.json"), fit("m6a")).unwrap();
        groups.push(serde_json::json!({ "group": g, "cells": 80, "status": "fitted",
            "dir": format!("within/m6a/{g}"), "record": "m6a_blockmodel.run.json", "converged": true }));
    }
    let record = serde_json::json!({ "job": "cluster", "results": { "within": {
        "modality": "m6a", "count_record": dir.join(crate::cluster::within::COUNT_RECORD), "groups": groups,
    } } });
    let path = dir.join(crate::cluster::within::record_name("m6a"));
    std::fs::write(&path, record.to_string()).unwrap();
    let out = dir.join("groups");
    let (path_s, out_s) = (path.to_string_lossy(), out.to_string_lossy());
    cluster_view(&[
        &path_s,
        "--screen",
        &out_s,
        "--each-group",
        "--format",
        "svg",
    ])
    .unwrap();
    assert_eq!(names(&out), ["C1_screen.svg", "C2_screen.svg"]);
    let one = dir.join("within/m6a/C1");
    let e =
        cluster_view(&[&one.to_string_lossy(), "--screen", &out_s, "--each-group"]).unwrap_err();
    assert!(e.to_string().contains("within record"), "{e}");
}

/// The screen's options come only with `--screen`; `--pin` and
/// `--pin-top` are one or the other; the figure's size and look are taken
/// with `--screen` as with `--figure`.
#[test]
fn the_screens_options_come_with_screen() {
    let parse = |args: &[&str]| {
        let line = ["faba", "cluster-view", "fit"].iter().chain(args).copied();
        crate::faba_command().try_get_matches_from(line).map(|_| ())
    };
    assert!(parse(&["--region", "grid"]).is_err());
    assert!(parse(&["--pin", "GENE1"]).is_err());
    assert!(parse(&["--screen", "s", "--pin", "GENE1", "--pin-top", "2"]).is_err());
    assert!(parse(&["--screen", "s", "--figure", "f"]).is_err());
    assert!(parse(&["--screen", "s", "--fig-width", "7", "--font-size", "6"]).is_ok());
    assert!(parse(&["--screen", "s", "--go-labels", "2", "--side", "go"]).is_ok());
    // `--go-none` with the GO beside the screen's grid, as with `--go-labels`.
    assert!(parse(&["--screen", "s", "--side", "go", "--go-none"]).is_ok());
    assert!(parse(&["--figure", "f", "--go-none"]).is_err());
    assert!(parse(&["--fig-width", "7"]).is_err());
}

/// The browser and the screen it saves, side by side, of the fit
/// `FABA_SCREEN_FIT` names (its run record), the GTF `FABA_SCREEN_GTF`
/// names, into `FABA_SCREENSHOTS_DIR`: for reading by eye.
#[test]
#[ignore = "reads a fit named by FABA_SCREEN_FIT"]
fn the_screen_reads_as_the_browser() {
    let (Some(fit), Some(dir)) = (
        std::env::var_os("FABA_SCREEN_FIT"),
        std::env::var_os(crate::tui::screenshot::DIR_VAR),
    ) else {
        return;
    };
    let fit = fit.to_string_lossy().into_owned();
    let (tables, modality, record) = super::super::resolve(&[fit.into()], None, None).unwrap();
    let m = Modules::load(&tables, &modality, record).unwrap();
    let gtf = std::env::var("FABA_SCREEN_GTF").ok();
    let models = crate::site_analysis::pileup::start_gene_models(gtf.as_deref());
    let mut v = super::super::opened_view(&tables, &m).with_genes(models);
    while v.pending.is_some() {
        std::thread::sleep(std::time::Duration::from_millis(200));
        v.tick();
    }
    let module = std::env::var("FABA_SCREEN_MODULE").unwrap_or("M1".into());
    v.select_module(super::super::text::gene_module(&m, &module).unwrap());
    let pins = std::env::var("FABA_SCREEN_PINS").unwrap_or_default();
    for name in pins.split(',').filter(|s| !s.is_empty()) {
        v.pins.push(screen::pin_of(&m, name).unwrap());
    }
    let (w, h) = (200, 56);
    let buf = drawn(&mut v, w, h);
    crate::tui::screenshot::save("screen_tui", &buf);
    let svg = v.screen_svg().unwrap();
    let png = crate::figure::png(&svg, 2.0).unwrap();
    std::fs::write(std::path::Path::new(&dir).join("screen_export.png"), png).unwrap();
}

/// The browser as a live view: its state is the screen's, set again on a
/// fresh browser it draws the same figure `s` saves, and a live export's
/// state (its view and command beside it) draws it again with `--screen
/// --state`; what the view cannot show is refused, saying why.
#[test]
fn a_live_browser_sets_and_exports_its_screen() {
    use crate::tui::live::{sheet_of, snapshot, Live, Look};
    let m = loaded();
    let mut v = ModulesView::new("t", &m);
    drawn(&mut v, 160, 40);
    let first = v.view_state();
    v.set_view_state(&serde_json::json!({
        "preset": "blocks",
        "module": "M1",
        "pins": ["GENE1"],
        "go_terms": 2,
        "unknown": true,
    }))
    .unwrap();
    let state = v.view_state();
    assert_eq!(state["preset"], "blocks");
    assert_eq!(state["module"], "M1");
    assert_eq!(state["pins"], serde_json::json!(["GENE1"]));
    assert_eq!(state["go_terms"], 2);
    assert_eq!(state["cell_order"], first["cell_order"]);
    let sheet = sheet_of(&v, &Look::default()).unwrap();
    let svg = v.saved_figure(sheet).unwrap();
    assert_eq!(svg, v.screen_svg().unwrap(), "the export is `s`'s");
    assert!(svg.contains("GENE1 M1") && svg.contains("by shift range"));

    let mut fresh = ModulesView::new("t", &m);
    drawn(&mut fresh, 160, 40);
    assert_ne!(fresh.screen_svg().unwrap(), svg);
    fresh.set_view_state(&state).unwrap();
    assert_eq!(fresh.view_state(), state);
    assert_eq!(fresh.screen_svg().unwrap(), svg);

    for (bad, says) in [
        (serde_json::json!({ "side": "none" }), "side none"),
        (serde_json::json!({ "side": "go" }), "no GO"),
        (serde_json::json!({ "pins": ["GENE9"] }), "GENE9"),
        (
            serde_json::json!({ "cell_order": ["C2"] }),
            "order it opened with",
        ),
    ] {
        let mut again = ModulesView::new("t", &m);
        drawn(&mut again, 160, 40);
        let e = format!("{:#}", again.set_view_state(&bad).unwrap_err());
        assert!(e.contains(says), "{bad}: {e}");
    }

    // The state a live export keeps draws the same screen with --screen.
    let tmp = tempfile::tempdir().unwrap();
    tables(tmp.path());
    let kept = tmp.path().join("live.view.json");
    std::fs::write(&kept, snapshot(&v).to_string()).unwrap();
    let out = tmp.path().join("again");
    let (dir, kept, out) = (
        tmp.path().to_string_lossy(),
        kept.to_string_lossy(),
        out.to_string_lossy(),
    );
    cluster_view(&[&dir, "--screen", &out, "--state", &kept, "--format", "svg"]).unwrap();
    let again = std::fs::read_to_string(format!("{out}.svg")).unwrap();
    assert!(again.contains("GENE1 M1") && again.contains("by shift range"));
}
