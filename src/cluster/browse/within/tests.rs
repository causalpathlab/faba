//! A within record opened: found from what is named, its groups by name,
//! the fit each opens, the printout, the pages and their links, the list.

use super::*;
use crate::cluster::browse::tests::{count_tables, tables};
use ratatui::crossterm::event::KeyModifiers;
use serde_json::json;

/// A fit's run record of `modality`, as little of one as opens.
fn fit_record(modality: &str) -> String {
    json!({ "job": "cluster", "results": { "model": {}, "tables": { "modality": modality } } })
        .to_string()
}

/// A within record in `dir`, its path: the count fit beside it (its
/// tables and record; the record names one that is gone, as after the
/// folder moved), C1 fitted under `within/m6a/C1`, C2 skipped, bg failed.
fn nested(dir: &Path) -> PathBuf {
    count_tables(dir);
    std::fs::write(dir.join(within::COUNT_RECORD), fit_record("count")).unwrap();
    let sub = dir.join("within/m6a/C1");
    std::fs::create_dir_all(&sub).unwrap();
    tables(&sub);
    std::fs::write(sub.join("m6a_blockmodel.run.json"), fit_record("m6a")).unwrap();
    let record = json!({ "job": "cluster", "results": { "within": {
        "modality": "m6a",
        "count_record": dir.join("gone").join(within::COUNT_RECORD),
        "groups": [
            { "group": "C1", "cells": 80, "status": "fitted", "dir": "within/m6a/C1",
              "record": "m6a_blockmodel.run.json", "rhat": 1.02, "converged": true },
            { "group": "C2", "cells": 9, "status": "skipped", "why": "9 cells to fit" },
            { "group": "bg", "cells": 60, "status": "failed", "why": "no sites past the filters" },
        ],
    } } });
    let path = dir.join(within::record_name("m6a"));
    std::fs::write(&path, record.to_string()).unwrap();
    path
}

fn names(paths: &[&Path]) -> Vec<Box<str>> {
    paths.iter().map(|p| p.to_string_lossy().into()).collect()
}

fn read(path: &Path) -> Nested {
    Nested::read(path).unwrap().unwrap()
}

fn press(list: &mut List, code: KeyCode) {
    list.handle_key(KeyEvent::new(code, KeyModifiers::NONE));
}

/// A within record is found named alone, or in the directory named when
/// that holds no site tables of its own; a count fit asked for there, or
/// another modality's fits, is not it.
#[test]
fn a_within_record_is_found_named_or_in_its_directory() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let record = nested(dir);
    for named in [record.as_path(), dir] {
        let found = Nested::find(&names(&[named]), None, None).unwrap().unwrap();
        assert_eq!(found.modality, "m6a");
        assert_eq!(found.dir(), dir);
        assert_eq!(found.fitted().count(), 1);
        // The count fit beside the record, the one it names being gone.
        assert_eq!(found.count, dir.join(within::COUNT_RECORD));
    }
    let found =
        |inputs: &[&Path], modality| Nested::find(&names(inputs), modality, None).map(Result::ok);
    assert!(found(&[dir], Some("m6a")).unwrap().is_some());
    assert!(found(&[record.as_path()], Some("M6A")).unwrap().is_some());
    assert!(found(&[dir], Some("count")).unwrap().is_none());
    assert!(found(&[dir], Some("atoi")).unwrap().is_none());
    assert!(found(&[], None).unwrap().is_none());
    // The count fit's own record is one fit, not a within record.
    let count = dir.join(within::COUNT_RECORD);
    assert!(found(&[count.as_path()], None).unwrap().is_none());

    let err = found(&[record.as_path()], Some("atoi")).unwrap_err();
    assert!(err.to_string().contains("not --modality atoi"), "{err}");
    let err = found(&[record.as_path(), dir], None).unwrap_err();
    assert!(err.to_string().contains("alone"), "{err}");
    let err = Nested::find(&names(&[record.as_path()]), None, Some("out")).unwrap_err();
    assert!(err.to_string().contains("no -o"), "{err}");

    // Site tables in the directory itself are what it opens.
    tables(dir);
    assert!(found(&[dir], None).unwrap().is_none());
    assert!(found(&[record.as_path()], None).unwrap().is_some());
}

/// A group is named as the count fit's view names its cell module, in any
/// case, or by its number; one that is not there is refused, those there
/// named.
#[test]
fn a_group_is_named_as_its_cell_module() {
    let tmp = tempfile::tempdir().unwrap();
    let nested = read(&nested(tmp.path()));
    for (name, group) in [
        ("C1", "C1"),
        ("c2", "C2"),
        ("1", "C1"),
        ("bg", "bg"),
        ("BG", "bg"),
    ] {
        assert_eq!(nested.group(name).unwrap().group, group, "{name}");
    }
    let err = nested.group("C9").unwrap_err().to_string();
    assert!(
        err.contains("no group C9 in m6a_within.run.json") && err.contains("C1, C2, bg"),
        "{err}"
    );
}

/// The count fit opens from its record, a fitted group from its own under
/// the within record's directory; a group without a fit says what became
/// of it, and a count fit that is gone says where it was.
#[test]
fn a_fit_opens_from_its_own_record() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let record = nested(dir);
    let nested = read(&record);
    let (at, modality, fit) = nested.open(None).unwrap();
    assert_eq!((at.as_path(), modality.as_str()), (dir, "count"));
    assert!(fit.is_some());
    let (at, modality, _) = nested.open(Some(nested.group("C1").unwrap())).unwrap();
    assert_eq!(at, dir.join("within/m6a/C1"));
    assert_eq!(modality, "m6a");
    let err = nested.open(Some(nested.group("C2").unwrap())).unwrap_err();
    assert_eq!(
        err.to_string(),
        "C2 has no fit to open (skipped: 9 cells to fit)"
    );

    std::fs::remove_file(dir.join(within::COUNT_RECORD)).unwrap();
    let err = read(&record).open(None).unwrap_err().to_string();
    assert!(
        err.contains("the count fit's run record") && err.contains("is not there"),
        "{err}"
    );
}

/// The printout: the count fit's modules, a line a group saying what
/// became of its fit, then each fitted group's modules under its name. A
/// fit that cannot be read says so in its place, the rest printed.
#[test]
fn the_printout_is_every_fits_modules() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let record = nested(dir);
    let lines = read(&record).printed();
    assert!(lines[0].starts_with("count fit: "), "{lines:?}");
    let heads: Vec<usize> = (0..lines.len())
        .filter(|&i| lines[i].starts_with("module\tmain genes"))
        .collect();
    assert_eq!(heads.len(), 2, "{lines:?}");
    assert_eq!(heads[0], 1);
    let summary = lines.iter().position(|l| l == "group\tcells\tfit").unwrap();
    assert!(
        lines[summary - 1].starts_with("m6a within its cell modules: 1 of 3 groups fitted ("),
        "{lines:?}"
    );
    assert_eq!(
        lines[summary + 1..summary + 4],
        [
            "C1\t80 cells\tfitted, R-hat 1.02, converged",
            "C2\t9 cells\tskipped: 9 cells to fit",
            "bg\t60 cells\tfailed: no sites past the filters",
        ]
    );
    assert_eq!(lines[heads[1] - 1], "m6a within C1");
    // The site fit's values are % converted, the count fit's UMIs per cell.
    assert!(lines[heads[0]].contains("UMI/cell") && lines[heads[1]].contains(" %"));

    std::fs::remove_file(dir.join("within/m6a/C1/m6a_block_sites.parquet")).unwrap();
    let lines = read(&record).printed();
    let at = lines.iter().position(|l| l == "m6a within C1").unwrap();
    assert!(lines[at + 1].starts_with("cannot be read: "), "{lines:?}");
    assert!(lines[1].starts_with("module\tmain genes"));
}

/// A group's page is beside the count fit's, named after it; the links on
/// a page are the count fit and every group, the page itself and the
/// groups without a fit not linked, each saying what became of its fit.
#[test]
fn the_pages_link_to_each_other() {
    assert_eq!(
        group_page(Path::new("out/my fit.html"), "C1"),
        Path::new("out/my fit.C1.html")
    );
    assert_eq!(
        group_page(Path::new("fit.html"), "bg"),
        Path::new("fit.bg.html")
    );
    assert_eq!(
        href(Path::new("out/my fit#1.C1.html")),
        "my%20fit%231.C1.html"
    );

    let tmp = tempfile::tempdir().unwrap();
    let nested = read(&nested(tmp.path()));
    let html = Path::new("out/fit.html");
    let on_count = nested.links(html, None);
    let labels: Vec<&str> = on_count.iter().map(|l| l.label.as_str()).collect();
    assert_eq!(labels, ["count", "C1", "C2", "bg"]);
    let hrefs =
        |links: &[Link]| -> Vec<Option<String>> { links.iter().map(|l| l.href.clone()).collect() };
    assert_eq!(
        hrefs(&on_count),
        [None, Some("fit.C1.html".into()), None, None]
    );
    assert!(on_count[0].here && on_count[1..].iter().all(|l| !l.here));
    assert_eq!(on_count[1].note, "80 cells, fitted, R-hat 1.02, converged");
    assert_eq!(on_count[2].note, "9 cells, skipped: 9 cells to fit");

    let on_group = nested.links(html, nested.group("C1").ok());
    assert_eq!(
        hrefs(&on_group),
        [Some("fit.html".into()), None, None, None]
    );
    assert!(on_group[1].here && !on_group[0].here);
}

/// `--html` writes the count fit's page and one beside it for each fitted
/// group, each named for its place and carrying the links; a group
/// without a fit gets none.
#[test]
fn a_page_is_written_for_every_fit() {
    let tmp = tempfile::tempdir().unwrap();
    let nested = read(&nested(tmp.path()));
    let html = tmp.path().join("pages").join("fit.html");
    std::fs::create_dir_all(html.parent().unwrap()).unwrap();
    nested.write_html(None, &html).unwrap();
    let count = std::fs::read_to_string(&html).unwrap();
    let group = std::fs::read_to_string(group_page(&html, "C1")).unwrap();
    assert!(count.contains("every cell") && count.contains(r#""href":"fit.C1.html""#));
    assert!(group.contains("within C1") && group.contains(r#""href":"fit.html""#));
    assert!(!group_page(&html, "C2").exists() && !group_page(&html, "bg").exists());
}

/// The list: the count fit, then a line a group. The arrows move within
/// it; Enter opens a fit and the list comes back as it was; a group
/// without a fit says why and stays; `q` leaves with nothing to open.
#[test]
fn the_list_moves_and_opens_only_fits() {
    let tmp = tempfile::tempdir().unwrap();
    let nested = read(&nested(tmp.path()));
    let mut list = List::new(&nested);
    let screen = crate::tui::screen_text(&mut list, 100, 12);
    for shown in [
        "1 of 3 groups fitted",
        "the count fit: gene modules",
        "80 cells   fitted, R-hat 1.02, converged",
        "9 cells   skipped: 9 cells to fit",
        "open the fit",
    ] {
        assert!(screen.contains(shown), "{shown}:\n{screen}");
    }

    press(&mut list, KeyCode::Up);
    assert_eq!(list.at, 0);
    press(&mut list, KeyCode::Down);
    press(&mut list, KeyCode::Down);
    press(&mut list, KeyCode::Enter);
    assert!(!list.done() && list.open.is_none());
    let screen = crate::tui::screen_text(&mut list, 100, 12);
    assert!(
        screen.contains("C2 has no fit to open (skipped"),
        "{screen}"
    );
    // The next key clears what was said.
    press(&mut list, KeyCode::Up);
    assert!(list.status.is_none());
    press(&mut list, KeyCode::Enter);
    assert!(list.done());
    assert_eq!(list.take_open(), Some(1));
    assert!(!list.done() && list.at == 1);
    assert_eq!(nested.row(1).unwrap().group, "C1");
    assert!(nested.row(0).is_none());

    press(&mut list, KeyCode::End);
    assert_eq!(list.at, 3);
    press(&mut list, KeyCode::PageDown);
    assert_eq!(list.at, 3);
    press(&mut list, KeyCode::Home);
    press(&mut list, KeyCode::Enter);
    assert_eq!(list.take_open(), Some(0));
    press(&mut list, KeyCode::Char('q'));
    assert!(list.done());
    assert_eq!(list.take_open(), None);
}
