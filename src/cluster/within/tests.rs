//! `faba cluster --within` up to its chains: a count fit's cells grouped,
//! the site cells matched to them, the nested cell table, and the record's
//! list of fits. No chain is run here.

use super::*;
use crate::cluster::read::Cell;
use crate::cluster::tests::parse_args;
use crate::cluster::write::{NestedRow, RoleRow};
use serde_json::json;

fn cell(batch: &str, barcode: &str) -> Cell {
    Cell {
        barcode: barcode.into(),
        batch: batch.into(),
        fit: true,
    }
}

/// A count fit in `dir`, its record's path: cells B0..B4 of batch `b`. B0
/// and B1 are in cell module 1, B2 and B3 in module 2, B4 not kept. Its
/// matrices were in the folders `a` (two of them) and `c`.
fn count_fit(dir: &Path) -> PathBuf {
    count_fit_of(dir, &["a/x_count.zarr", "a/y_count.zarr", "c/z_count.zarr"])
}

/// [`count_fit`], its matrices at `matrices`.
fn count_fit_of(dir: &Path, matrices: &[&str]) -> PathBuf {
    let cells: Vec<Cell> = (0..5).map(|i| cell("b", &format!("B{i}"))).collect();
    let member: [Option<[f32; 2]>; 5] = [
        Some([1.0, 0.0]),
        Some([0.9, 0.1]),
        Some([0.0, 1.0]),
        Some([0.2, 0.8]),
        None,
    ];
    let rows: Vec<RoleRow> = (cells.iter().zip(&member))
        .map(|(cell, m)| RoleRow {
            cell,
            membership: m.as_ref().map(|m| m.as_slice()),
            level: None,
        })
        .collect();
    write::cell_roles(&rows, 2, &dir.join("count_cell_roles.parquet")).unwrap();
    let record = json!({
        "job": "cluster",
        "inputs": { "matrices": matrices.iter().map(|m| json!({ "path": m })).collect::<Vec<_>>() },
        "results": { "tables": { "modality": "count", "cells": ["count_cell_roles.parquet"] } },
    });
    let path = dir.join(COUNT_RECORD);
    std::fs::write(&path, record.to_string()).unwrap();
    path
}

/// The count fit's main module and group of one site cell, read alone.
fn of(count: &CountCells, batch: &str, barcode: &str) -> Option<Group> {
    group(count, &[cell(batch, barcode)]).ok()?.of[0]
}

/// `cells` grouped as a run off a terminal groups them: their batches
/// paired by name, refused when that is ambiguous.
fn group(count: &CountCells, cells: &[Cell]) -> anyhow::Result<Grouped> {
    let pairing = pairing_of(count, cells, &[])?;
    count.group(cells, &pairing)
}

/// How the batches of `cells` pair with `count`'s, `named` (as
/// `--batch-map` gives them) first, nothing asked.
fn pairing_of(
    count: &CountCells,
    cells: &[Cell],
    named: &[&str],
) -> anyhow::Result<pairing::Pairing> {
    let ours: Vec<&str> = cells.iter().map(|c| &*c.batch).collect();
    let named: Vec<pairing::Named> = (named.iter())
        .map(|s| pairing::Named::parse(s))
        .collect::<anyhow::Result<_>>()?;
    pairing::resolve(count, &ours, &named, None)
}

/// A cell goes to its main cell module, numbered from 0; one the count
/// fit did not keep has no group.
#[test]
fn a_cell_is_grouped_by_its_main_module() {
    assert_eq!(Group::of(3), Some(Group(2)));
    assert_eq!(Group::of(1), Some(Group(0)));
    assert_eq!(Group::of(0), None);
    assert_eq!(Group(2).name(), "C3");
    assert!(Group(0) < Group(1));
}

/// The count fit is read from its record or its directory: the cells it
/// kept, each in its group, and the folders of its matrices, each once.
#[test]
fn a_count_fit_is_read_from_its_record_or_its_directory() {
    let tmp = tempfile::tempdir().unwrap();
    let record = count_fit(tmp.path());
    for named in [record.as_path(), tmp.path()] {
        let count = CountCells::read(&named.to_string_lossy()).unwrap();
        assert!(count.record.is_absolute() && count.record.ends_with(COUNT_RECORD));
        assert_eq!(of(&count, "b", "B0"), Some(Group(0)));
        assert_eq!(of(&count, "b", "B1"), Some(Group(0)));
        assert_eq!(of(&count, "b", "B2"), Some(Group(1)));
        assert_eq!(of(&count, "b", "B3"), Some(Group(1)));
        assert_eq!(of(&count, "b", "B4"), None);
    }
}

/// The sites are looked for where the count fit's matrices are, each
/// folder once; a matrix that is no longer there (the output was moved)
/// is taken to be beside the record.
#[test]
fn the_sites_default_to_the_folders_of_the_count_fits_matrices() {
    let tmp = tempfile::tempdir().unwrap();
    let at = |name: &str| tmp.path().join(name).to_string_lossy().into_owned();
    for name in ["a/x_count.zarr", "a/y_count.zarr", "c/z_count.zarr"] {
        std::fs::create_dir_all(at(name)).unwrap();
    }
    let here = tmp.path().join("fit");
    std::fs::create_dir_all(&here).unwrap();
    let matrices = [
        at("a/x_count.zarr"),
        at("a/y_count.zarr"),
        at("c/z_count.zarr"),
    ];
    let matrices: Vec<&str> = matrices.iter().map(String::as_str).collect();
    let record = count_fit_of(&here, &matrices);
    let count = CountCells::read(&record.to_string_lossy()).unwrap();
    assert_eq!(
        count.folders(),
        [at("a").into(), at("c").into()] as [Box<str>; 2]
    );
    let record = count_fit_of(&here, &["/gone/a/x_count.zarr", "/gone/a/y_count.zarr"]);
    let count = CountCells::read(&record.to_string_lossy()).unwrap();
    let here: Box<str> = here.to_string_lossy().into();
    assert_eq!(count.folders(), [here]);
}

/// What is not a count fit's record is refused, saying what to run first:
/// a directory without one, and the record of a site fit.
#[test]
fn what_is_not_a_count_fit_is_refused() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().to_string_lossy().into_owned();
    let err = CountCells::read(&dir).unwrap_err().to_string();
    assert!(err.contains("no count_blockmodel.run.json there"), "{err}");
    let site = tmp.path().join("m6a_blockmodel.run.json");
    let fit = json!({ "job": "cluster", "results": { "tables": { "modality": "m6a" } } });
    std::fs::write(&site, fit.to_string()).unwrap();
    let err = CountCells::read(&site.to_string_lossy()).unwrap_err();
    let err = err.to_string();
    assert!(
        err.contains("is not a count fit's run record")
            && err.contains("faba cluster --modality count"),
        "{err}"
    );
}

/// The site cells are grouped in the order read; those the count fit has
/// no group for are counted. A count fit of one batch matches by barcode
/// alone, the sites' batch named otherwise; with none matched it is
/// refused, a cell of each named.
#[test]
fn site_cells_are_matched_to_the_count_fits_by_batch_and_barcode() {
    let tmp = tempfile::tempdir().unwrap();
    let record = count_fit(tmp.path());
    let count = CountCells::read(&record.to_string_lossy()).unwrap();
    let cells: Vec<Cell> = ["B3", "B1", "B9", "B2", "B0", "B4"]
        .iter()
        .map(|b| cell("sites", b))
        .collect();
    let grouped = group(&count, &cells).unwrap();
    let expected: BTreeMap<Group, Vec<usize>> =
        [(Group(0), vec![1, 4]), (Group(1), vec![0, 3])].into();
    assert_eq!(grouped.groups, expected);
    assert_eq!((grouped.matched(), grouped.of.len()), (4, 6));
    assert_eq!(grouped.of[0], Some(Group(1)));
    assert_eq!(grouped.of[2], None);
    let err = group(&count, &[cell("sites", "X1")])
        .unwrap_err()
        .to_string();
    assert!(
        err.contains("none of the 1 cells") && err.contains("sites X1") && err.contains("b B"),
        "{err}"
    );
}

/// With several batches in the count fit, a site cell is matched in its
/// own batch only.
#[test]
fn a_count_fit_of_several_batches_matches_within_the_batch() {
    let tmp = tempfile::tempdir().unwrap();
    let record = count_fit(tmp.path());
    let cells = [cell("b1", "B0"), cell("b2", "B0")];
    let member = [[1.0f32, 0.0], [0.0, 1.0]];
    let rows: Vec<RoleRow> = (cells.iter().zip(&member))
        .map(|(cell, m)| RoleRow {
            cell,
            membership: Some(m.as_slice()),
            level: None,
        })
        .collect();
    write::cell_roles(&rows, 2, &tmp.path().join("count_cell_roles.parquet")).unwrap();
    let count = CountCells::read(&record.to_string_lossy()).unwrap();
    assert_eq!(of(&count, "b1", "B0"), Some(Group(0)));
    assert_eq!(of(&count, "b2", "B0"), Some(Group(1)));
    assert_eq!(of(&count, "b3", "B0"), None);
    // Sites of several batches are not matched by barcode alone either.
    let record = count_fit(tmp.path());
    let count = CountCells::read(&record.to_string_lossy()).unwrap();
    let grouped = group(&count, &[cell("b", "B0"), cell("other", "B0")]).unwrap();
    assert_eq!(grouped.of, [Some(Group(0)), None]);
}

/// The count fit of [`count_fit`] with cells of three batches, `a/b`,
/// `c/b` and `solo`: B0 of each, in module 1, 2 and 1.
fn three_batches(tmp: &Path) -> CountCells {
    let record = count_fit(tmp);
    let cells = [cell("a/b", "B0"), cell("c/b", "B0"), cell("solo", "B0")];
    let member = [[1.0f32, 0.0], [0.0, 1.0], [1.0, 0.0]];
    let rows: Vec<RoleRow> = (cells.iter().zip(&member))
        .map(|(cell, m)| RoleRow {
            cell,
            membership: Some(m.as_slice()),
            level: None,
        })
        .collect();
    write::cell_roles(&rows, 2, &tmp.join("count_cell_roles.parquet")).unwrap();
    CountCells::read(&record.to_string_lossy()).unwrap()
}

/// A batch of the same name wins over one named with its folder: `solo`
/// is the count fit's `solo`, and `/x/solo`, which could be it too, is
/// left unmatched, said why.
#[test]
fn a_batch_of_the_same_name_wins() {
    let tmp = tempfile::tempdir().unwrap();
    let count = three_batches(tmp.path());
    let sites = [cell("/x/solo", "B0"), cell("solo", "B0")];
    let p = pairing_of(&count, &sites, &[]).unwrap();
    assert_eq!(
        p.pairs,
        [pairing::Pair {
            sites: "solo".into(),
            count: "solo".into(),
            by: pairing::By::Same
        }]
    );
    assert_eq!(p.left.len(), 1);
    assert_eq!(&*p.left[0].0, "/x/solo");
    assert!(
        p.left[0]
            .1
            .contains("could be the count fit's solo, which another batch has"),
        "{:?}",
        p.left
    );
    assert_eq!(group(&count, &sites).unwrap().of, [None, Some(Group(0))]);
    let recorded = p.recorded();
    assert_eq!(recorded["pairs"][0]["by"], "same");
    assert_eq!(recorded["unmatched"][0]["sites"], "/x/solo");
}

/// `--batch-map` pairs come before any by name: a pairing that would be
/// ambiguous is made, a batch is left unmatched with `SITES=`, and a pair
/// naming a batch neither side has is refused, as is one not `SITES=COUNT`.
#[test]
fn batch_map_pairs_come_first() {
    let tmp = tempfile::tempdir().unwrap();
    let count = three_batches(tmp.path());
    let sites = [
        cell("/x/solo", "B0"),
        cell("/y/solo", "B0"),
        cell("b", "B0"),
    ];
    let named = ["/y/solo=solo", "b=c/b"];
    let p = pairing_of(&count, &sites, &named).unwrap();
    let pairs: Vec<(&str, &str, pairing::By)> = (p.pairs.iter())
        .map(|p| (&*p.sites, &*p.count, p.by))
        .collect();
    assert_eq!(
        pairs,
        [
            ("/y/solo", "solo", pairing::By::Named),
            ("b", "c/b", pairing::By::Named)
        ]
    );
    assert!(p.ambiguous.is_empty());
    let grouped = count.group(&sites, &p).unwrap();
    assert_eq!(grouped.of, [None, Some(Group(0)), Some(Group(1))]);
    // `SITES=` leaves a batch unmatched; the other then pairs alone.
    let p = pairing_of(&count, &sites, &["/x/solo=", "b=a/b"]).unwrap();
    let left: Vec<&str> = p.left.iter().map(|l| &*l.0).collect();
    assert_eq!(left, ["/x/solo"]);
    assert!(p
        .pairs
        .iter()
        .any(|p| &*p.sites == "/y/solo" && p.by == pairing::By::Folder));
    for (named, says) in [
        ("nope=solo", "the sites have no batch nope"),
        ("b=nope", "the count fit has no batch nope"),
        ("b", "write SITES=COUNT"),
    ] {
        let err = pairing_of(&count, &sites, &[named])
            .unwrap_err()
            .to_string();
        assert!(err.contains(says), "{named}: {err}");
    }
}

/// What is still ambiguous is refused off a terminal, each choice named
/// and how to make it, on the command line and through `faba mcp`.
#[test]
fn an_ambiguous_pairing_is_refused_naming_each_choice() {
    let tmp = tempfile::tempdir().unwrap();
    let count = three_batches(tmp.path());
    let sites = [
        cell("/x/solo", "B0"),
        cell("/y/solo", "B0"),
        cell("b", "B0"),
    ];
    let err = group(&count, &sites).unwrap_err().to_string();
    for says in [
        "pair with the count fit's in more than one way",
        "the sites' batch `b` could be the count fit's `a/b` or `c/b`",
        "the sites' batches `/x/solo` and `/y/solo` could each be the count fit's `solo`",
        "--batch-map SITES=COUNT",
        "--batch-map '/x/solo=solo'",
        "`batch_map`, e.g. {\"/x/solo\": \"solo\"}",
    ] {
        assert!(err.contains(says), "{says}: {err}");
    }
}

/// On a terminal each ambiguity is asked in turn: a choice is taken as if
/// named, `neither` leaves its site batches unmatched, and a question
/// cancelled cancels the run.
#[test]
fn an_ambiguous_pairing_is_asked_one_at_a_time() {
    let tmp = tempfile::tempdir().unwrap();
    let count = three_batches(tmp.path());
    let ours = ["/x/solo", "/y/solo", "b"];
    let mut asked: Vec<String> = Vec::new();
    let mut answers = vec![pairing::Answer::Neither, pairing::Answer::Take(1)];
    let mut ask = |a: &pairing::Ambiguity| {
        asked.push(a.said());
        Ok(answers.pop())
    };
    let p = pairing::resolve(&count, &ours, &[], Some(&mut ask)).unwrap();
    assert_eq!(asked.len(), 2, "{asked:?}");
    assert!(asked[0].contains("`/x/solo` and `/y/solo`"), "{asked:?}");
    assert!(asked[1].contains("`b` could be"), "{asked:?}");
    // The second choice of the first, `neither` of the second.
    let pairs: Vec<(&str, &str)> = p.pairs.iter().map(|p| (&*p.sites, &*p.count)).collect();
    assert_eq!(pairs, [("/y/solo", "solo")]);
    let left: Vec<&str> = p.left.iter().map(|l| &*l.0).collect();
    assert_eq!(left, ["b", "/x/solo"]);
    let mut cancel = |_: &pairing::Ambiguity| Ok(None);
    let err = pairing::resolve(&count, &ours, &[], Some(&mut cancel)).unwrap_err();
    assert!(err.is::<progress::Cancelled>(), "{err}");
}

/// `--batch-map` is taken once a pair, with `--within` only.
#[test]
fn batch_map_is_repeated_with_within() {
    let args = parse_args(&[
        "cluster",
        "--within",
        "c.run.json",
        "--batch-map",
        "/x/solo=solo",
        "--batch-map=b=",
    ]);
    assert_eq!(
        args.batch_map,
        ["/x/solo=solo".into(), "b=".into()] as [Box<str>; 2]
    );
    use clap::Args;
    let cmd = ClusterArgs::augment_args(clap::Command::new("cluster"));
    let err = cmd
        .try_get_matches_from(["cluster", "d", "--batch-map", "b=c"])
        .unwrap_err()
        .to_string();
    assert!(err.contains("--within"), "{err}");
}

/// The pop-up lists each pairing and `neither`; ↓ and Enter take one, the
/// last row leaves them unmatched, Esc cancels.
#[test]
fn the_question_takes_a_pairing_or_neither() {
    use ratatui::crossterm::event::{KeyCode, KeyEvent, KeyModifiers};
    let a = pairing::Ambiguity {
        choices: vec![("b".into(), "a/b".into()), ("b".into(), "c/b".into())],
    };
    let press = |q: &mut pairing::Question, code| {
        crate::tui::feed(q, KeyEvent::new(code, KeyModifiers::NONE))
    };
    let mut q = pairing::Question::new(&a);
    let screen = crate::tui::screen_text(&mut q, 100, 20);
    for says in [
        "could be the count fit's `a/b` or `c/b`",
        "count fit's c/b",
        "neither",
    ] {
        assert!(screen.contains(says), "{says}:\n{screen}");
    }
    press(&mut q, KeyCode::Down);
    press(&mut q, KeyCode::Enter);
    assert_eq!(q.decision, Some(Some(pairing::Answer::Take(1))));
    let mut q = pairing::Question::new(&a);
    (0..5).for_each(|_| press(&mut q, KeyCode::Down));
    press(&mut q, KeyCode::Enter);
    assert_eq!(q.decision, Some(Some(pairing::Answer::Neither)));
    let mut q = pairing::Question::new(&a);
    press(&mut q, KeyCode::Esc);
    assert_eq!(q.decision, Some(None));
    let mut q = pairing::Question::new(&a);
    let keys = crate::tui::help::shown(&mut q, 100, 30).expect("? opens the keys");
    for key in ["↑/↓", "PgUp/PgDn", "Home/End", "Enter", "Esc/q"] {
        assert!(keys.contains(key), "{key}:\n{keys}");
    }
    assert_eq!(q.decision, None, "? chooses nothing");
}

/// A batch met in several folders is named with its folder, as typed in
/// each round: `a/b` in the count fit is `/runs/a/b` of the sites. A name
/// that could be either of two batches is not taken as either.
#[test]
fn a_batch_named_with_its_folder_is_matched_however_the_folder_was_typed() {
    let tmp = tempfile::tempdir().unwrap();
    let count = three_batches(tmp.path());
    let sites = [
        cell("/runs/a/b", "B0"),
        cell("/runs/c/b", "B0"),
        cell("/runs/a/solo", "B0"),
        cell("xa/b", "B0"),
    ];
    let grouped = group(&count, &sites).unwrap();
    let (first, second) = (Some(Group(0)), Some(Group(1)));
    assert_eq!(grouped.of, [first, second, first, None]);
    // `b` alone is either folder's: not taken as either.
    assert_eq!(of(&count, "b", "B0"), None);
}

/// The nested table has a row a cell read: its count module, its group,
/// its site module within the group, and the label of the two; a cell
/// without a group has none, one whose group was not fitted the group's.
#[test]
fn the_nested_table_labels_each_cell_by_group_and_site_module() {
    let cells = [cell("b", "B0"), cell("b", "B1"), cell("b", "B2")];
    let rows = [
        NestedRow {
            cell: &cells[0],
            count_role: Some(3),
            group: Some("C3".into()),
            site_role: Some(2),
        },
        NestedRow {
            cell: &cells[1],
            count_role: Some(1),
            group: Some("C1".into()),
            site_role: None,
        },
        NestedRow {
            cell: &cells[2],
            count_role: None,
            group: None,
            site_role: None,
        },
    ];
    let labels: Vec<Option<String>> = rows.iter().map(NestedRow::label).collect();
    assert_eq!(labels, [Some("C3.2".into()), Some("C1".into()), None]);
    let tmp = tempfile::tempdir().unwrap();
    let path = tmp.path().join("m6a_within_cells.parquet");
    write::nested_cells(&rows, &path).unwrap();
    let t = read_record_batch(&path).unwrap();
    assert_eq!(t.num_rows(), 3);
    let label = column::<StringArray>(&t, "label").unwrap();
    let read: Vec<Option<&str>> = label.iter().collect();
    assert_eq!(read, [Some("C3.2"), Some("C1"), None]);
    let role = column::<Int32Array>(&t, "site_role").unwrap();
    assert_eq!(role.iter().collect::<Vec<_>>(), [Some(2), None, None]);
    let count = column::<Int32Array>(&t, "count_role").unwrap();
    assert_eq!(count.iter().collect::<Vec<_>>(), [Some(3), Some(1), None]);
    // A group's cells get their site modules from its cell table, in order.
    assert!(main_roles(&path, 3).is_err(), "no `role` column there");
}

/// A group's fit is listed in the record as it is read back, and said in
/// words: its status, its chains' R-hat and verdict, why it has no fit.
#[test]
fn a_groups_fit_is_listed_and_said_in_words() {
    let mut fitted = GroupFit::new(Group(2), 120, Status::Fitted);
    fitted.dir = Some("within/m6a/C3".into());
    fitted.record = Some("m6a_blockmodel.run.json".into());
    fitted.converged = Some(true);
    let mut early = fitted.clone();
    early.converged = Some(false);
    early.stopped_early = true;
    let mut skipped = GroupFit::new(Group(4), 7, Status::Skipped);
    skipped.why = Some("7 cells to fit, fewer than --within-min-cells 50".into());
    let waiting = GroupFit::new(Group(3), 60, Status::NotRun);
    assert_eq!(fitted.what(), "fitted, converged");
    assert_eq!(fitted.line(), "C3\t120 cells\tfitted, converged");
    assert_eq!(early.what(), "fitted, not converged, stopped early");
    assert_eq!(
        skipped.what(),
        "skipped: 7 cells to fit, fewer than --within-min-cells 50"
    );
    assert_eq!(waiting.what(), "not run");

    let fits = vec![fitted, early, skipped, waiting];
    let record = json!({ "job": "cluster", "results": { "within": {
        "modality": "m6a", "count_record": "/fits/count_blockmodel.run.json", "groups": fits,
    } } });
    assert_eq!(
        record["results"]["within"]["groups"][3]["status"],
        "not_run"
    );
    assert!(record["results"]["within"]["groups"][2]
        .get("dir")
        .is_none());
    let (modality, count, read) = listed(&record).unwrap();
    assert_eq!(modality, "m6a");
    assert_eq!(count, Path::new("/fits/count_blockmodel.run.json"));
    assert_eq!(read, fits);
    // Any other record lists none.
    assert!(listed(&json!({ "job": "cluster", "results": { "model": {} } })).is_none());
}

/// `--within` fits sites with the blockmodel, within groups worth
/// fitting: a count modality, the mixture and a group of one cell are
/// refused before anything is read.
#[test]
fn settings_within_cannot_fit_with_are_refused() {
    let check = |words: &[&str]| {
        let line = ["cluster", "--within", "count_blockmodel.run.json"];
        let all: Vec<&str> = line.iter().chain(words).copied().collect();
        args::check(&parse_args(&all)).map_err(|e| e.to_string())
    };
    assert!(check(&[]).is_ok());
    assert!(check(&["--model", "blocks"]).is_ok());
    for (words, says) in [
        (&["--modality", "count"][..], "not count"),
        (&["--within-min-cells", "1"], "at least 2"),
    ] {
        let err = check(words).unwrap_err();
        assert!(err.contains(says), "{words:?}: {err}");
    }
}

/// Two `--batch-map` pairs may not name one count batch: refused, not two
/// site batches grouped against it.
#[test]
fn a_count_batch_named_twice_is_refused() {
    let tmp = tempfile::tempdir().unwrap();
    let count = three_batches(tmp.path());
    let sites = [cell("x", "B0"), cell("y", "B0")];
    let err = pairing_of(&count, &sites, &["x=solo", "y=solo"]).unwrap_err();
    assert!(
        format!("{err:#}").contains("names the count fit's batch solo twice"),
        "{err:#}"
    );
}

/// A site batch whose same-named count batch `--batch-map` gave to another
/// is left unmatched and said why, not dropped without a word.
#[test]
fn a_same_name_taken_by_the_map_is_said() {
    let tmp = tempfile::tempdir().unwrap();
    let count = three_batches(tmp.path());
    let sites = [cell("/x/solo", "B0"), cell("solo", "B0")];
    let p = pairing_of(&count, &sites, &["/x/solo=solo"]).unwrap();
    assert!(
        p.left
            .iter()
            .any(|(s, why)| &**s == "solo" && why.contains("given to another batch")),
        "{:?}",
        p.left
    );
}

/// A site batch with `=` in its name splits at the last one.
#[test]
fn batch_map_splits_at_the_last_equals() {
    let n = pairing::Named::parse("run=1=b1").unwrap();
    assert_eq!((&*n.sites, n.count.as_deref()), ("run=1", Some("b1")));
}

/// Answers given when asked are kept, as the flags that give them again.
#[test]
fn answers_are_kept_as_flags() {
    let tmp = tempfile::tempdir().unwrap();
    let count = three_batches(tmp.path());
    let ours = ["b"];
    let mut ask = |_: &pairing::Ambiguity| Ok(Some(pairing::Answer::Take(0)));
    let p = pairing::resolve(&count, &ours, &[], Some(&mut ask)).unwrap();
    assert_eq!(p.answered, ["b=a/b"]);
}
