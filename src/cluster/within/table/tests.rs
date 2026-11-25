//! A cell table's groups in place of a count fit's: read from a TSV or a
//! Parquet file, its columns chosen, its cells matched to the sites', and
//! a run that writes what a count fit's does.

use super::super::{pairing, CountCells, Grouped, Status};
use super::*;
use crate::cluster::read::Cell;
use crate::cluster::tests::parse_args;
use arrow::array::{ArrayRef, Int64Array};
use arrow::record_batch::RecordBatch;
use std::path::PathBuf;
use std::sync::Arc;

fn cell(batch: &str, barcode: &str) -> Cell {
    Cell {
        barcode: barcode.into(),
        batch: batch.into(),
        fit: true,
    }
}

/// `--within path` with `flags` after it, read as a run reads it.
fn named(path: &Path, flags: &[&str]) -> anyhow::Result<CountCells> {
    let path = path.to_string_lossy().into_owned();
    let mut line = vec!["cluster", "--within", path.as_str()];
    line.extend(flags);
    CountCells::named(&path, &parse_args(&line))
}

/// `cells` grouped as a run off a terminal groups them.
fn group(count: &CountCells, cells: &[Cell]) -> anyhow::Result<Grouped> {
    let ours: Vec<&str> = cells.iter().map(|c| &*c.batch).collect();
    let pairing = pairing::resolve(count, &ours, &[], None)?;
    count.group(cells, &pairing)
}

/// Each site cell's group by name, none when it has none.
fn names(count: &CountCells, cells: &[Cell]) -> Vec<Option<String>> {
    let grouped = group(count, cells).unwrap();
    (grouped.of.iter())
        .map(|g| g.map(|g| count.name(g)))
        .collect()
}

fn write_parquet(path: &Path, columns: Vec<(&str, ArrayRef)>) {
    let batch = RecordBatch::try_from_iter(columns).unwrap();
    let file = std::fs::File::create(path).unwrap();
    let mut w = parquet::arrow::ArrowWriter::try_new(file, batch.schema(), None).unwrap();
    w.write(&batch).unwrap();
    w.close().unwrap();
}

fn strings(v: &[&str]) -> ArrayRef {
    Arc::new(StringArray::from(v.to_vec()))
}

/// A TSV, with or without its header, groups its cells as it says: each
/// group named by its value, the groups in order, the cells it names in
/// each counted. Its cells are matched by barcode to sites of one batch.
#[test]
fn a_tsv_groups_cells_as_given() {
    let tmp = tempfile::tempdir().unwrap();
    for (name, text) in [
        ("bare.tsv", "c1\tCT2\nc2\tCT1\nc3\tCT1\n"),
        ("head.tsv", "cell\tgroup\nc1\tCT2\nc2\tCT1\nc3\tCT1\n"),
        ("head.csv", "barcode,group\nc1,CT2\nc2,CT1\nc3,CT1\n"),
    ] {
        let path = tmp.path().join(name);
        std::fs::write(&path, text).unwrap();
        let count = named(&path, &[]).unwrap();
        let t = count.table.as_ref().unwrap();
        assert_eq!(t.names, ["CT1", "CT2"], "{name}");
        assert_eq!(t.named, [2, 1], "{name}");
        let cells = [cell("b1", "c3"), cell("b1", "c1"), cell("b1", "c9")];
        assert_eq!(
            names(&count, &cells),
            [Some("CT1".into()), Some("CT2".into()), None],
            "{name}"
        );
        assert_eq!(count.role(Group(0)), None);
        assert_eq!(count.source(), "cell_table");
    }
    // Without a header, no column can be named.
    let err = named(&tmp.path().join("bare.tsv"), &["--within-column", "group"]).unwrap_err();
    assert!(format!("{err:#}").contains("no header"), "{err:#}");
}

/// A Parquet table's group column is `label` by default, else `cluster`,
/// or the one `--within-column` names; its cell column is found by its
/// usual names or named. One that fits none is refused, naming the
/// columns found. Any type is read as text.
#[test]
fn the_column_choice_works() {
    let tmp = tempfile::tempdir().unwrap();
    let path = tmp.path().join("t.cell_annotation.parquet");
    write_parquet(
        &path,
        vec![
            ("obs_names", strings(&["c1", "c2", "c3"])),
            (
                "cluster",
                Arc::new(Int64Array::from(vec![2, 10, 2])) as ArrayRef,
            ),
            ("label", strings(&["CT1", "CT2", ""])),
        ],
    );
    let by = |flags: &[&str]| {
        let count = named(&path, flags).unwrap();
        let cells = [cell("b", "c1"), cell("b", "c2"), cell("b", "c3")];
        (
            count.table.as_ref().unwrap().names.clone(),
            names(&count, &cells),
        )
    };
    // A cell with an empty label has no group.
    let (groups, of) = by(&[]);
    assert_eq!(groups, ["CT1", "CT2"]);
    assert_eq!(of, [Some("CT1".into()), Some("CT2".into()), None]);
    // Numbers in their order, as text.
    let (groups, of) = by(&["--within-column", "cluster"]);
    assert_eq!(groups, ["2", "10"]);
    assert_eq!(of, [Some("2".into()), Some("10".into()), Some("2".into())]);
    // The cell column named; one that is not there refused.
    let (groups, _) = by(&["--within-cell-column", "OBS_NAMES"]);
    assert_eq!(groups, ["CT1", "CT2"]);
    let err = named(&path, &["--within-column", "kind"]).unwrap_err();
    let err = format!("{err:#}");
    assert!(
        err.contains("no column kind") && err.contains("obs_names, cluster, label"),
        "{err}"
    );
    let other = tmp.path().join("other.parquet");
    write_parquet(
        &other,
        vec![("id", strings(&["c1"])), ("label", strings(&["CT1"]))],
    );
    let err = format!("{:#}", named(&other, &[]).unwrap_err());
    assert!(
        err.contains("none of cell, barcode, cell_barcode, obs_names")
            && err.contains("--within-cell-column")
            && err.contains("its columns: id, label"),
        "{err}"
    );
    let cells = named(&other, &["--within-cell-column", "id"]).unwrap();
    assert_eq!(cells.table.unwrap().names, ["CT1"]);
}

/// `{barcode}@{batch}` is matched in that batch, its matrix kind aside
/// (`b1_count` is `b1`); a bare barcode in each site batch nothing pairs
/// with, left out where it is in more than one.
#[test]
fn a_cells_batch_is_matched_as_a_count_fits_is() {
    let tmp = tempfile::tempdir().unwrap();
    let path = tmp.path().join("t.tsv");
    std::fs::write(
        &path,
        "c1@b1_count\tCT1\nc1@b2_count\tCT2\nc7\tCT1\nc8\tCT2\n",
    )
    .unwrap();
    let count = named(&path, &[]).unwrap();
    let cells = [
        cell("b1", "c1"),
        cell("b2", "c1"),
        cell("b3", "c1"),
        cell("b3", "c7"),
        cell("b4", "c8"),
        cell("b5", "c8"),
    ];
    assert_eq!(
        names(&count, &cells),
        [
            Some("CT1".into()),
            Some("CT2".into()),
            None,
            Some("CT1".into()),
            None,
            None
        ]
    );
    // A batch named with its folder pairs too.
    let cells = [cell("runs/b1", "c1"), cell("b2", "c1")];
    assert_eq!(
        names(&count, &cells),
        [Some("CT1".into()), Some("CT2".into())]
    );
}

/// Names that match none of the site cells are refused, a cell of each
/// side named, with how a table's cells are matched.
#[test]
fn unmatched_names_are_refused_with_examples() {
    let tmp = tempfile::tempdir().unwrap();
    for (name, text, theirs) in [
        ("tagged.tsv", "x1@b1\tCT1\n", "b1 x1"),
        ("bare.tsv", "x1\tCT1\n", "x1 (no batch)"),
    ] {
        let path = tmp.path().join(name);
        std::fs::write(&path, text).unwrap();
        let count = named(&path, &[]).unwrap();
        let err = group(&count, &[cell("b1", "c1"), cell("b2", "c2")])
            .unwrap_err()
            .to_string();
        assert!(
            err.contains("none of the 2 cells")
                && err.contains("the cell table's")
                && err.contains("a site cell is b1 c1")
                && err.contains(theirs)
                && err.contains("{barcode}@{batch}"),
            "{name}: {err}"
        );
    }
}

/// A group's value is made a file stem; two made the same are refused,
/// as is one cell named in two groups.
#[test]
fn group_names_are_file_stems() {
    assert_eq!(sanitised("CT 1/a"), "CT_1_a");
    assert_eq!(sanitised(".x"), "_.x");
    assert_eq!(sanitised("c1-a.b_c"), "c1-a.b_c");
    let tmp = tempfile::tempdir().unwrap();
    let path = tmp.path().join("t.tsv");
    std::fs::write(&path, "c1\tCT 1\nc2\tCT/1\n").unwrap();
    let err = format!("{:#}", named(&path, &[]).unwrap_err());
    assert!(
        err.contains("`CT 1` and `CT/1` are both named CT_1"),
        "{err}"
    );
    std::fs::write(&path, "c1\tCT1\nc1\tCT2\n").unwrap();
    let err = format!("{:#}", named(&path, &[]).unwrap_err());
    assert!(err.contains("c1 is named twice"), "{err}");
    // A table's columns are not taken for a count fit's record.
    let err = named(tmp.path(), &["--within-column", "label"]).unwrap_err();
    assert!(err.to_string().contains("is no table"), "{err}");
}

/// Write `{dir}/{name}.zarr`: `rows` by `cols`, from triplets.
fn write_matrix(dir: &Path, name: &str, rows: &[String], cols: &[&str], e: Vec<(u64, u64, f32)>) {
    use data_beans::sparse_io::{create_sparse_from_triplets_owned, SparseIoBackend};
    let path = dir
        .join(format!("{name}.zarr"))
        .to_string_lossy()
        .into_owned();
    let n = e.len();
    let mut data = create_sparse_from_triplets_owned(
        e,
        (rows.len(), cols.len(), n),
        Some(&path),
        Some(&SparseIoBackend::Zarr),
    )
    .unwrap();
    let rows: Vec<Box<str>> = rows.iter().map(|s| s.as_str().into()).collect();
    let cols: Vec<Box<str>> = cols.iter().map(|s| (*s).into()).collect();
    data.register_row_names_vec(&rows);
    data.register_column_names_vec(&cols);
}

/// A run within a table's groups writes what a run within a count fit's
/// does, `{m}_within.run.json` and `{m}_within_cells.parquet`, a fit under
/// `within/{m}/{group}/` for each group; `cluster-view` opens a group by
/// its name.
#[test]
fn a_run_within_a_tables_groups_writes_the_count_fits_files() {
    let tmp = tempfile::tempdir().unwrap();
    let data = tmp.path().join("data");
    std::fs::create_dir_all(&data).unwrap();
    let cols = ["c1", "c2", "c3", "c4", "c5", "c6", "c7", "c8"];
    let mut rows = Vec::new();
    for s in 0..6 {
        let gene = if s < 3 { "GENE1" } else { "GENE2" };
        for channel in ["methylated", "unmethylated"] {
            rows.push(format!("{gene}/m6a/chr1:{}/{channel}", 100 + s));
        }
    }
    let mut e = Vec::new();
    for c in 0..cols.len() as u64 {
        for s in 0..6u64 {
            let on = ((c + s) % 3 + 1) as f32;
            e.push((2 * s, c, on));
            e.push((2 * s + 1, c, 4.0 - on));
        }
    }
    write_matrix(&data, "b1_m6a_site", &rows, &cols, e);
    let table = tmp.path().join("t.tsv");
    std::fs::write(
        &table,
        "cell\tlabel\nc1@b1_count\tCT1\nc2@b1_count\tCT1\nc3@b1_count\tCT1\nc4@b1_count\tCT1\n\
         c5@b1_count\tCT2\nc6@b1_count\tCT2\nc7@b1_count\tCT2\nc8@b1_count\tCT2\nc9@b1_count\tCT3\n",
    )
    .unwrap();
    let out = tmp.path().join("out");
    let (table_s, data_s, out_s) = (
        table.to_string_lossy().into_owned(),
        data.to_string_lossy().into_owned(),
        out.to_string_lossy().into_owned(),
    );
    let args = parse_args(&[
        "cluster",
        &data_s,
        "--within",
        &table_s,
        "-o",
        &out_s,
        "--within-min-cells",
        "2",
        "--chains",
        "1",
        "--min-cell-sites",
        "1",
        "--min-site-cells",
        "1",
        "--cell-roles",
        "2",
        "--gene-roles",
        "1",
        "--epochs",
        "3",
        "--batch-process",
    ]);
    assert_eq!(super::super::run(&args).unwrap().as_deref(), Some("m6a"));
    let record = out.join("m6a_within.run.json");
    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&record).unwrap()).unwrap();
    let (modality, source, groups) = super::super::listed(&json).unwrap();
    assert_eq!(modality, "m6a");
    assert_eq!(source, super::super::Source::Table(table.clone()));
    let said: Vec<(&str, usize, Status)> = (groups.iter())
        .map(|g| (g.group.as_str(), g.cells, g.status))
        .collect();
    assert_eq!(
        said,
        [
            ("CT1", 4, Status::Fitted),
            ("CT2", 4, Status::Fitted),
            ("CT3", 0, Status::Skipped),
        ]
    );
    for g in ["CT1", "CT2"] {
        assert!(out
            .join(format!("within/m6a/{g}/m6a_blockmodel.run.json"))
            .is_file());
    }
    let cells = read_record_batch(out.join("m6a_within_cells.parquet")).unwrap();
    let group = crate::common::column::<StringArray>(&cells, "group").unwrap();
    let groups: Vec<&str> = group.iter().map(Option::unwrap).collect();
    assert_eq!(
        groups,
        ["CT1", "CT1", "CT1", "CT1", "CT2", "CT2", "CT2", "CT2"]
    );
    let label = crate::common::column::<StringArray>(&cells, "label").unwrap();
    assert!(label.value(0).starts_with("CT1."), "{}", label.value(0));
    // cluster-view: a group by its name, the record whole, a group's page.
    let path: Box<str> = record.to_string_lossy().into();
    let opens = |group, html| crate::cluster::opens(std::slice::from_ref(&path), None, group, html);
    assert!(!opens(Some("ct2"), None).unwrap().whole);
    assert!(opens(None, None).unwrap().whole);
    let html = PathBuf::from("fit.html");
    assert_eq!(
        opens(None, Some(html.as_path())).unwrap().pages,
        [PathBuf::from("fit.CT1.html"), "fit.CT2.html".into()]
    );
    let err = opens(Some("CT3"), None).unwrap_err().to_string();
    assert!(err.contains("CT3 has no fit to open"), "{err}");
}
