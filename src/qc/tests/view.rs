//! `faba qc-view` on a synthetic qc output: its kept and dropped sites
//! read back together at the thresholds recorded, each batch's cell QC off
//! its report; a plain directory with none; the page and the figures
//! never over a file.

use super::*;
use crate::qc::sites::{read_site_table, write_site_tables};
use arrow::array::{ArrayRef, Float32Array, Int64Array, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use std::sync::Arc;

/// `n` m6A sites over three genes on chr1: p-values from 0.5 down, a
/// spread of depths, a clean control every fourth.
fn sites(n: usize) -> RecordBatch {
    let i: Vec<usize> = (0..n).collect();
    let pos: Vec<i64> = i.iter().map(|&i| 1000 + 10 * i as i64).collect();
    let cov: Vec<u64> = i.iter().map(|&i| 1 + (i * 7 % 30) as u64).collect();
    let conv: Vec<u64> = (i.iter().zip(&cov))
        .map(|(&i, &c)| (i as u64 % 5).min(c))
        .collect();
    let ctl_conv: Vec<u64> = i.iter().map(|&i| u64::from(i % 4 != 0)).collect();
    let pv: Vec<f32> = i.iter().map(|&i| 0.5 / (1 + i * i) as f32).collect();
    let text = |v: Vec<String>| Arc::new(StringArray::from(v)) as ArrayRef;
    let schema = Schema::new(vec![
        Field::new("chr", DataType::Utf8, false),
        Field::new("gene", DataType::Utf8, false),
        Field::new("strand", DataType::Utf8, false),
        Field::new("primary_pos", DataType::Int64, false),
        Field::new("conversion_pos", DataType::Int64, true),
        Field::new("pv", DataType::Float32, false),
        Field::new("coverage", DataType::UInt64, false),
        Field::new("converted", DataType::UInt64, false),
        Field::new("control_coverage", DataType::UInt64, false),
        Field::new("control_converted", DataType::UInt64, false),
    ]);
    let columns: Vec<ArrayRef> = vec![
        text(i.iter().map(|_| "chr1".into()).collect()),
        text(
            i.iter()
                .map(|i| format!("ENSG{}_GENE{}", i % 3, i % 3))
                .collect(),
        ),
        text(i.iter().map(|_| "+".into()).collect()),
        Arc::new(Int64Array::from(
            pos.iter().map(|p| p - 1).collect::<Vec<_>>(),
        )),
        Arc::new(Int64Array::from(pos)),
        Arc::new(Float32Array::from(pv)),
        Arc::new(UInt64Array::from(cov)),
        Arc::new(UInt64Array::from(conv)),
        Arc::new(UInt64Array::from(vec![10; n])),
        Arc::new(UInt64Array::from(ctl_conv)),
    ];
    RecordBatch::try_new(Arc::new(schema), columns).unwrap()
}

/// What `faba qc` applied in the synthetic run: a real, a count, and the
/// log odds off (an infinity, which the record writes as null).
fn applied() -> SiteFilterArgs {
    let mut f = SiteFilterArgs::parse_flags([
        "--site-max-pv",
        "0.05",
        "--site-min-coverage",
        "3",
        "--site-min-edit-ratio",
        "0.1",
    ]);
    Criterion::MinLogOdds.set(&mut f, f64::NEG_INFINITY);
    f
}

/// One batch's cell QC report: one kept cell, one near empty, one outlier.
const REPORT: &str = "#cell\tn_genes\ttotal_counts\tnear_empty\ttrain_keep\n\
                      c1\t40\t900\t0\t1\nc2\t1\t2\t1\t1\nc3\t80\t9000\t0\t0\n";

/// A qc output in `dir` of `n` sites cut at [`applied`]: the kept and
/// dropped tables, the record and a batch's cell QC report; the sites
/// kept and dropped.
fn qc_output(dir: &Path, n: usize) -> (usize, usize) {
    let f = applied();
    let t = site_table_of(sites(n), "m6a").unwrap();
    let reasons = f.reasons(&t, None);
    let at = |name: &str| dir.join(name).to_string_lossy().into_owned();
    let cut = write_site_tables(
        &t,
        &reasons,
        &at("m6a_sites.parquet"),
        &at("m6a_sites_dropped.parquet"),
    )
    .unwrap();
    let record = serde_json::json!({
        "job": "qc", "status": "ok",
        "inputs": { "fileset": { "path": at("no-such-input") } },
        "options": { "site": f, "qc_mads": 4.0, "no_cell_qc": false },
    });
    std::fs::write(dir.join(RECORD), record.to_string()).unwrap();
    std::fs::write(dir.join("b1_cell_qc_report.tsv"), REPORT).unwrap();
    cut
}

/// The counts table's `sites` rows, as `(flag, removes, only, kept, of)`.
fn site_rows(table: &str) -> Vec<(String, usize, usize, usize, usize)> {
    let n = |s: &str| s.parse::<usize>().unwrap();
    (table.lines())
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|r| r[0] == "sites")
        .map(|r| (r[2].to_string(), n(r[4]), n(r[5]), n(r[6]), n(r[7])))
        .collect()
}

#[test]
fn a_qc_output_opens_at_the_thresholds_recorded_with_its_dropped_sites() {
    let dir = tempfile::tempdir().unwrap();
    let (n_kept, n_dropped) = qc_output(dir.path(), 60);
    assert!(
        n_kept > 0 && n_dropped > 0,
        "{n_kept} kept, {n_dropped} dropped"
    );
    let o = open(&dir.path().to_string_lossy(), None).unwrap();
    assert!(o.applied);
    let shown = |f: &SiteFilterArgs| serde_json::to_value(f).unwrap();
    assert_eq!(shown(&o.filter), shown(&applied()));
    assert!(Criterion::MinLogOdds.is_off(&o.filter), "null reads as off");
    assert_eq!(o.tables["m6a"].len(), n_kept + n_dropped);
    assert_eq!(o.dropped["m6a"], n_dropped);
    assert_eq!(o.cells.len(), 1);
    assert_eq!(o.cells[0].tally(), [1, 1, 1]);
    assert_eq!(o.cell_rule, "4 MADs");
    let batch = Batch::new(
        &o.dir,
        &o.tables,
        &o.site_cells,
        &o.filter,
        (None, &[]),
        (o.cells.clone(), None),
    )
    .unwrap();
    let kept: Vec<(&str, usize, usize)> = vec![("b1", 1, 3)];
    let table = counts_table((&o.cell_rule, &kept), &batch);
    assert!(
        table.contains("cells\tb1\tcell QC\t4 MADs\t2\t\t1\t3"),
        "{table}"
    );
    let rows = site_rows(&table);
    assert!(!rows.is_empty(), "{table}");
    for (flag, _, _, kept, of) in &rows {
        assert_eq!((*kept, *of), (n_kept, n_kept + n_dropped), "{flag}");
    }
    // Every dropped site fails one threshold or more, as recorded.
    let pv = rows.iter().find(|r| r.0 == "--site-max-pv").unwrap();
    let t = read_site_table(
        &dir.path()
            .join("m6a_sites_dropped.parquet")
            .to_string_lossy(),
        "m6a",
    )
    .unwrap();
    let failing_pv = (0..t.len()).filter(|&i| t.pv[i] > 0.05).count();
    assert_eq!(pv.1, failing_pv);
}

/// A real threshold is read back from the record's file as the `f32` it
/// was, not rounded to keep a site at it.
#[test]
fn a_recorded_real_reads_back_exactly() {
    let f = applied();
    let written = serde_json::json!({ "options": { "site": f } }).to_string();
    let record: Value = serde_json::from_str(&written).unwrap();
    let read = recorded_thresholds(&record).unwrap();
    assert_eq!(read.site_min_edit_ratio, 0.1f32);
    assert_eq!(read.site_max_pv, 0.05f32);
    assert_eq!(read.site_min_coverage, 3);
    let err = recorded_thresholds(&serde_json::json!({ "options": {} })).unwrap_err();
    assert!(err.to_string().contains("no site thresholds"), "{err}");
}

/// A directory no qc run wrote opens with no thresholds: every site kept.
#[test]
fn a_plain_directory_opens_with_nothing_removed() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("m6a_sites.parquet");
    crate::editing::io::write_record_batch(&sites(40), &path).unwrap();
    let o = open(&dir.path().to_string_lossy(), None).unwrap();
    assert!(!o.applied && o.dropped.is_empty() && o.cells.is_empty());
    let shown = |f: &SiteFilterArgs| serde_json::to_value(f).unwrap();
    assert_eq!(shown(&o.filter), shown(&SiteFilterArgs::permissive()));
    let batch = Batch::new(
        &o.dir,
        &o.tables,
        &o.site_cells,
        &o.filter,
        (None, &[]),
        (o.cells.clone(), None),
    )
    .unwrap();
    for (flag, removes, _, kept, of) in site_rows(&counts_table(("", &[]), &batch)) {
        assert_eq!((removes, kept, of), (0, 40, 40), "{flag}");
    }
}

fn parse(line: &str) -> QcViewArgs {
    #[derive(clap::Parser)]
    struct Cli {
        #[command(flatten)]
        view: QcViewArgs,
    }
    let words = line.split(' ').skip(2); // `faba qc-view`
    <Cli as clap::Parser>::parse_from(std::iter::once("faba").chain(words)).view
}

/// `--html` writes the page instead of the view, never over a file: one
/// that exists is refused before anything is read (the input here does
/// not exist) and left as it was.
#[test]
fn html_writes_the_page_and_never_over_a_file() {
    let dir = tempfile::tempdir().unwrap();
    let taken = dir.path().join("taken.html");
    std::fs::write(&taken, "mine").unwrap();
    let args = parse(&format!(
        "faba qc-view /no/such/dir --html {}",
        taken.display()
    ));
    let err = run_qc_view(&args).unwrap_err().to_string();
    assert!(
        err.contains("already exists") && err.contains("another path"),
        "{err}"
    );
    assert_eq!(std::fs::read_to_string(&taken).unwrap(), "mine");
    assert!(args.view_session().is_none(), "--html opens no view");
    let out = dir.path().join("out");
    std::fs::create_dir(&out).unwrap();
    qc_output(&out, 30);
    let page = dir.path().join("pages/qc.html");
    let args = parse(&format!(
        "faba qc-view {} --html {}",
        out.display(),
        page.display()
    ));
    run_qc_view(&args).unwrap();
    let html = std::fs::read_to_string(&page).unwrap();
    assert!(
        html.contains("<html") || html.contains("<!doctype"),
        "a page"
    );
    assert!(run_qc_view(&args)
        .unwrap_err()
        .to_string()
        .contains("already exists"));
}

/// `--plots` writes the figures under -o DIR/qc_plots, DIR new or empty,
/// the site panels under qc_plots/panels; a DIR with files is refused
/// before anything is read.
#[test]
fn plots_go_under_a_new_folder_only() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("out");
    std::fs::create_dir(&out).unwrap();
    qc_output(&out, 30);
    let figs = dir.path().join("figs");
    let args = parse(&format!(
        "faba qc-view {} --plots -o {} --hist pv",
        out.display(),
        figs.display()
    ));
    assert!(args.view_session().is_none(), "--plots opens no view");
    run_qc_view(&args).unwrap();
    let plots = figs.join("qc_plots");
    assert!(plots.join("b1_cells.png").is_file());
    assert!(plots.join("panels/m6a_hist_pv.png").is_file());
    let again = run_qc_view(&args).unwrap_err().to_string();
    assert!(again.contains("never writes over files"), "{again}");
    let view = parse(&format!("faba qc-view {}", out.display()));
    assert!(view.view_session().is_some(), "the view, by default");
}
