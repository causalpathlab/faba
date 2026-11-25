//! End-to-end checks of `faba qc --batch-process` on a synthetic output
//! directory: a count matrix, an m6A site matrix, a gene-level m6A matrix and
//! an m6A site table, all in the producers' layout. Drives the built binary.

use std::process::Command;
use std::sync::Arc;

use arrow::array::{Array, ArrayRef, Float32Array, Int64Array, StringArray, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use data_beans::sparse_io::{
    create_sparse_from_triplets_owned, open_sparse_matrix, SparseIoBackend,
};
use parquet::arrow::arrow_reader::ParquetRecordBatchReaderBuilder;
use parquet::arrow::ArrowWriter;

fn faba() -> Command {
    Command::new(env!("CARGO_BIN_EXE_faba"))
}

/// Write `{dir}/{name}.zarr` from named triplets.
fn write_matrix(
    dir: &str,
    name: &str,
    rows: &[&str],
    cols: &[&str],
    entries: &[(usize, usize, f32)],
) {
    let path = format!("{dir}/{name}.zarr");
    let triplets: Vec<(u64, u64, f32)> = entries
        .iter()
        .map(|&(r, c, x)| (r as u64, c as u64, x))
        .collect();
    let mut data = create_sparse_from_triplets_owned(
        triplets,
        (rows.len(), cols.len(), entries.len()),
        Some(&path),
        Some(&SparseIoBackend::Zarr),
    )
    .unwrap();
    let rows: Vec<Box<str>> = rows.iter().map(|s| (*s).into()).collect();
    let cols: Vec<Box<str>> = cols.iter().map(|s| (*s).into()).collect();
    data.register_row_names_vec(&rows);
    data.register_column_names_vec(&cols);
}

struct Site {
    chr: &'static str,
    gene: &'static str,
    pos: i64,
    pv: f32,
    log_odds: Option<f32>,
    cov: u64,
    conv: u64,
    ctl_cov: u64,
    ctl_conv: u64,
}

/// The m6A site table with exactly the columns `faba qc` reads.
fn write_sites(path: &str, sites: &[Site]) {
    let schema = Arc::new(Schema::new(vec![
        Field::new("chr", DataType::Utf8, false),
        Field::new("gene", DataType::Utf8, false),
        Field::new("primary_pos", DataType::Int64, false),
        Field::new("conversion_pos", DataType::Int64, true),
        Field::new("pv", DataType::Float32, false),
        Field::new("log_odds", DataType::Float32, true),
        Field::new("coverage", DataType::UInt64, false),
        Field::new("converted", DataType::UInt64, false),
        Field::new("control_coverage", DataType::UInt64, false),
        Field::new("control_converted", DataType::UInt64, false),
    ]));
    let batch = RecordBatch::try_new(
        schema.clone(),
        vec![
            Arc::new(StringArray::from(
                sites.iter().map(|s| s.chr).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(StringArray::from(
                sites.iter().map(|s| s.gene).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int64Array::from(
                sites.iter().map(|s| s.pos - 1).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Int64Array::from(
                sites.iter().map(|s| Some(s.pos)).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Float32Array::from(
                sites.iter().map(|s| s.pv).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(Float32Array::from(
                sites.iter().map(|s| s.log_odds).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(UInt64Array::from(
                sites.iter().map(|s| s.cov).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(UInt64Array::from(
                sites.iter().map(|s| s.conv).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(UInt64Array::from(
                sites.iter().map(|s| s.ctl_cov).collect::<Vec<_>>(),
            )) as ArrayRef,
            Arc::new(UInt64Array::from(
                sites.iter().map(|s| s.ctl_conv).collect::<Vec<_>>(),
            )) as ArrayRef,
        ],
    )
    .unwrap();
    let file = std::fs::File::create(path).unwrap();
    let mut w = ArrowWriter::try_new(file, schema, None).unwrap();
    w.write(&batch).unwrap();
    w.close().unwrap();
}

fn read_parquet(path: &str) -> RecordBatch {
    let file = std::fs::File::open(path).unwrap();
    let reader = ParquetRecordBatchReaderBuilder::try_new(file)
        .unwrap()
        .build()
        .unwrap();
    let batches: Vec<RecordBatch> = reader.map(|b| b.unwrap()).collect();
    let schema = batches[0].schema();
    arrow::compute::concat_batches(&schema, &batches).unwrap()
}

fn strings(batch: &RecordBatch, col: &str) -> Vec<String> {
    let a = batch
        .column_by_name(col)
        .unwrap()
        .as_any()
        .downcast_ref::<StringArray>()
        .unwrap();
    (0..a.len()).map(|i| a.value(i).to_string()).collect()
}

/// One batch, six cells. Cells C1..C4 carry counts; C5 is empty in the count
/// matrix; C6 has a single gene. Sites: S1 strong and seen in 4 cells, S2
/// strong but seen in 1 cell, S3 weak (p = 0.6), S4 thin coverage.
fn build_fixture(dir: &str) {
    let cells = ["C1", "C2", "C3", "C4", "C5", "C6"];
    let genes = [
        "G1/count/spliced",
        "G1/count/unspliced",
        "G2/count/spliced",
        "G3/count/spliced",
    ];
    let mut entries = Vec::new();
    for c in 0..4 {
        for r in 0..4 {
            entries.push((r, c, 5.0 + c as f32));
        }
    }
    entries.push((0, 5, 1.0)); // C6: one gene only
    write_matrix(dir, "b1_count", &genes, &cells, &entries);

    // Site matrix: {gene}/m6a/{chr}:{pos}/{channel}; pos = conversion_pos.
    let site_rows = [
        "G1/m6a/chr1:100/methylated",
        "G1/m6a/chr1:100/unmethylated",
        "G1/m6a/chr1:200/methylated",
        "G1/m6a/chr1:200/unmethylated",
        "G2/m6a/chr1:300/methylated",
        "G2/m6a/chr1:300/unmethylated",
        "G3/m6a/chr1:400/methylated",
        "G3/m6a/chr1:400/unmethylated",
    ];
    let mut site_entries = Vec::new();
    for c in 0..4 {
        site_entries.push((0, c, 2.0)); // S1 methylated in C1..C4
        site_entries.push((1, c, 3.0));
        site_entries.push((5, c, 4.0)); // S3 unmethylated everywhere
    }
    site_entries.push((2, 0, 3.0)); // S2 methylated in C1 only
    site_entries.push((3, 0, 1.0));
    site_entries.push((4, 1, 1.0)); // S3 methylated in C2
    site_entries.push((6, 4, 1.0)); // S4 methylated in the dropped cell C5 only
    site_entries.push((7, 4, 1.0));
    write_matrix(dir, "b1_m6a_site", &site_rows, &cells, &site_entries);

    // Producer's gene-level pool (all sites), which qc must re-derive.
    write_matrix(
        dir,
        "b1_m6a",
        &[
            "G1/m6a/methylated",
            "G1/m6a/unmethylated",
            "G2/m6a/methylated",
            "G2/m6a/unmethylated",
        ],
        &cells,
        &[(0, 0, 99.0), (1, 0, 99.0), (2, 1, 99.0), (3, 1, 99.0)],
    );

    write_sites(
        &format!("{dir}/m6a_sites.parquet"),
        &[
            Site {
                chr: "chr1",
                gene: "G1",
                pos: 100,
                pv: 1e-4,
                log_odds: Some(3.0),
                cov: 20,
                conv: 8,
                ctl_cov: 20,
                ctl_conv: 0,
            },
            Site {
                chr: "chr1",
                gene: "G1",
                pos: 200,
                pv: 1e-3,
                log_odds: Some(2.0),
                cov: 10,
                conv: 3,
                ctl_cov: 10,
                ctl_conv: 0,
            },
            Site {
                chr: "chr1",
                gene: "G2",
                pos: 300,
                pv: 0.6,
                log_odds: Some(0.1),
                cov: 20,
                conv: 1,
                ctl_cov: 20,
                ctl_conv: 1,
            },
            Site {
                chr: "chr1",
                gene: "G3",
                pos: 400,
                pv: 1e-2,
                log_odds: Some(1.0),
                cov: 2,
                conv: 1,
                ctl_cov: 0,
                ctl_conv: 0,
            },
        ],
    );
    std::fs::write(format!("{dir}/pipeline_summary.json"), "{}").unwrap();
}

#[test]
fn qc_filters_cells_sites_and_repools_gene_level() {
    let tmp = tempfile::tempdir().unwrap();
    let input = format!("{}/in", tmp.path().display());
    let output = format!("{}/out", tmp.path().display());
    std::fs::create_dir_all(&input).unwrap();
    build_fixture(&input);

    let status = faba()
        .args([
            "qc",
            &input,
            "--batch-process",
            "-o",
            &output,
            "--no-cell-qc",
            "--column-nnz-cutoff",
            "2",
            "--site-max-pv",
            "0.05",
            "--site-min-coverage",
            "3",
            "--site-min-cells",
            "2",
            "--no-zip",
        ])
        .status()
        .unwrap();
    assert!(status.success(), "faba qc failed");

    // Cells: C1..C4 (nnz 4 >= 2); C5 empty and C6 (nnz 1) drop.
    let count =
        open_sparse_matrix(&format!("{output}/b1_count.zarr"), &SparseIoBackend::Zarr).unwrap();
    let cols = count.column_names().unwrap();
    assert_eq!(
        cols.iter().map(|c| c.as_ref()).collect::<Vec<_>>(),
        ["C1", "C2", "C3", "C4"]
    );

    // Sites: S1 kept; S2 dropped (1 cell < 2); S3 dropped (equal rates in both
    // arms: raw log odds exactly 0, checked before its p of 0.6); S4 dropped
    // (coverage 2 < 3, checked before cells).
    let kept = read_parquet(&format!("{output}/m6a_sites.parquet"));
    assert_eq!(kept.num_rows(), 1);
    assert_eq!(strings(&kept, "gene"), ["G1"]);
    let dropped = read_parquet(&format!("{output}/m6a_sites_dropped.parquet"));
    assert_eq!(dropped.num_rows(), 3);
    let reasons = strings(&dropped, "reason");
    assert_eq!(reasons, ["cells", "log_odds", "coverage"]);

    // Site matrix: both channels of S1 only, over the kept cells.
    let site = open_sparse_matrix(
        &format!("{output}/b1_m6a_site.zarr"),
        &SparseIoBackend::Zarr,
    )
    .unwrap();
    let rows = site.row_names().unwrap();
    assert_eq!(
        rows.iter().map(|r| r.as_ref()).collect::<Vec<_>>(),
        ["G1/m6a/chr1:100/methylated", "G1/m6a/chr1:100/unmethylated"]
    );
    assert_eq!(site.num_columns(), Some(4));

    // Gene level re-pooled from the kept site rows, not the producer's 99s:
    // G1 methylated = 2 per cell, unmethylated = 3 per cell; G2 is gone.
    let gene =
        open_sparse_matrix(&format!("{output}/b1_m6a.zarr"), &SparseIoBackend::Zarr).unwrap();
    let rows = gene.row_names().unwrap();
    assert_eq!(
        rows.iter().map(|r| r.as_ref()).collect::<Vec<_>>(),
        ["G1/m6a/methylated", "G1/m6a/unmethylated"]
    );
    let (_, _, triplets) = gene.read_triplets_by_columns((0..4).collect()).unwrap();
    let mut sums = [0.0f32; 2];
    for (r, _, x) in triplets {
        sums[r as usize] += x;
    }
    assert_eq!(sums, [8.0, 12.0]);

    // The input's run records describe the input: qc writes its own instead.
    assert!(!std::path::Path::new(&format!("{output}/pipeline_summary.json")).exists());
    let record: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(format!("{output}/qc.run.json")).unwrap())
            .unwrap();
    assert_eq!(record["status"], "ok");
    assert!(record["inputs"]["fileset"]["path"].is_string());
    assert!(std::path::Path::new(&format!("{output}/b1_cells.tsv.gz")).exists());
    assert!(std::path::Path::new(&format!("{output}/qc_summary.tsv")).exists());
}

const CUT: [&str; 9] = [
    "--no-cell-qc",
    "--column-nnz-cutoff",
    "2",
    "--site-max-pv",
    "0.05",
    "--site-min-coverage",
    "3",
    "--site-min-cells",
    "2",
];

/// `--dry-run` prints what the cut above keeps and removes, and writes
/// nothing.
#[test]
fn qc_dry_run_counts_without_writing() {
    let tmp = tempfile::tempdir().unwrap();
    let input = format!("{}/in", tmp.path().display());
    std::fs::create_dir_all(&input).unwrap();
    build_fixture(&input);
    let before = std::fs::read_dir(tmp.path()).unwrap().count();
    let out = faba()
        .args(["qc", &input, "--dry-run"])
        .args(CUT)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(std::fs::read_dir(tmp.path()).unwrap().count(), before);
    let text = String::from_utf8(out.stdout).unwrap();
    let rows: Vec<Vec<&str>> = text.lines().map(|l| l.split('\t').collect()).collect();
    assert_eq!(
        rows[0],
        [
            "level",
            "name",
            "rule",
            "threshold",
            "removes",
            "only_this",
            "kept",
            "of"
        ]
    );
    // C5 and C6 drop, as in the cut itself.
    let cells = rows.iter().find(|r| r[0] == "cells").expect("a cells row");
    assert_eq!(
        (cells[1], cells[4], cells[6], cells[7]),
        ("b1", "2", "4", "6")
    );
    // One site of four is kept, whichever threshold the row is for.
    let sites: Vec<_> = rows.iter().filter(|r| r[0] == "sites").collect();
    assert!(!sites.is_empty());
    for r in &sites {
        assert_eq!((r[6], r[7]), ("1", "4"), "{r:?}");
    }
    let cells_rule = sites
        .iter()
        .find(|r| r[2] == "--site-min-cells")
        .expect("a --site-min-cells row");
    assert_eq!(cells_rule[5], "1", "only S2 falls to the cells rule alone");
}

/// `--plots` draws the view's figures in batch: one PDF and PNG per
/// modality and knob.
#[test]
fn qc_plots_are_drawn_in_batch() {
    let tmp = tempfile::tempdir().unwrap();
    let input = format!("{}/in", tmp.path().display());
    let output = format!("{}/out", tmp.path().display());
    std::fs::create_dir_all(&input).unwrap();
    build_fixture(&input);
    let status = faba()
        .args([
            "qc",
            &input,
            "--batch-process",
            "--plots",
            "--no-zip",
            "-o",
            &output,
        ])
        .args(CUT)
        .status()
        .unwrap();
    assert!(status.success());
    let mut plots: Vec<String> = std::fs::read_dir(format!("{output}/qc_plots"))
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    plots.sort();
    assert!(plots.iter().any(|p| p.ends_with(".png")), "{plots:?}");
    for png in plots.iter().filter(|p| p.ends_with(".png")) {
        let pdf = png.replace(".png", ".pdf");
        assert!(plots.contains(&pdf), "{png} without {pdf}");
    }
}

/// `--plots --panel-only` draws each figure's plot alone at the size asked,
/// however small; without --panel-only a page that narrow is refused,
/// naming it.
#[test]
fn qc_panels_are_drawn_alone_at_any_size() {
    let tmp = tempfile::tempdir().unwrap();
    let input = format!("{}/in", tmp.path().display());
    std::fs::create_dir_all(&input).unwrap();
    build_fixture(&input);
    let run = |output: &str, alone: &[&str]| {
        faba()
            .args(["qc", &input, "--batch-process", "--plots", "--no-zip"])
            .args(["-o", output, "--fig-width", "2.5", "--fig-height", "1.8"])
            .args(["--font-size", "6", "--palette", "gray"])
            .args(alone)
            .args(CUT)
            .output()
            .unwrap()
    };
    let refused = run(&format!("{}/narrow", tmp.path().display()), &[]);
    assert!(!refused.status.success());
    let said = String::from_utf8_lossy(&refused.stderr);
    assert!(said.contains("--panel-only"), "{said}");
    let output = format!("{}/out", tmp.path().display());
    assert!(run(&output, &["--panel-only"]).status.success());
    let pdfs: Vec<std::path::PathBuf> = std::fs::read_dir(format!("{output}/qc_plots"))
        .unwrap()
        .map(|e| e.unwrap().path())
        .filter(|p| p.extension().is_some_and(|x| x == "pdf"))
        .collect();
    assert!(!pdfs.is_empty());
    // The pairs alone are the site panels, a figure each: no 2x2 grid.
    let pairs = format!("{output}/qc_plots/m6a_pairs.pdf");
    assert!(!std::path::Path::new(&pairs).exists(), "{pairs}");
    for path in &pdfs {
        assert_eq!(pdf_size(path), (180.0, 129.6), "{}", path.display());
    }
}

/// The page of the PDF at `path`, width by height in points, to a tenth.
fn pdf_size(path: &std::path::Path) -> (f64, f64) {
    let pdf = std::fs::read(path).unwrap();
    let key = b"/MediaBox";
    let at = (pdf.windows(key.len()))
        .position(|w| w == key)
        .unwrap_or_else(|| panic!("{}: no page size", path.display()));
    let rest = String::from_utf8_lossy(&pdf[at + key.len()..at + key.len() + 80]).into_owned();
    let inner = rest.split('[').nth(1).and_then(|r| r.split(']').next());
    let v: Vec<f64> = (inner.unwrap_or_default().split_whitespace())
        .map(|x| x.parse().unwrap())
        .collect();
    let tenth = |x: f64| (x * 10.0).round() / 10.0;
    (tenth(v[2] - v[0]), tenth(v[3] - v[1]))
}

/// A knob's figure is drawn whatever scale `--scale` gives its quantity
/// (which bins it reads is checked on the drawing itself, in the unit
/// tests).
#[test]
fn qc_knob_figures_read_the_scale_asked() {
    let tmp = tempfile::tempdir().unwrap();
    let input = format!("{}/in", tmp.path().display());
    std::fs::create_dir_all(&input).unwrap();
    build_fixture(&input);
    for (k, scale) in [None, Some("pv=linear"), Some("pv=log")]
        .into_iter()
        .enumerate()
    {
        let output = format!("{}/out{k}", tmp.path().display());
        let mut run = faba();
        run.args(["qc", &input, "--dry-run", "--plots", "-o", &output]);
        if let Some(s) = scale {
            run.args(["--scale", s]);
        }
        let out = run.args(CUT).output().unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let pdf = format!("{output}/qc_plots/m6a_max-pv.pdf");
        assert!(std::path::Path::new(&pdf).exists(), "{scale:?}");
    }
}

/// The files directly in `dir`, sorted.
fn listed(dir: &str) -> Vec<String> {
    let mut v: Vec<String> = std::fs::read_dir(dir)
        .unwrap()
        .map(|e| e.unwrap().file_name().to_string_lossy().into_owned())
        .collect();
    v.sort();
    v
}

/// `faba site-plot` writes exactly the figures asked for (PDF and PNG
/// each, and its run record), prints their paths, and leaves the input as
/// it was; it never writes into a directory holding files.
#[test]
fn site_plot_writes_exactly_the_figures_asked() {
    let tmp = tempfile::tempdir().unwrap();
    let input = format!("{}/in", tmp.path().display());
    let output = format!("{}/out", tmp.path().display());
    std::fs::create_dir_all(&input).unwrap();
    build_fixture(&input);
    let before = listed(&input);
    let out = faba()
        .args(["site-plot", &input, "-o", &output])
        .args([
            "--hist",
            "pv,cells",
            "--pairs",
            "pv:cells",
            "--marginals",
            "top,right",
        ])
        .args(["--site-max-pv", "0.05"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let mut want: Vec<String> = ["m6a_hist_pv", "m6a_hist_cells", "m6a_cells_vs_pv"]
        .iter()
        .flat_map(|s| ["pdf", "png"].map(|x| format!("{s}.{x}")))
        .collect();
    let printed: Vec<String> = String::from_utf8(out.stdout)
        .unwrap()
        .lines()
        .map(|l| l.rsplit('/').next().unwrap().to_string())
        .collect();
    let mut sorted = printed.clone();
    sorted.sort();
    want.sort();
    assert_eq!(sorted, want, "printed");
    want.push("site-plot.run.json".into());
    want.sort();
    assert_eq!(listed(&output), want);
    assert_eq!(listed(&input), before);
    let again = faba()
        .args(["site-plot", &input, "-o", &output, "--hist", "pv"])
        .output()
        .unwrap();
    assert!(!again.status.success());
    assert_eq!(listed(&output), want, "nothing written over");
    let bad = faba()
        .args([
            "site-plot",
            &input,
            "-o",
            &format!("{output}2"),
            "--pairs",
            "pv:nope",
        ])
        .output()
        .unwrap();
    assert!(!bad.status.success());
    assert!(String::from_utf8_lossy(&bad.stderr).contains("nope"));
}

/// A site-plot panel is the size asked, `--width` and `--height` too, with
/// `--panel-only` and `--font-size` (each label's size is checked on the
/// drawing itself, in the unit tests).
#[test]
fn site_plot_panels_are_the_size_and_font_asked() {
    let tmp = tempfile::tempdir().unwrap();
    let input = format!("{}/in", tmp.path().display());
    let output = format!("{}/out", tmp.path().display());
    std::fs::create_dir_all(&input).unwrap();
    build_fixture(&input);
    let out = faba()
        .args(["site-plot", &input, "-o", &output, "--panel-only"])
        .args(["--width", "1.8", "--fig-height", "1.4", "--font-size", "6"])
        .args([
            "--pairs",
            "coverage:cells",
            "--marginals",
            "right",
            "--palette",
            "gray",
        ])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let pdfs: Vec<String> = listed(&output)
        .into_iter()
        .filter(|f| f.ends_with(".pdf"))
        .collect();
    assert!(
        pdfs.contains(&"m6a_cells_vs_coverage.pdf".to_string()),
        "{pdfs:?}"
    );
    assert!(
        pdfs.contains(&"m6a_hist_fold.pdf".to_string()),
        "every histogram: {pdfs:?}"
    );
    for name in &pdfs {
        let path = std::path::Path::new(&output).join(name);
        assert_eq!(pdf_size(&path), (129.6, 100.8), "{name}");
    }
}

/// `faba qc --plots` writes the same panels at its thresholds under
/// `qc_plots/panels/`, in a dry run too.
#[test]
fn qc_plots_write_the_site_panels() {
    let tmp = tempfile::tempdir().unwrap();
    let input = format!("{}/in", tmp.path().display());
    std::fs::create_dir_all(&input).unwrap();
    build_fixture(&input);
    let panels = [
        "--hist",
        "pv,coverage",
        "--pairs",
        "pv:cells",
        "--marginals",
        "top",
        "--scale",
        "pv=log",
    ];
    for dry in [false, true] {
        let output = format!("{}/out{dry}", tmp.path().display());
        let mode = if dry { "--dry-run" } else { "--batch-process" };
        let out = faba()
            .args(["qc", &input, mode, "--plots", "--no-zip", "-o", &output])
            .args(panels)
            .args(CUT)
            .output()
            .unwrap();
        assert!(
            out.status.success(),
            "{}",
            String::from_utf8_lossy(&out.stderr)
        );
        let got: Vec<String> = listed(&format!("{output}/qc_plots/panels"))
            .into_iter()
            .filter(|f| f.ends_with(".pdf"))
            .collect();
        assert_eq!(
            got,
            [
                "m6a_cells_vs_pv.pdf",
                "m6a_hist_coverage.pdf",
                "m6a_hist_pv.pdf"
            ],
            "dry run {dry}"
        );
    }
    // The panel flags need --plots.
    let refused = faba()
        .args(["qc", &input, "--dry-run", "--hist", "pv"])
        .output()
        .unwrap();
    assert!(!refused.status.success());
}

/// The `sites` rows of a printed counts table, by flag: `(removes,
/// only_this, kept, of)`.
fn site_counts(text: &str) -> Vec<(String, [String; 4])> {
    text.lines()
        .map(|l| l.split('\t').collect::<Vec<_>>())
        .filter(|r| r[0] == "sites")
        .map(|r| (r[2].to_string(), [4, 5, 6, 7].map(|i| r[i].to_string())))
        .collect()
}

/// `faba qc-view --batch-process` on what `faba qc` wrote prints the same
/// counts the dry run of that cut does: the thresholds recorded, the
/// dropped sites read back with the kept, each one's cells counted as qc
/// counted them (off the input), and each batch's cells off its report.
/// At one cell a site, S3 (dropped for its log odds) passes the cells
/// check only when its cell is counted. With the input gone, it still
/// opens, the dropped sites' cells unknown.
#[test]
fn qc_view_reads_back_the_cut_qc_wrote() {
    let tmp = tempfile::tempdir().unwrap();
    let input = format!("{}/in", tmp.path().display());
    let output = format!("{}/out", tmp.path().display());
    std::fs::create_dir_all(&input).unwrap();
    build_fixture(&input);
    let mut at_one = CUT;
    assert_eq!(at_one[7], "--site-min-cells");
    at_one[8] = "1";
    let cut = faba()
        .args(["qc", &input, "--batch-process", "--no-zip", "-o", &output])
        .args(at_one)
        .output()
        .unwrap();
    assert!(
        cut.status.success(),
        "{}",
        String::from_utf8_lossy(&cut.stderr)
    );
    let dry = faba()
        .args(["qc", &input, "--dry-run"])
        .args(at_one)
        .output()
        .unwrap();
    let viewed = faba()
        .args(["-v", "qc-view", &output, "--batch-process"])
        .output()
        .unwrap();
    let said = String::from_utf8_lossy(&viewed.stderr);
    assert!(viewed.status.success(), "{said}");
    assert!(!said.contains("not on this machine"), "{said}");
    let (dry, viewed) = (
        String::from_utf8(dry.stdout).unwrap(),
        String::from_utf8(viewed.stdout).unwrap(),
    );
    let rows = site_counts(&viewed);
    assert_eq!(rows, site_counts(&dry), "{viewed}");
    let by_cells = rows.iter().find(|r| r.0 == "--site-min-cells").unwrap();
    assert_eq!(by_cells.1[0], "1", "only S4, with no kept cell: {viewed}");
    let cells = |t: &str| {
        t.lines()
            .find(|l| l.starts_with("cells\t"))
            .map(String::from)
    };
    assert_eq!(cells(&viewed), cells(&dry), "{viewed}");

    // Its input moved away: the kept sites' cells off the output's own
    // `_site` matrix, the dropped ones' unknown, as it warns.
    std::fs::rename(&input, format!("{input}.moved")).unwrap();
    let alone = faba()
        .args(["-v", "qc-view", &output, "--batch-process"])
        .output()
        .unwrap();
    assert!(alone.status.success());
    let warned = String::from_utf8_lossy(&alone.stderr);
    assert!(warned.contains("not on this machine"), "{warned}");
    let rows = site_counts(&String::from_utf8(alone.stdout).unwrap());
    assert!(
        rows.iter().all(|(_, r)| r[2] == "2" && r[3] == "4"),
        "{rows:?}"
    );

    // The page, never over a file.
    let page = format!("{}/qc.html", tmp.path().display());
    std::fs::write(&page, "mine").unwrap();
    let refused = faba()
        .args(["qc-view", &output, "--html", &page])
        .output()
        .unwrap();
    assert!(!refused.status.success());
    assert_eq!(std::fs::read_to_string(&page).unwrap(), "mine");
}
