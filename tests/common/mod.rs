//! Fixtures the session tests share: a faba output directory for the qc
//! views.

use std::path::Path;
use std::sync::Arc;

use arrow::array::{ArrayRef, Float32Array, Int64Array, StringArray, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use data_beans::sparse_io::{create_sparse_from_triplets_owned, SparseIoBackend};
use parquet::arrow::ArrowWriter;

/// `{dir}/{name}.zarr` from named triplets.
fn write_matrix(dir: &Path, name: &str, rows: &[&str], cols: &[&str], entries: &[(u64, u64)]) {
    let path = dir.join(format!("{name}.zarr"));
    let triplets: Vec<(u64, u64, f32)> = entries.iter().map(|&(r, c)| (r, c, 3.0)).collect();
    let mut data = create_sparse_from_triplets_owned(
        triplets,
        (rows.len(), cols.len(), entries.len()),
        Some(&path.to_string_lossy()),
        Some(&SparseIoBackend::Zarr),
    )
    .unwrap();
    let rows: Vec<Box<str>> = rows.iter().map(|s| (*s).into()).collect();
    let cols: Vec<Box<str>> = cols.iter().map(|s| (*s).into()).collect();
    data.register_row_names_vec(&rows);
    data.register_column_names_vec(&cols);
}

/// A faba output directory in `dir` for `faba qc`: a count matrix, an m6A
/// site matrix and its site table, for one batch of four cells.
pub fn qc_fixture(dir: &Path) {
    let cells = ["C1", "C2", "C3", "C4"];
    let counts: Vec<(u64, u64)> = (0..2).flat_map(|r| (0..4).map(move |c| (r, c))).collect();
    write_matrix(
        dir,
        "b1_count",
        &["G1/count/spliced", "G2/count/spliced"],
        &cells,
        &counts,
    );
    let sites: Vec<(u64, u64)> = (0..4).flat_map(|r| (0..4).map(move |c| (r, c))).collect();
    write_matrix(
        dir,
        "b1_m6a_site",
        &[
            "G1/m6a/chr1:100/methylated",
            "G1/m6a/chr1:100/unmethylated",
            "G2/m6a/chr1:300/methylated",
            "G2/m6a/chr1:300/unmethylated",
        ],
        &cells,
        &sites,
    );
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
            Arc::new(StringArray::from(vec!["chr1", "chr1"])) as ArrayRef,
            Arc::new(StringArray::from(vec!["G1", "G2"])) as ArrayRef,
            Arc::new(Int64Array::from(vec![99, 299])) as ArrayRef,
            Arc::new(Int64Array::from(vec![Some(100), Some(300)])) as ArrayRef,
            Arc::new(Float32Array::from(vec![1e-4, 1e-3])) as ArrayRef,
            Arc::new(Float32Array::from(vec![Some(3.0), Some(2.0)])) as ArrayRef,
            Arc::new(UInt64Array::from(vec![20, 10])) as ArrayRef,
            Arc::new(UInt64Array::from(vec![8, 3])) as ArrayRef,
            Arc::new(UInt64Array::from(vec![20, 10])) as ArrayRef,
            Arc::new(UInt64Array::from(vec![0, 0])) as ArrayRef,
        ],
    )
    .unwrap();
    let file = std::fs::File::create(dir.join("m6a_sites.parquet")).unwrap();
    let mut w = ArrowWriter::try_new(file, schema, None).unwrap();
    w.write(&batch).unwrap();
    w.close().unwrap();
    std::fs::write(dir.join("pipeline_summary.json"), "{}").unwrap();
}
