use super::*;
use arrow::array::{ArrayRef, Int64Array, StringArray, UInt64Array};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use std::sync::Arc;

/// A site table of `(chr, pos, coverage, converted)` rows; without
/// `with_reads`, no `coverage` or `converted` column.
fn write(path: &std::path::Path, rows: &[(&str, i64, u64, u64)], with_reads: bool) {
    let mut fields = vec![
        Field::new("chr", DataType::Utf8, false),
        Field::new("primary_pos", DataType::Int64, false),
        Field::new("strand", DataType::Utf8, false),
    ];
    let mut cols: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from(
            rows.iter().map(|r| r.0).collect::<Vec<_>>(),
        )),
        Arc::new(Int64Array::from(
            rows.iter().map(|r| r.1).collect::<Vec<_>>(),
        )),
        Arc::new(StringArray::from(vec!["+"; rows.len()])),
    ];
    if with_reads {
        fields.push(Field::new("coverage", DataType::UInt64, false));
        fields.push(Field::new("converted", DataType::UInt64, false));
        cols.push(Arc::new(UInt64Array::from(
            rows.iter().map(|r| r.2).collect::<Vec<_>>(),
        )));
        cols.push(Arc::new(UInt64Array::from(
            rows.iter().map(|r| r.3).collect::<Vec<_>>(),
        )));
    }
    let schema = Arc::new(Schema::new(fields));
    let batch = RecordBatch::try_new(schema.clone(), cols).unwrap();
    let mut w = ArrowWriter::try_new(File::create(path).unwrap(), schema, None).unwrap();
    w.write(&batch).unwrap();
    w.close().unwrap();
}

#[test]
fn a_sites_reads_are_its_rows_summed_or_none_without_the_columns() {
    let tmp = tempfile::tempdir().unwrap();
    let path = tmp.path().join("sites.parquet");
    let rows = [("chr1", 10, 8, 3), ("chr1", 10, 2, 1), ("chr2", 5, 4, 4)];
    write(&path, &rows, true);
    let file = path.to_string_lossy();
    let sites = read_sites(&file).unwrap();
    let mut asked = read_sites(&file).unwrap();
    asked.push(GenomicSite {
        chr: "chr3".into(),
        position: 1,
        strand: Strand::Forward,
    });
    let reads = read_site_reads(&file, &asked)
        .unwrap()
        .expect("columns there");
    let by_site: Vec<(String, i64, (f64, f64))> = asked
        .iter()
        .zip(&reads)
        .map(|(s, &r)| (s.chr.to_string(), s.position, r))
        .collect();
    assert!(
        by_site.contains(&("chr1".into(), 10, (4.0, 6.0))),
        "{by_site:?}"
    );
    assert!(
        by_site.contains(&("chr2".into(), 5, (4.0, 0.0))),
        "{by_site:?}"
    );
    assert!(
        by_site.contains(&("chr3".into(), 1, (0.0, 0.0))),
        "no row, no reads"
    );

    write(&path, &rows, false);
    assert!(read_site_reads(&file, &sites).unwrap().is_none());
}

#[test]
fn apa_sites_keep_their_utr_strand() {
    use crate::apa::cell_assign::ApaSiteAnnotation;
    let site = |id: &str, alpha: i64, strand: Strand| ApaSiteAnnotation {
        site_id: id.into(),
        gene_name: "G".into(),
        chr: "chr1".into(),
        strand,
        genomic_alpha: alpha,
        beta: 30.0,
        genomic_start: alpha - 50,
        genomic_stop: alpha + 50,
        pi_weight: 0.5,
        expected_tail_length: 100.0,
        utr_length: 1000,
    };
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("apa_components.parquet");
    let path = path.to_str().unwrap();
    let written = [
        site("a", 100, Strand::Forward),
        site("b", 200, Strand::Backward),
    ];
    crate::apa::pipeline::write_apa_annotations(&written, path).unwrap();
    let read = read_sites(path).unwrap();
    let got: Vec<(i64, bool)> = read
        .iter()
        .map(|s| (s.position, matches!(s.strand, Strand::Forward)))
        .collect();
    assert_eq!(got, [(100, true), (200, false)]);
}
