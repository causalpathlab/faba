//! The clustering's tables.

use super::fit::Tiles;
use super::read::{Cell, Site};
use arrow::array::{
    ArrayRef, BooleanArray, Float32Array, Int32Array, Int64Array, StringArray, UInt64Array,
};
use arrow::datatypes::{DataType, Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::arrow::ArrowWriter;
use parquet::file::properties::WriterProperties;
use std::path::Path;
use std::sync::Arc;

/// Rows per row group, so a reader after one gene decodes little else.
const ROW_GROUP: usize = 65_536;

fn write(batch: &RecordBatch, path: &Path) -> anyhow::Result<()> {
    let props = WriterProperties::builder()
        .set_max_row_group_row_count(Some(ROW_GROUP))
        .build();
    let mut writer =
        ArrowWriter::try_new(std::fs::File::create(path)?, batch.schema(), Some(props))?;
    writer.write(batch)?;
    writer.close()?;
    Ok(())
}

fn batch(fields: Vec<(&str, DataType, bool, ArrayRef)>) -> anyhow::Result<RecordBatch> {
    let schema = Schema::new(
        fields
            .iter()
            .map(|(n, t, null, _)| Field::new(*n, t.clone(), *null))
            .collect::<Vec<_>>(),
    );
    let columns = fields.into_iter().map(|f| f.3).collect();
    Ok(RecordBatch::try_new(Arc::new(schema), columns)?)
}

/// One row per input cell: its cluster (1.., none when filtered out), whether
/// it was fitted or assigned, its confidence and its level offset.
pub struct CellRow<'a> {
    pub cell: &'a Cell,
    pub kept: bool,
    pub cluster: Option<u16>,
    pub confidence: Option<f64>,
    pub offset: Option<f64>,
}

pub fn cells(rows: &[CellRow], path: &Path) -> anyhow::Result<()> {
    fn str_col<'a>(rows: &'a [CellRow], f: impl Fn(&'a CellRow) -> &'a str) -> ArrayRef {
        Arc::new(StringArray::from(rows.iter().map(f).collect::<Vec<_>>()))
    }
    let f32_col = |f: &dyn Fn(&CellRow) -> Option<f64>| -> ArrayRef {
        Arc::new(Float32Array::from(
            rows.iter()
                .map(|r| f(r).map(|v| v as f32))
                .collect::<Vec<_>>(),
        ))
    };
    let b = batch(vec![
        (
            "barcode",
            DataType::Utf8,
            false,
            str_col(rows, |r| &r.cell.barcode),
        ),
        (
            "batch",
            DataType::Utf8,
            false,
            str_col(rows, |r| &r.cell.batch),
        ),
        (
            "fitted",
            DataType::Boolean,
            false,
            Arc::new(BooleanArray::from(
                rows.iter().map(|r| r.cell.fit).collect::<Vec<_>>(),
            )),
        ),
        (
            "kept",
            DataType::Boolean,
            false,
            Arc::new(BooleanArray::from(
                rows.iter().map(|r| r.kept).collect::<Vec<_>>(),
            )),
        ),
        (
            "cluster",
            DataType::Int32,
            true,
            Arc::new(Int32Array::from(
                rows.iter()
                    .map(|r| r.cluster.map(|k| k as i32 + 1))
                    .collect::<Vec<_>>(),
            )),
        ),
        (
            "confidence",
            DataType::Float32,
            true,
            f32_col(&|r| r.confidence),
        ),
        ("offset", DataType::Float32, true, f32_col(&|r| r.offset)),
    ])?;
    write(&b, path)
}

/// One row per site and cluster with reads, sorted by gene and position:
/// the fitted cells' converted and total reads, and the posterior mean and
/// standard deviation of the converted fraction.
pub fn sites(sites: &[Site], tiles: &Tiles, path: &Path) -> anyhow::Result<()> {
    let mut order: Vec<usize> = (0..sites.len()).collect();
    order.sort_by(|&a, &b| {
        let (x, y) = (&sites[a], &sites[b]);
        (&x.gene, &x.chr, x.pos).cmp(&(&y.gene, &y.chr, y.pos))
    });
    let (mut gene, mut chr, mut pos, mut cluster) = (vec![], vec![], vec![], vec![]);
    let (mut on, mut total, mut mean, mut sd) = (vec![], vec![], vec![], vec![]);
    for s in order {
        for k in 0..tiles.k {
            let (m, n) = (tiles.on(s, k), tiles.total(s, k));
            if n == 0 {
                continue;
            }
            let (mu, sigma) = tiles.posterior(s, k);
            gene.push(sites[s].gene.as_ref());
            chr.push(sites[s].chr.as_ref());
            pos.push(sites[s].pos);
            cluster.push(k as i32 + 1);
            on.push(m as u64);
            total.push(n as u64);
            mean.push(mu as f32);
            sd.push(sigma as f32);
        }
    }
    let b = batch(vec![
        (
            "gene",
            DataType::Utf8,
            false,
            Arc::new(StringArray::from(gene)),
        ),
        (
            "chr",
            DataType::Utf8,
            false,
            Arc::new(StringArray::from(chr)),
        ),
        (
            "pos",
            DataType::Int64,
            false,
            Arc::new(Int64Array::from(pos)),
        ),
        (
            "cluster",
            DataType::Int32,
            false,
            Arc::new(Int32Array::from(cluster)),
        ),
        (
            "converted",
            DataType::UInt64,
            false,
            Arc::new(UInt64Array::from(on)),
        ),
        (
            "total",
            DataType::UInt64,
            false,
            Arc::new(UInt64Array::from(total)),
        ),
        (
            "mean",
            DataType::Float32,
            false,
            Arc::new(Float32Array::from(mean)),
        ),
        (
            "sd",
            DataType::Float32,
            false,
            Arc::new(Float32Array::from(sd)),
        ),
    ])?;
    write(&b, path)
}
