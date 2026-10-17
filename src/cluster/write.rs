//! The clustering's tables.

use super::fit::{Placement, Tiles};
use super::read::{Cell, Site};
use arrow::array::{
    ArrayRef, BooleanArray, Float32Array, Int32Array, Int64Array, StringArray, UInt64Array,
};
use arrow::datatypes::{Field, Schema};
use arrow::record_batch::RecordBatch;
use parquet::file::properties::WriterProperties;
use std::path::Path;
use std::sync::Arc;

/// Rows per row group, so a reader after one gene decodes little else.
const ROW_GROUP: usize = 65_536;

/// Named columns, each field typed by its array, written as one row group
/// per [`ROW_GROUP`] rows.
fn write(columns: Vec<(&str, ArrayRef)>, path: &Path) -> anyhow::Result<()> {
    let fields: Vec<Field> = columns
        .iter()
        .map(|(n, a)| Field::new(*n, a.data_type().clone(), true))
        .collect();
    let batch = RecordBatch::try_new(
        Arc::new(Schema::new(fields)),
        columns.into_iter().map(|c| c.1).collect(),
    )?;
    let props = WriterProperties::builder()
        .set_max_row_group_row_count(Some(ROW_GROUP))
        .build();
    crate::editing::io::write_record_batch_with(&batch, path, props)
}

/// One row per input cell: where it was placed, none when the filters
/// dropped it.
pub struct CellRow<'a> {
    pub cell: &'a Cell,
    pub placed: Option<&'a Placement>,
}

/// Every input cell: barcode, batch, whether it was fitted and kept, and its
/// cluster (1..), confidence and level offset.
pub fn cells(rows: &[CellRow], path: &Path) -> anyhow::Result<()> {
    let opt = |f: fn(&Placement) -> f64| -> ArrayRef {
        Arc::new(Float32Array::from(
            rows.iter()
                .map(|r| r.placed.map(|p| f(p) as f32))
                .collect::<Vec<_>>(),
        ))
    };
    let strings = |f: fn(&Cell) -> &str| -> ArrayRef {
        Arc::new(StringArray::from(
            rows.iter().map(|r| f(r.cell)).collect::<Vec<_>>(),
        ))
    };
    write(
        vec![
            ("barcode", strings(|c| &c.barcode)),
            ("batch", strings(|c| &c.batch)),
            (
                "fitted",
                Arc::new(BooleanArray::from(
                    rows.iter().map(|r| r.cell.fit).collect::<Vec<_>>(),
                )),
            ),
            (
                "kept",
                Arc::new(BooleanArray::from(
                    rows.iter().map(|r| r.placed.is_some()).collect::<Vec<_>>(),
                )),
            ),
            (
                "cluster",
                Arc::new(Int32Array::from(
                    rows.iter()
                        .map(|r| r.placed.map(|p| p.cluster as i32 + 1))
                        .collect::<Vec<_>>(),
                )),
            ),
            ("confidence", opt(|p| p.probability)),
            ("offset", opt(|p| p.offset)),
        ],
        path,
    )
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
            let n = tiles.total(s, k);
            if n == 0 {
                continue;
            }
            let (mu, sigma) = tiles.posterior(s, k);
            gene.push(sites[s].gene.as_ref());
            chr.push(sites[s].chr.as_ref());
            pos.push(sites[s].pos);
            cluster.push(k as i32 + 1);
            on.push(tiles.on(s, k) as u64);
            total.push(n as u64);
            mean.push(mu as f32);
            sd.push(sigma as f32);
        }
    }
    write(
        vec![
            ("gene", Arc::new(StringArray::from(gene))),
            ("chr", Arc::new(StringArray::from(chr))),
            ("pos", Arc::new(Int64Array::from(pos))),
            ("cluster", Arc::new(Int32Array::from(cluster))),
            ("converted", Arc::new(UInt64Array::from(on))),
            ("total", Arc::new(UInt64Array::from(total))),
            ("mean", Arc::new(Float32Array::from(mean))),
            ("sd", Arc::new(Float32Array::from(sd))),
        ],
        path,
    )
}
