//! The clustering's tables.

use super::fit::{Placement, Tiles};
use super::link::Profile;
use super::read::{Cell, Site, SiteStats};
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

/// Site indices by gene, chromosome and position.
fn site_order(sites: &[Site]) -> Vec<usize> {
    let mut order: Vec<usize> = (0..sites.len()).collect();
    order.sort_by(|&a, &b| site_key(&sites[a]).cmp(&site_key(&sites[b])));
    order
}

/// The key [`site_order`] sorts sites by.
fn site_key(s: &Site) -> (&str, &str, i64) {
    (&s.gene, &s.chr, s.pos)
}

/// The `gene`, `chr` and `pos` columns of the sites at `index`, one row each.
fn site_columns(
    sites: &[Site],
    index: impl Iterator<Item = usize> + Clone,
) -> [(&'static str, ArrayRef); 3] {
    [
        (
            "gene",
            Arc::new(StringArray::from(
                index
                    .clone()
                    .map(|s| sites[s].gene.as_ref())
                    .collect::<Vec<_>>(),
            )),
        ),
        (
            "chr",
            Arc::new(StringArray::from(
                index
                    .clone()
                    .map(|s| sites[s].chr.as_ref())
                    .collect::<Vec<_>>(),
            )),
        ),
        (
            "pos",
            Arc::new(Int64Array::from(
                index.map(|s| sites[s].pos).collect::<Vec<_>>(),
            )),
        ),
    ]
}

/// One row per site and cluster with reads, sorted by gene and position:
/// the fitted cells' converted and total reads, and the posterior mean and
/// standard deviation of the converted fraction.
pub fn sites(sites: &[Site], tiles: &Tiles, path: &Path) -> anyhow::Result<()> {
    let order = site_order(sites);
    let (mut site, mut cluster) = (vec![], vec![]);
    let (mut on, mut total, mut mean, mut sd) = (vec![], vec![], vec![], vec![]);
    for s in order {
        for k in 0..tiles.k {
            let n = tiles.total(s, k);
            if n == 0 {
                continue;
            }
            let (mu, sigma) = tiles.posterior(s, k);
            site.push(s);
            cluster.push(k as i32 + 1);
            on.push(tiles.on(s, k) as u64);
            total.push(n as u64);
            mean.push(mu as f32);
            sd.push(sigma as f32);
        }
    }
    let mut columns = Vec::from(site_columns(sites, site.into_iter()));
    columns.extend([
        ("cluster", Arc::new(Int32Array::from(cluster)) as ArrayRef),
        ("converted", Arc::new(UInt64Array::from(on))),
        ("total", Arc::new(UInt64Array::from(total))),
        ("mean", Arc::new(Float32Array::from(mean))),
        ("sd", Arc::new(Float32Array::from(sd))),
    ]);
    write(columns, path)
}

/// One row per input cell under the link or block model: its membership
/// over the communities or roles and, under the block model, its own
/// methylation level; none when it was not kept.
pub struct CommunityRow<'a> {
    pub cell: &'a Cell,
    pub membership: Option<&'a [f32]>,
    pub level: Option<f32>,
}

/// Every input cell: barcode, batch, its main `name` (a community or role,
/// 1..), that one's weight, its level when any cell has one, and its weight
/// in each.
pub fn communities(rows: &[CommunityRow], k: usize, name: &str, path: &Path) -> anyhow::Result<()> {
    let labels = [
        (
            "barcode",
            rows.iter().map(|r| r.cell.barcode.as_ref()).collect(),
        ),
        (
            "batch",
            rows.iter().map(|r| r.cell.batch.as_ref()).collect(),
        ),
    ];
    let level = rows
        .iter()
        .any(|r| r.level.is_some())
        .then(|| rows.iter().map(|r| r.level).collect());
    let membership: Vec<Option<&[f32]>> = rows.iter().map(|r| r.membership).collect();
    memberships(&labels, &membership, level, k, name, path)
}

/// One row per gene: its main gene role (1..), that role's weight, and its
/// weight in each role (`weights`, genes x `roles`, row-major).
pub fn gene_roles(
    names: &[&str],
    weights: &[f32],
    roles: usize,
    path: &Path,
) -> anyhow::Result<()> {
    let membership: Vec<Option<&[f32]>> = weights.chunks(roles).map(Some).collect();
    memberships(
        &[("gene", names.to_vec())],
        &membership,
        None,
        roles,
        "role",
        path,
    )
}

/// One row per member: its `labels` columns, its main `name` (1..), that
/// one's weight, its `level` when given, and its weight in each of `k`.
fn memberships(
    labels: &[(&str, Vec<&str>)],
    membership: &[Option<&[f32]>],
    level: Option<Vec<Option<f32>>>,
    k: usize,
    name: &str,
    path: &Path,
) -> anyhow::Result<()> {
    let main = |m: &[f32]| super::summary::main_labels(m, k)[0] as usize;
    let mut columns: Vec<(String, ArrayRef)> = labels
        .iter()
        .map(|(n, v)| {
            (
                n.to_string(),
                Arc::new(StringArray::from(v.clone())) as ArrayRef,
            )
        })
        .collect();
    columns.push((
        name.into(),
        Arc::new(Int32Array::from(
            membership
                .iter()
                .map(|m| m.map(|m| main(m) as i32 + 1))
                .collect::<Vec<_>>(),
        )),
    ));
    columns.push((
        "weight".into(),
        Arc::new(Float32Array::from(
            membership
                .iter()
                .map(|m| m.map(|m| m[main(m)]))
                .collect::<Vec<_>>(),
        )),
    ));
    if let Some(level) = level {
        columns.push(("level".into(), Arc::new(Float32Array::from(level))));
    }
    for j in 0..k {
        columns.push((
            format!("{name}_{}", j + 1),
            Arc::new(Float32Array::from(
                membership
                    .iter()
                    .map(|m| m.map(|m| m[j]))
                    .collect::<Vec<_>>(),
            )),
        ));
    }
    write(
        columns
            .iter()
            .map(|(n, a)| (n.as_str(), a.clone()))
            .collect(),
        path,
    )
}

/// One row per site and community with reads, sorted by gene and position:
/// the community (1..), its weighted reads and posterior-mean fraction.
pub fn community_sites(sites: &[Site], profile: &Profile, path: &Path) -> anyhow::Result<()> {
    let order = site_order(sites);
    let k = profile.k;
    let (mut site, mut label) = (vec![], vec![]);
    let (mut on, mut total, mut mean) = (vec![], vec![], vec![]);
    for s in order {
        for j in 0..k {
            let at = s * k + j;
            if profile.total[at] <= 0.0 {
                continue;
            }
            site.push(s);
            label.push(j as i32 + 1);
            on.push(profile.on[at] as f32);
            total.push(profile.total[at] as f32);
            mean.push(profile.mean[at] as f32);
        }
    }
    let mut columns = Vec::from(site_columns(sites, site.into_iter()));
    columns.extend([
        ("community", Arc::new(Int32Array::from(label)) as ArrayRef),
        ("converted", Arc::new(Float32Array::from(on))),
        ("total", Arc::new(Float32Array::from(total))),
        ("mean", Arc::new(Float32Array::from(mean))),
    ]);
    write(columns, path)
}

/// One row per gene and community holding any of its links: the gene, the
/// community (1..), its links, and their share of the gene's.
/// Genes whose links fall in the same communities across cells form modules.
pub fn gene_communities(
    gene_of_link: &[u32],
    names: &[&str],
    z: &[u16],
    k: usize,
    path: &Path,
) -> anyhow::Result<()> {
    let mut count = vec![0u32; names.len() * k];
    for (&g, &j) in gene_of_link.iter().zip(z) {
        count[g as usize * k + j as usize] += 1;
    }
    let mut order: Vec<usize> = (0..names.len()).collect();
    order.sort_by_key(|&g| names[g]);
    let (mut gene, mut label, mut links, mut share) = (vec![], vec![], vec![], vec![]);
    for g in order {
        let row = &count[g * k..(g + 1) * k];
        let total: u32 = row.iter().sum();
        for (j, &n) in row.iter().enumerate().filter(|(_, &n)| n > 0) {
            gene.push(names[g]);
            label.push(j as i32 + 1);
            links.push(n);
            share.push(n as f32 / total as f32);
        }
    }
    write(
        vec![
            ("gene", Arc::new(StringArray::from(gene))),
            ("community", Arc::new(Int32Array::from(label))),
            ("links", Arc::new(arrow::array::UInt32Array::from(links))),
            ("share", Arc::new(Float32Array::from(share))),
        ],
        path,
    )
}

/// The (cell role, gene role) of block `k` of `gene_roles` per cell role,
/// both numbered from 1.
fn roles(k: usize, gene_roles: usize) -> (i32, i32) {
    ((k / gene_roles) as i32 + 1, (k % gene_roles) as i32 + 1)
}

/// One row per (cell role, gene role) block: both roles (1..), its
/// posterior-mean methylation shift (logit scale) and sd, and its share of
/// the links.
pub fn blocks(
    shift: &[f64],
    sd: &[f64],
    links: &[f64],
    gene_roles: usize,
    path: &Path,
) -> anyhow::Result<()> {
    let n = shift.len();
    let total: f64 = links.iter().sum::<f64>().max(1e-12);
    let (cell, gene): (Vec<i32>, Vec<i32>) = (0..n).map(|k| roles(k, gene_roles)).unzip();
    let f32s = |x: &[f64]| {
        Arc::new(Float32Array::from(
            x.iter().map(|&v| v as f32).collect::<Vec<_>>(),
        ))
    };
    write(
        vec![
            ("cell_role", Arc::new(Int32Array::from(cell))),
            ("gene_role", Arc::new(Int32Array::from(gene))),
            ("shift", f32s(shift)),
            ("sd", f32s(sd)),
            (
                "link_share",
                f32s(&links.iter().map(|l| l / total).collect::<Vec<_>>()),
            ),
        ],
        path,
    )
}

/// One row per block and member (a gene or a cell) with links there, by
/// block and then enrichment, highest first: both roles, the member's
/// `labels` columns, its posterior-mean links in the block, and its
/// enrichment (its share of links in the block over the block's share).
pub fn block_members(
    labels: &[(&str, Vec<&str>)],
    links: &[f64],
    enrichment: &[f64],
    blocks: usize,
    gene_roles: usize,
    path: &Path,
) -> anyhow::Result<()> {
    let mut rows: Vec<(usize, usize)> = (0..links.len())
        .filter(|&at| links[at] >= 0.5)
        .map(|at| (at % blocks, at / blocks))
        .collect();
    rows.sort_by(|x, y| {
        x.0.cmp(&y.0)
            .then(enrichment[y.1 * blocks + y.0].total_cmp(&enrichment[x.1 * blocks + x.0]))
    });
    let (cell, gene): (Vec<i32>, Vec<i32>) =
        rows.iter().map(|&(k, _)| roles(k, gene_roles)).unzip();
    let mut columns: Vec<(&str, ArrayRef)> = vec![
        ("cell_role", Arc::new(Int32Array::from(cell))),
        ("gene_role", Arc::new(Int32Array::from(gene))),
    ];
    for (name, values) in labels {
        columns.push((
            name,
            Arc::new(StringArray::from(
                rows.iter().map(|&(_, r)| values[r]).collect::<Vec<_>>(),
            )),
        ));
    }
    let at = |&(k, r): &(usize, usize)| r * blocks + k;
    columns.push((
        "links",
        Arc::new(Float32Array::from(
            rows.iter().map(|x| links[at(x)] as f32).collect::<Vec<_>>(),
        )),
    ));
    columns.push((
        "enrichment",
        Arc::new(Float32Array::from(
            rows.iter()
                .map(|x| enrichment[at(x)] as f32)
                .collect::<Vec<_>>(),
        )),
    ));
    write(columns, path)
}

/// One row per block and site with reads there, by block, gene and
/// position: both roles, the site, its posterior-mean reads in the block,
/// their methylated fraction, and its lift over the site's baseline (logit
/// scale).
pub fn block_sites(
    sites: &[Site],
    on: &[f64],
    off: &[f64],
    lift: &[f64],
    blocks: usize,
    gene_roles: usize,
    path: &Path,
) -> anyhow::Result<()> {
    let mut rows: Vec<(usize, usize)> = (0..on.len())
        .filter(|&at| on[at] + off[at] >= 1.0)
        .map(|at| (at % blocks, at / blocks))
        .collect();
    rows.sort_by(|x, y| {
        x.0.cmp(&y.0)
            .then(site_key(&sites[x.1]).cmp(&site_key(&sites[y.1])))
    });
    let at = |&(k, s): &(usize, usize)| s * blocks + k;
    let (cell, gene): (Vec<i32>, Vec<i32>) =
        rows.iter().map(|&(k, _)| roles(k, gene_roles)).unzip();
    let f32s = |f: &dyn Fn(usize) -> f64| {
        Arc::new(Float32Array::from(
            rows.iter().map(|x| f(at(x)) as f32).collect::<Vec<_>>(),
        ))
    };
    let mut columns: Vec<(&str, ArrayRef)> = vec![
        ("cell_role", Arc::new(Int32Array::from(cell))),
        ("gene_role", Arc::new(Int32Array::from(gene))),
    ];
    columns.extend(site_columns(sites, rows.iter().map(|&(_, s)| s)));
    columns.extend([
        ("converted", f32s(&|i| on[i]) as ArrayRef),
        ("total", f32s(&|i| on[i] + off[i])),
        ("fraction", f32s(&|i| on[i] / (on[i] + off[i]))),
        ("lift", f32s(&|i| lift[i])),
    ]);
    write(columns, path)
}

/// One row per site, by gene and position: its raw reads over every cell
/// read, converted and unconverted, and the cells with any there.
pub fn site_stats(sites: &[Site], stats: &SiteStats, path: &Path) -> anyhow::Result<()> {
    let order = site_order(sites);
    let u64s = |x: &[u64]| -> ArrayRef {
        Arc::new(UInt64Array::from(
            order.iter().map(|&s| x[s]).collect::<Vec<_>>(),
        ))
    };
    let mut columns = Vec::from(site_columns(sites, order.iter().copied()));
    columns.extend([
        ("converted", u64s(&stats.on)),
        ("unconverted", u64s(&stats.off)),
        ("cells", u64s(&stats.cells)),
    ]);
    write(columns, path)
}
