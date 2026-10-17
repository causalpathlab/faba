//! The module tables' columns as the browser reads them: the blocks'
//! shifts, the roles, the cell modules' posterior sizes, the per-site reads
//! and each site's counts.

use super::Count;
use crate::cluster::read::{Site, SiteStats};
use crate::common::column;
use anyhow::Context;
use arrow::array::{
    Array, BooleanArray, Float32Array, Int32Array, Int64Array, StringArray, UInt64Array,
};
use arrow::record_batch::RecordBatch;
use rustc_hash::{FxHashMap, FxHashSet};

/// The columns of a table of posterior-mean reads by site and cell role
/// (`{m}_site_cells`, and `{m}_block_sites` but for its gene role).
pub(super) struct SiteRows<'a> {
    pub(super) cell_role: &'a Int32Array,
    pub(super) gene: &'a StringArray,
    chr: &'a StringArray,
    pos: &'a Int64Array,
    converted: &'a Float32Array,
    total: &'a Float32Array,
}

impl<'a> SiteRows<'a> {
    /// The columns of `t`.
    pub(super) fn of(t: &'a RecordBatch) -> anyhow::Result<Self> {
        Ok(Self {
            cell_role: column(t, "cell_role")?,
            gene: column(t, "gene")?,
            chr: column(t, "chr")?,
            pos: column(t, "pos")?,
            converted: column(t, "converted")?,
            total: column(t, "total")?,
        })
    }

    /// How many rows there are.
    pub(super) fn len(&self) -> usize {
        self.pos.len()
    }

    /// Row `i`'s site: its gene, chromosome and position.
    pub(super) fn key(&self, i: usize) -> (&'a str, &'a str, i64) {
        (self.gene.value(i), self.chr.value(i), self.pos.value(i))
    }

    /// Row `i`'s posterior-mean reads.
    pub(super) fn reads(&self, i: usize) -> Count {
        Count {
            on: self.converted.value(i) as f64,
            total: self.total.value(i) as f64,
        }
    }
}

/// The columns of `{m}_block_sites` the browser reads.
pub(super) struct BlockSites<'a> {
    pub(super) rows: SiteRows<'a>,
    pub(super) gene_role: &'a Int32Array,
}

impl<'a> BlockSites<'a> {
    /// The columns of `t`.
    pub(super) fn of(t: &'a RecordBatch) -> anyhow::Result<Self> {
        Ok(Self {
            rows: SiteRows::of(t)?,
            gene_role: column(t, "gene_role")?,
        })
    }
}

/// Row `i`'s role in `roles` (from 1) as an index, when it is one of `n`.
pub(in crate::cluster::browse) fn role_at(
    roles: &Int32Array,
    i: usize,
    n: usize,
) -> anyhow::Result<usize> {
    let r = if roles.is_valid(i) { roles.value(i) } else { 0 };
    let index = usize::try_from(r).ok().and_then(|r| r.checked_sub(1));
    index
        .filter(|&k| k < n)
        .ok_or_else(|| anyhow::anyhow!("row {i} has role {r}, not one of 1 to {n}"))
}

/// The roles and each block's shift from `{m}_blocks`.
pub(super) fn read_shifts(blocks: &RecordBatch) -> anyhow::Result<(usize, usize, Vec<f64>)> {
    let cell_role: &Int32Array = column(blocks, "cell_role")?;
    let gene_role: &Int32Array = column(blocks, "gene_role")?;
    let shift_col: &Float32Array = column(blocks, "shift")?;
    let roles = |a: &Int32Array| a.iter().flatten().max().unwrap_or(0).max(0) as usize;
    let (cell_roles, gene_roles) = (roles(cell_role), roles(gene_role));
    anyhow::ensure!(cell_roles > 0 && gene_roles > 0, "no blocks in the table");
    let mut shift = vec![0.0; cell_roles * gene_roles];
    for i in 0..blocks.num_rows() {
        let c = role_at(cell_role, i, cell_roles).context("blocks cell_role")?;
        let g = role_at(gene_role, i, gene_roles).context("blocks gene_role")?;
        shift[c * gene_roles + g] = shift_col.value(i) as f64;
    }
    Ok((cell_roles, gene_roles, shift))
}

/// Each row's main role (the `role` column, from 1, one of `n`) by its
/// `key`; rows without a role are left out.
pub(super) fn main_roles<K: Eq + std::hash::Hash>(
    table: &RecordBatch,
    key: impl Fn(usize) -> K,
    n: usize,
) -> anyhow::Result<FxHashMap<K, usize>> {
    let role: &Int32Array = column(table, "role")?;
    (0..table.num_rows())
        .filter(|&i| role.is_valid(i))
        .map(|i| Ok((key(i), role_at(role, i, n)?)))
        .collect()
}

/// How many of `roles` each of `n` roles has.
pub(super) fn role_sizes(roles: impl Iterator<Item = usize>, n: usize) -> Vec<usize> {
    let mut sizes = vec![0; n];
    for r in roles {
        if let Some(s) = sizes.get_mut(r) {
            *s += 1;
        }
    }
    sizes
}

/// Each of `n` cell roles' posterior size: the memberships (`role_1` ..
/// in `{m}_cell_roles`) of the fitted cells, summed; and whether they are
/// approximate. A cell is fitted by the table's `fitted` column; older
/// tables have none, and their fitted cells are taken to be those listed
/// in `{m}_block_cells`, read by `block_cells` (see [`listed_rows`]).
pub(super) fn posterior_sizes(
    cells: &RecordBatch,
    block_cells: impl FnOnce() -> anyhow::Result<RecordBatch>,
    n: usize,
) -> anyhow::Result<(Vec<f64>, bool)> {
    let approximate = cells.column_by_name("fitted").is_none();
    let rows = if approximate {
        listed_rows(cells, &block_cells()?)?
    } else {
        let fitted: &BooleanArray = column(cells, "fitted")?;
        (0..cells.num_rows())
            .filter(|&i| fitted.is_valid(i) && fitted.value(i))
            .collect()
    };
    let sizes = (0..n)
        .map(|c| {
            let weight: &Float32Array = column(cells, &format!("role_{}", c + 1))?;
            Ok(rows
                .iter()
                .filter(|&&i| weight.is_valid(i))
                .map(|&i| weight.value(i) as f64)
                .sum())
        })
        .collect::<anyhow::Result<_>>()?;
    Ok((sizes, approximate))
}

/// The rows of `cells` that `block_cells` lists, by batch and barcode:
/// older tables' fitted cells, less those with under half a link in every
/// block.
fn listed_rows(cells: &RecordBatch, block_cells: &RecordBatch) -> anyhow::Result<Vec<usize>> {
    let (batch, barcode): (&StringArray, &StringArray) = (
        column(block_cells, "batch")?,
        column(block_cells, "barcode")?,
    );
    let listed: FxHashSet<(&str, &str)> = (0..block_cells.num_rows())
        .map(|i| (batch.value(i), barcode.value(i)))
        .collect();
    let (batch, barcode): (&StringArray, &StringArray) =
        (column(cells, "batch")?, column(cells, "barcode")?);
    Ok((0..cells.num_rows())
        .filter(|&i| listed.contains(&(batch.value(i), barcode.value(i))))
        .collect())
}

/// Each site's counts from a `{m}_site_stats` table; a site it lacks
/// counts nothing.
pub(super) fn site_counts(sites: &[Site], t: &RecordBatch) -> anyhow::Result<SiteStats> {
    let (gene, chr): (&StringArray, &StringArray) = (column(t, "gene")?, column(t, "chr")?);
    let pos: &Int64Array = column(t, "pos")?;
    let (on, off, cells): (&UInt64Array, &UInt64Array, &UInt64Array) = (
        column(t, "converted")?,
        column(t, "unconverted")?,
        column(t, "cells")?,
    );
    let index: FxHashMap<(&str, &str, i64), usize> = sites
        .iter()
        .enumerate()
        .map(|(i, s)| ((s.gene.as_ref(), s.chr.as_ref(), s.pos), i))
        .collect();
    let mut out = SiteStats::zeros(sites.len());
    for i in 0..t.num_rows() {
        if let Some(&s) = index.get(&(gene.value(i), chr.value(i), pos.value(i))) {
            out.on[s] = on.value(i);
            out.off[s] = off.value(i);
            out.cells[s] = cells.value(i);
        }
    }
    Ok(out)
}
