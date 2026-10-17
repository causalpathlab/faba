//! What `faba cluster --model blocks` wrote, read for the browser: each
//! gene module's sites with their reads in each cell module, the shifts,
//! the cells' and genes' main roles, and each site's counts.

use super::order::{diagonal, standardized};
use crate::cluster::read::{Site, SiteStats};
use crate::common::{column, read_record_batch};
use crate::site_analysis::miami::genemodel::GeneModel;
use crate::site_analysis::pileup::query_symbol;
use anyhow::Context;
use arrow::array::{Array, Float32Array, Int32Array, Int64Array, StringArray, UInt64Array};
use arrow::record_batch::RecordBatch;
use genomic_data::coordinates::chr_eq;
use rustc_hash::{FxHashMap, FxHashSet};
use std::ops::Range;
use std::path::Path;

/// Converted reads and all reads.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Count {
    pub on: f64,
    pub total: f64,
}

impl Count {
    pub fn add(&mut self, other: Count) {
        self.on += other.on;
        self.total += other.total;
    }

    pub fn off(self) -> f64 {
        self.total - self.on
    }

    /// The methylated fraction; none without reads.
    pub fn fraction(self) -> Option<f64> {
        (self.total > 0.0).then(|| self.on / self.total)
    }

    /// The methylated fraction as a whole %, `-` without reads.
    pub fn pct(self) -> String {
        self.fraction()
            .map_or("-".into(), |f| format!("{:.0}%", 100.0 * f))
    }
}

/// One site of a gene module, with its reads in each cell module, in the
/// panes' column order (see [`Modules::columns`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Member {
    pub site: usize,
    pub reads: Vec<Count>,
}

/// A cell's `(batch, barcode)`.
type CellKey = (Box<str>, Box<str>);

/// What `faba cluster --model blocks` wrote, as the browser reads it.
#[derive(Debug)]
pub struct Modules {
    /// The sites' modality, e.g. `m6a`.
    pub modality: String,
    pub cell_roles: usize,
    pub gene_roles: usize,
    /// Fitted cells whose main role each cell role is.
    pub cell_sizes: Vec<usize>,
    /// Genes whose main role each gene role is.
    pub gene_sizes: Vec<usize>,
    /// Each block's shift, `c * gene_roles + g`.
    shift: Vec<f64>,
    /// The sites of the genes in their main role, by gene and position.
    pub sites: Vec<Site>,
    /// Per gene role, the sites of the genes whose main role it is, with
    /// their reads there: genes with the most reads first, each gene's
    /// sites by position.
    pub members: Vec<Vec<Member>>,
    /// Per gene role, each gene's span of its members.
    genes: Vec<Vec<Range<usize>>>,
    /// Per gene role, its reads in each cell module (in column order): its
    /// members' summed.
    blocks: Vec<Vec<Count>>,
    /// Gene roles with sites, by the range of their shifts across the cell
    /// roles, widest first.
    pub order: Vec<usize>,
    /// The same gene roles in diagonal order (see [`diagonal`]).
    pub diagonal: Vec<usize>,
    /// The cell roles in diagonal order, the order every pane shows them in:
    /// column `k` is cell role `columns[k]`, and every member's and block's
    /// reads are in this order.
    pub columns: Vec<usize>,
    /// Every cell with a main role, placed cells too.
    pub cells: FxHashMap<CellKey, usize>,
    /// Per site, its raw reads over every cell read and its cells, from
    /// `{m}_site_stats`; none when the tables predate it.
    pub stats: Option<SiteStats>,
}

/// The columns of `{m}_block_sites` the browser reads.
struct BlockSites<'a> {
    cell_role: &'a Int32Array,
    gene_role: &'a Int32Array,
    gene: &'a StringArray,
    chr: &'a StringArray,
    pos: &'a Int64Array,
    converted: &'a Float32Array,
    total: &'a Float32Array,
}

impl<'a> BlockSites<'a> {
    fn of(t: &'a RecordBatch) -> anyhow::Result<Self> {
        Ok(Self {
            cell_role: column(t, "cell_role")?,
            gene_role: column(t, "gene_role")?,
            gene: column(t, "gene")?,
            chr: column(t, "chr")?,
            pos: column(t, "pos")?,
            converted: column(t, "converted")?,
            total: column(t, "total")?,
        })
    }

    fn key(&self, i: usize) -> (&'a str, &'a str, i64) {
        (self.gene.value(i), self.chr.value(i), self.pos.value(i))
    }
}

/// The sites' part of the modules, from `{m}_block_sites`.
struct FromSites {
    sites: Vec<Site>,
    members: Vec<Vec<Member>>,
    genes: Vec<Vec<Range<usize>>>,
    blocks: Vec<Vec<Count>>,
}

impl Modules {
    /// The tables of `modality` in `dir`.
    pub fn load(dir: &Path, modality: &str) -> anyhow::Result<Self> {
        let table = |name: &str| {
            let path = dir.join(format!("{modality}_{name}.parquet"));
            read_record_batch(&path).map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))
        };
        let (cell_roles, gene_roles, shift) = read_shifts(&table("blocks")?)?;
        let cells_table = table("cell_roles")?;
        let (batch, barcode): (&StringArray, &StringArray) = (
            column(&cells_table, "batch")?,
            column(&cells_table, "barcode")?,
        );
        let cells = main_roles(&cells_table, |i| {
            (batch.value(i).into(), barcode.value(i).into())
        })?;
        let cell_sizes = fitted_sizes(&cells, &table("block_cells")?, cell_roles)?;
        let genes_table = table("gene_roles")?;
        let names: &StringArray = column(&genes_table, "gene")?;
        let main = main_roles(&genes_table, |i| names.value(i))?;
        let gene_sizes = role_sizes(main.values().copied(), gene_roles);
        let rows = table("block_sites")?;
        let from = from_sites(&BlockSites::of(&rows)?, &main, (cell_roles, gene_roles))?;
        let stats_path = dir.join(format!("{modality}_site_stats.parquet"));
        let stats = if stats_path.is_file() {
            Some(site_counts(&from.sites, &read_record_batch(&stats_path)?)?)
        } else {
            log::warn!(
                "no {}: rerun `faba cluster --model blocks` for each site's counts",
                stats_path.display()
            );
            None
        };
        let mut modules = Self {
            modality: modality.to_string(),
            cell_roles,
            gene_roles,
            cell_sizes,
            gene_sizes,
            shift,
            sites: from.sites,
            members: from.members,
            genes: from.genes,
            blocks: from.blocks,
            order: Vec::new(),
            diagonal: Vec::new(),
            columns: (0..cell_roles).collect(),
            cells,
            stats,
        };
        modules.set_orders();
        modules.put_in_columns();
        Ok(modules)
    }

    /// Every member's and block's reads put in column order, once, so
    /// read `k` is column `k`'s.
    fn put_in_columns(&mut self) {
        let columns = &self.columns;
        let permute = |reads: &mut Vec<Count>| {
            *reads = columns.iter().map(|&c| reads[c]).collect();
        };
        self.members
            .iter_mut()
            .flatten()
            .for_each(|m| permute(&mut m.reads));
        self.blocks.iter_mut().for_each(permute);
    }

    /// The gene roles with sites by shift range, and the diagonal orders.
    fn set_orders(&mut self) {
        let mut order: Vec<usize> = (0..self.gene_roles)
            .filter(|&g| !self.members[g].is_empty())
            .collect();
        order.sort_by(|&a, &b| {
            self.shift_range(b)
                .total_cmp(&self.shift_range(a))
                .then(a.cmp(&b))
        });
        let z: Vec<Vec<Option<f64>>> = order.iter().map(|&g| standardized(self.block(g))).collect();
        let (rows, columns) = diagonal(&z);
        self.diagonal = rows.iter().map(|&r| order[r]).collect();
        if !z.is_empty() {
            self.columns = columns;
        }
        self.order = order;
    }

    /// Each cell role's column: the inverse of [`Modules::columns`].
    pub fn role_columns(&self) -> Vec<usize> {
        let mut column = vec![0; self.columns.len()];
        for (k, &c) in self.columns.iter().enumerate() {
            column[c] = k;
        }
        column
    }

    /// Gene role `g`'s reads in each cell module, in column order.
    pub fn block(&self, g: usize) -> &[Count] {
        &self.blocks[g]
    }

    /// The range of gene role `g`'s shifts across the cell roles.
    pub fn shift_range(&self, g: usize) -> f64 {
        let s = (0..self.cell_roles).map(|c| self.shift[c * self.gene_roles + g]);
        s.clone().fold(f64::MIN, f64::max) - s.fold(f64::MAX, f64::min)
    }

    /// Gene role `g`'s genes in its order, each a span of its members.
    pub fn genes(&self, g: usize) -> &[Range<usize>] {
        &self.genes[g]
    }

    /// The gene of gene role `g`'s member `i`.
    pub fn gene_of(&self, g: usize, i: usize) -> &str {
        &self.sites[self.members[g][i].site].gene
    }
}

/// Row `i`'s role in `roles` (from 1) as an index, when it is one of `n`.
pub(super) fn role_at(roles: &Int32Array, i: usize, n: usize) -> anyhow::Result<usize> {
    let r = if roles.is_valid(i) { roles.value(i) } else { 0 };
    let index = usize::try_from(r).ok().and_then(|r| r.checked_sub(1));
    index
        .filter(|&k| k < n)
        .ok_or_else(|| anyhow::anyhow!("row {i} has role {r}, not one of 1 to {n}"))
}

/// The roles and each block's shift from `{m}_blocks`.
fn read_shifts(blocks: &RecordBatch) -> anyhow::Result<(usize, usize, Vec<f64>)> {
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

/// Each row's main role (the `role` column, from 1) by its `key`; rows
/// without a role are left out.
fn main_roles<K: Eq + std::hash::Hash>(
    table: &RecordBatch,
    key: impl Fn(usize) -> K,
) -> anyhow::Result<FxHashMap<K, usize>> {
    let role: &Int32Array = column(table, "role")?;
    (0..table.num_rows())
        .filter(|&i| role.is_valid(i))
        .map(|i| Ok((key(i), role_at(role, i, usize::MAX)?)))
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

/// The fitted cells (those the model gave links, in `{m}_block_cells`) in
/// each cell role.
fn fitted_sizes(
    cells: &FxHashMap<CellKey, usize>,
    fitted: &RecordBatch,
    n: usize,
) -> anyhow::Result<Vec<usize>> {
    let (batch, barcode): (&StringArray, &StringArray) =
        (column(fitted, "batch")?, column(fitted, "barcode")?);
    let fitted: FxHashSet<(&str, &str)> = (0..fitted.num_rows())
        .map(|i| (batch.value(i), barcode.value(i)))
        .collect();
    let roles = cells
        .iter()
        .filter(|((b, c), _)| fitted.contains(&(b.as_ref(), c.as_ref())))
        .map(|(_, &r)| r);
    Ok(role_sizes(roles, n))
}

/// The modules' sites from `{m}_block_sites`' rows: each gene a member of
/// its `main` gene role only, its reads there in each cell module.
fn from_sites(
    rows: &BlockSites,
    main: &FxHashMap<&str, usize>,
    (cell_roles, gene_roles): (usize, usize),
) -> anyhow::Result<FromSites> {
    // Each row's cell and gene role, checked once.
    let roles: Vec<(usize, usize)> = (0..rows.pos.len())
        .map(|i| {
            let c = role_at(rows.cell_role, i, cell_roles).context("block_sites cell_role")?;
            let g = role_at(rows.gene_role, i, gene_roles).context("block_sites gene_role")?;
            Ok((c, g))
        })
        .collect::<anyhow::Result<_>>()?;
    // Only the rows of a gene's main role are kept: index their sites alone.
    let kept: Vec<usize> = (0..rows.pos.len())
        .filter(|&i| main.get(rows.gene.value(i)) == Some(&roles[i].1))
        .collect();
    let mut keys: Vec<(&str, &str, i64)> = kept.iter().map(|&i| rows.key(i)).collect();
    keys.sort_unstable();
    keys.dedup();
    let index: FxHashMap<(&str, &str, i64), usize> =
        keys.iter().enumerate().map(|(i, &k)| (k, i)).collect();
    let sites: Vec<Site> = keys
        .iter()
        .map(|&(gene, chr, pos)| Site {
            gene: gene.into(),
            chr: chr.into(),
            pos,
        })
        .collect();
    let mut by_role: Vec<FxHashMap<usize, Vec<Count>>> = vec![FxHashMap::default(); gene_roles];
    for i in kept {
        let (c, g) = roles[i];
        by_role[g]
            .entry(index[&rows.key(i)])
            .or_insert_with(|| vec![Count::default(); cell_roles])[c]
            .add(Count {
                on: rows.converted.value(i) as f64,
                total: rows.total.value(i) as f64,
            });
    }
    let members: Vec<Vec<Member>> = by_role
        .into_iter()
        .map(|of| deepest_first(of, &sites))
        .collect();
    let genes = members
        .iter()
        .map(|list| gene_spans(list, &sites))
        .collect();
    let blocks = members
        .iter()
        .map(|list| block_sums(list, cell_roles))
        .collect();
    Ok(FromSites {
        sites,
        members,
        genes,
        blocks,
    })
}

/// A gene role's members, its genes with the most reads first, each gene's
/// sites by position (sites are numbered so).
fn deepest_first(of: FxHashMap<usize, Vec<Count>>, sites: &[Site]) -> Vec<Member> {
    let mut gene_reads: FxHashMap<&str, f64> = FxHashMap::default();
    for (&s, reads) in &of {
        *gene_reads.entry(&sites[s].gene).or_default() +=
            reads.iter().map(|r| r.total).sum::<f64>();
    }
    // Each member's gene's reads looked up once, not in every comparison.
    let mut keyed: Vec<(f64, Member)> = of
        .into_iter()
        .map(|(site, reads)| {
            (
                gene_reads[sites[site].gene.as_ref()],
                Member { site, reads },
            )
        })
        .collect();
    keyed.sort_unstable_by(|(ra, a), (rb, b)| {
        let (ga, gb) = (&sites[a.site].gene, &sites[b.site].gene);
        rb.total_cmp(ra).then(ga.cmp(gb)).then(a.site.cmp(&b.site))
    });
    keyed.into_iter().map(|(_, m)| m).collect()
}

/// Each gene's span of `list`, whose genes' members are together.
fn gene_spans(list: &[Member], sites: &[Site]) -> Vec<Range<usize>> {
    let mut spans: Vec<Range<usize>> = Vec::new();
    for (i, m) in list.iter().enumerate() {
        match spans.last_mut() {
            Some(span) if sites[list[span.start].site].gene == sites[m.site].gene => {
                span.end = i + 1
            }
            _ => spans.push(i..i + 1),
        }
    }
    spans
}

/// `list`'s reads in each cell module, summed.
fn block_sums(list: &[Member], cell_roles: usize) -> Vec<Count> {
    let mut out = vec![Count::default(); cell_roles];
    for m in list {
        for (o, r) in out.iter_mut().zip(&m.reads) {
            o.add(*r);
        }
    }
    out
}

/// Each site's counts from a `{m}_site_stats` table; a site it lacks
/// counts nothing.
fn site_counts(sites: &[Site], t: &RecordBatch) -> anyhow::Result<SiteStats> {
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

/// The models of the modules' genes among `all`: by key, else by symbol
/// (ASCII case aside, as the pileup matches them), on the gene's chromosome.
pub fn gene_models(m: &Modules, all: &[GeneModel]) -> FxHashMap<Box<str>, GeneModel> {
    let by_key: FxHashMap<&str, &GeneModel> = all.iter().map(|g| (g.key.as_ref(), g)).collect();
    let by_symbol: FxHashMap<String, &GeneModel> = all
        .iter()
        .map(|g| (g.symbol.to_ascii_lowercase(), g))
        .collect();
    let mut out = FxHashMap::default();
    for s in &m.sites {
        if out.contains_key(&s.gene) {
            continue;
        }
        let found = by_key
            .get(s.gene.as_ref())
            .or_else(|| by_symbol.get(&query_symbol(&s.gene).to_ascii_lowercase()));
        if let Some(&g) = found.filter(|g| chr_eq(&g.chr, &s.chr)) {
            out.insert(s.gene.clone(), g.clone());
        }
    }
    out
}
