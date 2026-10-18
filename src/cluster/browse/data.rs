//! What `faba cluster`'s blockmodel wrote, read for the browser: each
//! gene module's sites with their posterior-mean reads in each cell module,
//! each site's over every gene module (the full pileups), the shifts, the
//! cell modules' posterior sizes, the genes' main roles, and each site's
//! counts.

use super::order::{diagonal, standardized};
use crate::cluster::read::{Site, SiteStats};
use crate::cluster::unified::Recorded;
use crate::common::{column, read_record_batch};
use crate::site_analysis::miami::genemodel::GeneModel;
use crate::site_analysis::pileup::query_symbol;
use anyhow::Context;
use arrow::array::StringArray;
use genomic_data::coordinates::chr_eq;
use rustc_hash::FxHashMap;
use sites::from_sites;
use std::ops::Range;
use std::path::Path;
use tables::{main_roles, posterior_sizes, read_shifts, role_sizes, site_counts};
use tables::{BlockSites, SiteRows};

mod run;
mod sites;
mod tables;
mod trace;

pub(super) use run::model as run_model;
pub use trace::Fit;

#[cfg(test)]
pub(super) use tables::role_at;

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

impl std::iter::Sum for Count {
    fn sum<I: Iterator<Item = Count>>(counts: I) -> Count {
        counts.fold(Count::default(), |mut sum, r| {
            sum.add(r);
            sum
        })
    }
}

/// One site of a gene module, with its posterior-mean reads in its gene's
/// main gene role in each cell module, in the panes' column order (see
/// [`Modules::columns`]).
#[derive(Debug, Clone, PartialEq)]
pub struct Member {
    pub site: usize,
    pub reads: Vec<Count>,
}

/// What tables older than this version leave approximate.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct Approximate {
    /// The cell modules' sizes: with no `fitted` column in the cell roles,
    /// the fitted cells are those `{m}_block_cells` lists, which leaves out
    /// cells with under half a link in every block.
    pub sizes: bool,
    /// The full pileups: with no `{m}_site_cells`, each site's reads over
    /// every gene role sum `{m}_block_sites`, which holds no entry under a
    /// read.
    pub pileup: bool,
}

impl Approximate {
    /// What approximate sizes in `modality`'s tables mean; none when exact.
    pub fn sizes_note(&self, modality: &str) -> Option<String> {
        self.sizes.then(|| {
            format!(
                "approximate cell-module sizes: no `fitted` column in {modality}_cell_roles, \
                 as in older tables (rerun with --refit)"
            )
        })
    }

    /// What approximate full pileups in `modality`'s tables mean; none when
    /// exact.
    pub fn pileup_note(&self, modality: &str) -> Option<String> {
        self.pileup.then(|| {
            format!(
                "approximate full pileups: no {modality}_site_cells, as in older tables \
                 (rerun with --refit)"
            )
        })
    }
}

/// A gene's full pileup (level 4): its sites on the chromosome with most of
/// them (the first by name on a tie), and the chromosomes left out.
#[derive(Debug)]
pub struct Pileup<'a> {
    pub chr: &'a str,
    pub sites: Range<usize>,
    pub dropped: Vec<&'a str>,
}

/// What `faba cluster`'s blockmodel wrote, as the browser reads it.
#[derive(Debug)]
pub struct Modules {
    /// The sites' modality, e.g. `m6a`.
    pub modality: String,
    pub cell_roles: usize,
    pub gene_roles: usize,
    /// Each cell role's posterior size: its fitted cells' memberships
    /// summed (placed cells are not in the posterior tally).
    pub cell_sizes: Vec<f64>,
    /// What the tables, when older, leave approximate.
    pub approximate: Approximate,
    /// Genes whose main role each gene role is.
    pub gene_sizes: Vec<usize>,
    /// Each block's shift, `c * gene_roles + g`.
    shift: Vec<f64>,
    /// The sites of the genes with a main role, by gene, chromosome and
    /// position.
    pub sites: Vec<Site>,
    /// Per site, its posterior-mean reads over every gene role in each cell
    /// module, in column order: its gene's full pileup. Site `s`'s are at
    /// `s * cell_roles ..`.
    whole: Vec<Count>,
    /// Per gene role, the sites of the genes whose main role it is with
    /// reads there, and those reads: genes with the most first, each gene's
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
    /// Per site, its raw reads over every cell read and its cells, from
    /// `{m}_site_stats`; none when the tables predate it.
    pub stats: Option<SiteStats>,
    /// The model the tables come from, from their run record; why not
    /// when there is none (see [`run::model`]).
    pub run: Result<Recorded, String>,
    /// Whether the fit converged, and its trace (see [`Fit`]).
    pub fit: Fit,
}

impl Modules {
    /// The tables of `modality` in `dir`.
    pub fn load(dir: &Path, modality: &str) -> anyhow::Result<Self> {
        let path = |name: &str| dir.join(format!("{modality}_{name}.parquet"));
        let table = |name: &str| {
            let path = path(name);
            read_record_batch(&path).map_err(|e| anyhow::anyhow!("{}: {e}", path.display()))
        };
        // A table older runs did not write: none when it is missing.
        let optional = |name: &str| path(name).is_file().then(|| table(name)).transpose();
        let (cell_roles, gene_roles, shift) = read_shifts(&table("blocks")?)?;
        let (cell_sizes, listed) =
            posterior_sizes(&table("cell_roles")?, || table("block_cells"), cell_roles)?;
        let genes_table = table("gene_roles")?;
        let names: &StringArray = column(&genes_table, "gene")?;
        let main = main_roles(&genes_table, |i| names.value(i), gene_roles)
            .with_context(|| format!("{modality}_gene_roles"))?;
        let gene_sizes = role_sizes(main.values().copied(), gene_roles);
        let rows = table("block_sites")?;
        let cells_rows = optional("site_cells")?;
        let site_cells = cells_rows.as_ref().map(SiteRows::of).transpose()?;
        let approximate = Approximate {
            sizes: listed,
            pileup: site_cells.is_none(),
        };
        let notes = [
            approximate.sizes_note(modality),
            approximate.pileup_note(modality),
        ];
        notes.into_iter().flatten().for_each(|n| log::warn!("{n}"));
        let from = from_sites(
            &BlockSites::of(&rows)?,
            site_cells.as_ref(),
            &main,
            (cell_roles, gene_roles),
        )?;
        let stats = optional("site_stats")?
            .map(|t| site_counts(&from.sites, &t))
            .transpose()?;
        if stats.is_none() {
            log::warn!(
                "no {}: rerun `faba cluster --model blocks` for each site's counts",
                path("site_stats").display()
            );
        }
        let record = run::record(dir, modality);
        let mut modules = Self {
            modality: modality.to_string(),
            cell_roles,
            gene_roles,
            cell_sizes,
            approximate,
            gene_sizes,
            shift,
            sites: from.sites,
            whole: from.whole,
            members: from.members,
            genes: from.genes,
            blocks: from.blocks,
            order: Vec::new(),
            diagonal: Vec::new(),
            columns: (0..cell_roles).collect(),
            stats,
            run: record
                .as_ref()
                .map_err(Clone::clone)
                .and_then(run::settings),
            fit: Fit::read(dir, modality, &record),
        };
        modules.set_orders();
        modules.put_in_columns();
        Ok(modules)
    }

    /// Every member's, site's and block's reads put in column order, once,
    /// so read `k` is column `k`'s.
    fn put_in_columns(&mut self) {
        let columns = &self.columns;
        let permute = |reads: &mut Vec<Count>| {
            *reads = columns.iter().map(|&c| reads[c]).collect();
        };
        for m in self.members.iter_mut().flatten() {
            permute(&mut m.reads);
        }
        self.blocks.iter_mut().for_each(permute);
        let mut scratch = Vec::with_capacity(self.cell_roles);
        for reads in self.whole.chunks_exact_mut(self.cell_roles) {
            scratch.clear();
            scratch.extend(columns.iter().map(|&c| reads[c]));
            reads.copy_from_slice(&scratch);
        }
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

    /// Gene role `g`'s reads in each cell module, in column order.
    pub fn block(&self, g: usize) -> &[Count] {
        &self.blocks[g]
    }

    /// The range of gene role `g`'s shifts across the cell roles.
    pub fn shift_range(&self, g: usize) -> f64 {
        crate::cluster::unified::shift_range(&self.shift, self.cell_roles, self.gene_roles, g)
    }

    /// Gene role `g`'s genes in its order, each a span of its members.
    pub fn genes(&self, g: usize) -> &[Range<usize>] {
        &self.genes[g]
    }

    /// The gene of gene role `g`'s member `i`.
    pub fn gene_of(&self, g: usize, i: usize) -> &str {
        &self.sites[self.members[g][i].site].gene
    }

    /// Site `s`'s reads over every gene role in each cell module, in
    /// column order.
    pub fn whole(&self, s: usize) -> &[Count] {
        let k = self.cell_roles;
        &self.whole[s * k..(s + 1) * k]
    }

    /// `gene`'s full pileup; none when it has no sites.
    pub fn pileup(&self, gene: &str) -> Option<Pileup<'_>> {
        let from = self.sites.partition_point(|s| s.gene.as_ref() < gene);
        let to = self.sites.partition_point(|s| s.gene.as_ref() <= gene);
        // A gene's sites are by chromosome: one run each.
        let runs = runs(from..to, |a, b| self.sites[a].chr == self.sites[b].chr);
        let most = runs
            .iter()
            .enumerate()
            .max_by_key(|(i, run)| (run.len(), std::cmp::Reverse(*i)))?
            .0;
        let chr = |run: &Range<usize>| self.sites[run.start].chr.as_ref();
        Some(Pileup {
            chr: chr(&runs[most]),
            sites: runs[most].clone(),
            dropped: (runs.iter().enumerate())
                .filter(|&(i, _)| i != most)
                .map(|(_, run)| chr(run))
                .collect(),
        })
    }

    /// A posterior number of cells, to a tenth, compact past a thousand;
    /// `≈` before it when the sizes are approximate.
    pub fn cells_label(&self, n: f64) -> String {
        let n = if n < 999.95 {
            format!("{n:.1}")
        } else {
            data_beans::interactive::ui::compact(n)
        };
        match self.approximate.sizes {
            true => format!("≈{n}"),
            false => n,
        }
    }
}

/// `range` in runs: each the first of it and those after it that are
/// `same` as that first.
fn runs(range: Range<usize>, same: impl Fn(usize, usize) -> bool) -> Vec<Range<usize>> {
    let mut out: Vec<Range<usize>> = Vec::new();
    for i in range {
        match out.last_mut() {
            Some(run) if same(run.start, i) => run.end = i + 1,
            _ => out.push(i..i + 1),
        }
    }
    out
}

/// The first and last of `positions`, `(0, 0)` when there are none.
pub(super) fn bounds(positions: impl Iterator<Item = i64> + Clone) -> (i64, i64) {
    (
        positions.clone().min().unwrap_or(0),
        positions.max().unwrap_or(0),
    )
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
