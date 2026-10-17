//! The modules' sites from `{m}_block_sites` and `{m}_site_cells`: each
//! site's reads over every gene role, and each gene role's members.

use super::tables::{role_at, BlockSites, SiteRows};
use super::{runs, Count, Member};
use crate::cluster::read::Site;
use anyhow::Context;
use rustc_hash::FxHashMap;
use std::ops::Range;

/// A site as the tables key it: its gene, chromosome and position.
type Key<'a> = (&'a str, &'a str, i64);

/// The sites' part of the modules.
pub(super) struct FromSites {
    pub(super) sites: Vec<Site>,
    /// Per site `s`, its reads over every gene role in cell role `c` at
    /// `s * cell_roles + c`.
    pub(super) whole: Vec<Count>,
    pub(super) members: Vec<Vec<Member>>,
    pub(super) genes: Vec<Vec<Range<usize>>>,
    pub(super) blocks: Vec<Vec<Count>>,
}

/// The sites met so far, numbered as met, each with its reads over every
/// gene role in each cell role.
struct Found<'a> {
    cell_roles: usize,
    keys: Vec<Key<'a>>,
    at: FxHashMap<Key<'a>, usize>,
    /// Per site met, its reads in each cell role, as [`FromSites::whole`].
    whole: Vec<Count>,
    /// The last site met, found again without a lookup: the tables are by
    /// site.
    last: Option<(Key<'a>, usize)>,
}

impl<'a> Found<'a> {
    /// None met yet, over `cell_roles` cell roles.
    fn new(cell_roles: usize) -> Self {
        Self {
            cell_roles,
            keys: Vec::new(),
            at: FxHashMap::default(),
            whole: Vec::new(),
            last: None,
        }
    }

    /// `key`'s number, a new one when it is met first.
    fn site(&mut self, key: Key<'a>) -> usize {
        if let Some((last, s)) = self.last {
            if last == key {
                return s;
            }
        }
        let s = *self.at.entry(key).or_insert_with(|| {
            self.keys.push(key);
            self.whole
                .extend(std::iter::repeat_n(Count::default(), self.cell_roles));
            self.keys.len() - 1
        });
        self.last = Some((key, s));
        s
    }

    /// `reads` added to `key`'s in cell role `c`.
    fn add(&mut self, key: Key<'a>, c: usize, reads: Count) {
        let s = self.site(key);
        self.whole[s * self.cell_roles + c].add(reads);
    }

    /// The sites by gene, chromosome and position, their reads in that
    /// order, and the place in it of each site by the number it was met as.
    fn sorted(self) -> (Vec<Site>, Vec<Count>, Vec<usize>) {
        let mut order: Vec<usize> = (0..self.keys.len()).collect();
        order.sort_unstable_by_key(|&s| self.keys[s]);
        let mut rank = vec![0; order.len()];
        for (k, &s) in order.iter().enumerate() {
            rank[s] = k;
        }
        let sites = order
            .iter()
            .map(|&s| {
                let (gene, chr, pos) = self.keys[s];
                Site {
                    gene: gene.into(),
                    chr: chr.into(),
                    pos,
                }
            })
            .collect();
        let k = self.cell_roles;
        let whole = order
            .iter()
            .flat_map(|&s| self.whole[s * k..(s + 1) * k].iter().copied())
            .collect();
        (sites, whole, rank)
    }
}

/// The modules' sites from `{m}_block_sites`' rows and, when there is
/// one, `{m}_site_cells`' (`whole`): every site of the genes with a `main`
/// gene role, each with its reads over every gene role in each cell module
/// (from `whole`, else summed over the block sites' gene roles), and as a
/// member of its gene's main gene role, with its reads there, when it has
/// any.
pub(super) fn from_sites(
    rows: &BlockSites,
    whole: Option<&SiteRows>,
    main: &FxHashMap<&str, usize>,
    (cell_roles, gene_roles): (usize, usize),
) -> anyhow::Result<FromSites> {
    let roles = row_roles(rows, (cell_roles, gene_roles))?;
    // Per row, whether it is in its gene's main role; none when its gene
    // has no main role.
    let in_main: Vec<Option<bool>> = (0..rows.rows.len())
        .map(|i| main.get(rows.rows.gene.value(i)).map(|&g| g == roles[i].1))
        .collect();
    let mut found = Found::new(cell_roles);
    match whole {
        Some(w) => add_site_cells(&mut found, w, main, cell_roles)?,
        None => {
            for i in (0..rows.rows.len()).filter(|&i| in_main[i].is_some()) {
                found.add(rows.rows.key(i), roles[i].0, rows.rows.reads(i));
            }
        }
    }
    let mut by_role: Vec<FxHashMap<usize, Vec<Count>>> = vec![FxHashMap::default(); gene_roles];
    for i in (0..rows.rows.len()).filter(|&i| in_main[i] == Some(true)) {
        let (c, g) = roles[i];
        let s = found.site(rows.rows.key(i));
        let reads = by_role[g]
            .entry(s)
            .or_insert_with(|| vec![Count::default(); cell_roles]);
        reads[c].add(rows.rows.reads(i));
    }
    let (sites, whole, rank) = found.sorted();
    let members: Vec<Vec<Member>> = by_role
        .into_iter()
        .map(|of| deepest_first(of.into_iter().map(|(s, r)| (rank[s], r)), &sites))
        .collect();
    let genes = members
        .iter()
        .map(|list| {
            runs(0..list.len(), |a, b| {
                sites[list[a].site].gene == sites[list[b].site].gene
            })
        })
        .collect();
    let blocks = members
        .iter()
        .map(|list| block_sums(list, cell_roles))
        .collect();
    Ok(FromSites {
        sites,
        whole,
        members,
        genes,
        blocks,
    })
}

/// Each block site row's cell and gene role, checked once.
fn row_roles(
    rows: &BlockSites,
    (cell_roles, gene_roles): (usize, usize),
) -> anyhow::Result<Vec<(usize, usize)>> {
    (0..rows.rows.len())
        .map(|i| {
            let c = role_at(rows.rows.cell_role, i, cell_roles).context("block_sites cell_role")?;
            let g = role_at(rows.gene_role, i, gene_roles).context("block_sites gene_role")?;
            Ok((c, g))
        })
        .collect()
}

/// `w`'s rows of the genes with a `main` role added to `found`.
fn add_site_cells<'a>(
    found: &mut Found<'a>,
    w: &SiteRows<'a>,
    main: &FxHashMap<&str, usize>,
    cell_roles: usize,
) -> anyhow::Result<()> {
    // A gene's rows are together: its lookup is reused.
    let mut last: Option<(&str, bool)> = None;
    for i in 0..w.len() {
        let gene = w.gene.value(i);
        let kept = match last {
            Some((g, kept)) if g == gene => kept,
            _ => {
                let kept = main.contains_key(gene);
                last = Some((gene, kept));
                kept
            }
        };
        if kept {
            let c = role_at(w.cell_role, i, cell_roles).context("site_cells cell_role")?;
            found.add(w.key(i), c, w.reads(i));
        }
    }
    Ok(())
}

/// A gene role's members `of` (each a site and its reads), its genes with
/// the most reads first, each gene's sites by position (sites are numbered
/// so).
fn deepest_first(of: impl Iterator<Item = (usize, Vec<Count>)>, sites: &[Site]) -> Vec<Member> {
    let of: Vec<(usize, Vec<Count>)> = of.collect();
    let mut gene_reads: FxHashMap<&str, f64> = FxHashMap::default();
    for (s, reads) in &of {
        *gene_reads.entry(&sites[*s].gene).or_default() +=
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

/// `list`'s reads in each cell module, summed.
fn block_sums(list: &[Member], cell_roles: usize) -> Vec<Count> {
    (0..cell_roles)
        .map(|c| list.iter().map(|m| m.reads[c]).sum())
        .collect()
}
