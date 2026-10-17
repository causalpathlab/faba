//! `faba modules --html`, and `e` in the browser: the modules as one page,
//! the browser's levels and keys in it, its data the module tables (see
//! [`crate::report`]).

use super::data::bounds;
use super::{Count, Modules};
use crate::report::{tenth, thousandth};
use crate::site_analysis::miami::genemodel::GeneModel;
use crate::site_analysis::pileup::{channel_names, query_symbol};
use rustc_hash::{FxHashMap, FxHashSet};
use serde::Serialize;
use std::path::Path;

mod plate;

/// The page's own markup, styles and code.
pub(super) const FRAGMENT: &str = include_str!("report.html");

/// What the page reads from `#data`.
#[derive(Serialize)]
#[serde(rename_all = "camelCase")]
pub(super) struct Data {
    title: String,
    modality: String,
    on: &'static str,
    off: &'static str,
    /// The cell modules' names and posterior sizes (their fitted cells'
    /// memberships summed), in the panes' column order, the sizes labelled
    /// as the terminal labels them (see [`Modules::cells_label`]), and all
    /// of theirs, the exact sum labelled once.
    cell_names: Vec<String>,
    cell_labels: Vec<String>,
    all_label: String,
    /// What older tables leave approximate in the full pileups (see
    /// [`super::data::Approximate`]); sent only when they do.
    #[serde(skip_serializing_if = "Option::is_none")]
    pileup_note: Option<String>,
    /// The modules' diagonal order for the heatmaps.
    diagonal: Vec<String>,
    modules: Vec<Module>,
    /// The model the tables come from, drawn, or why it is not.
    model: plate::Section,
}

/// A gene module, its genes in its order.
#[derive(Serialize)]
struct Module {
    name: String,
    size: usize,
    range: f64,
    genes: Vec<Gene>,
}

/// A gene on `chr` (its full pileup's) over `len` bases from `lo` (its
/// model's, else its sites' span with a margin), its exons within it: its
/// sites in its module (levels 2 and 3, and 4 unless `o`), and those in its
/// full pileup only (`more`).
#[derive(Serialize)]
struct Gene {
    name: String,
    chr: String,
    lo: i64,
    len: i64,
    exons: Vec<(i64, i64)>,
    sites: Vec<Site>,
    #[serde(skip_serializing_if = "Vec::is_empty")]
    more: Vec<Whole>,
}

/// A site at `pos` within its gene (the page adds the gene's `lo` for its
/// place on the chromosome): each cell module's posterior-mean `[converted,
/// all]` reads in column order in its gene's main gene role (`c`), and over
/// every gene role (`w`, level 4) when they differ; all cell modules' of
/// each (`ca`, `wa`, the exact sums rounded once) where the page, adding
/// the cell modules' up, would not get the same (see [`unless`]); `o` when
/// it is off its full pileup's chromosome; and its counts over every cell
/// read, `[converted, unconverted, cells]` (`s`), when the tables have them.
#[derive(Serialize)]
struct Site {
    pos: i64,
    c: Vec<[f64; 2]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    w: Option<Vec<[f64; 2]>>,
    #[serde(skip_serializing_if = "Option::is_none")]
    ca: Option<[f64; 2]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    wa: Option<[f64; 2]>,
    #[serde(skip_serializing_if = "std::ops::Not::not")]
    o: bool,
    s: Option<[u64; 3]>,
}

/// A site in its gene's full pileup only (no reads in its main gene role),
/// at `pos` as [`Site`]'s, with its `w` and `wa`.
#[derive(Serialize)]
struct Whole {
    pos: i64,
    w: Vec<[f64; 2]>,
    #[serde(skip_serializing_if = "Option::is_none")]
    wa: Option<[f64; 2]>,
}

/// The page's data: the modules in their order, as the browser shows them.
pub(super) fn data(m: &Modules, models: &FxHashMap<Box<str>, GeneModel>) -> Data {
    let modules = m
        .order
        .iter()
        .map(|&g| Module {
            name: format!("M{}", g + 1),
            size: m.gene_sizes[g],
            range: thousandth(m.shift_range(g)),
            genes: m
                .genes(g)
                .iter()
                .map(|span| gene(m, g, span.clone(), models))
                .collect(),
        })
        .collect();
    let (on, off) = channel_names(&m.modality);
    Data {
        title: format!("gene modules × cell modules · {}", m.modality),
        modality: m.modality.clone(),
        on,
        off,
        cell_names: m.columns.iter().map(|c| format!("C{}", c + 1)).collect(),
        cell_labels: (m.columns.iter())
            .map(|&c| m.cells_label(m.cell_sizes[c]))
            .collect(),
        all_label: m.cells_label(m.cell_sizes.iter().sum()),
        pileup_note: m.approximate.pileup_note(&m.modality),
        diagonal: m.diagonal.iter().map(|g| format!("M{}", g + 1)).collect(),
        modules,
        model: plate::section(&m.run),
    }
}

/// The sites `span` of gene module `g` (one gene's run of its members), on
/// that gene's model in `models` when it has one, their reads in column
/// order, with the rest of its full pileup.
fn gene(
    m: &Modules,
    g: usize,
    span: std::ops::Range<usize>,
    models: &FxHashMap<Box<str>, GeneModel>,
) -> Gene {
    let members = &m.members[g][span];
    let name = &m.sites[members[0].site].gene;
    let pileup = m.pileup(name).expect("a member's gene has sites");
    let (lo, hi, exons) = match models.get(name) {
        Some(model) => {
            let exons = model
                .exons
                .iter()
                .map(|&(s, e)| (s - model.lo, e - model.lo))
                .collect();
            (model.lo, model.hi - 1, exons)
        }
        None => {
            let at = members.iter().map(|x| x.site).chain(pileup.sites.clone());
            let (lo, hi) = bounds(at.map(|s| m.sites[s].pos));
            let pad = ((hi - lo) / 20).max(50);
            (lo - pad, hi + pad, Vec::new())
        }
    };
    // At the model's ends when outside it, as the strip draws them, so no
    // reads drop from the gene's tracks.
    let pos = |s: usize| m.sites[s].pos.clamp(lo, hi) - lo;
    let sites = members
        .iter()
        .map(|member| {
            let s = member.site;
            let c = tenths(&member.reads);
            let ca = unless(all_tenths(&member.reads), page_sum(&c));
            let in_pileup = pileup.sites.contains(&s);
            let w = Some(tenths(m.whole(s))).filter(|w| in_pileup && *w != c);
            // Without `w` or `wa` the page takes `c`'s and `ca`'s.
            let fallback = match &w {
                Some(w) => page_sum(w),
                None => ca.unwrap_or(page_sum(&c)),
            };
            let wa = unless(all_tenths(m.whole(s)), fallback).filter(|_| in_pileup);
            Site {
                pos: pos(s),
                c,
                w,
                ca,
                wa,
                o: !in_pileup,
                s: m.stats.as_ref().map(|t| [t.on[s], t.off[s], t.cells[s]]),
            }
        })
        .collect();
    let listed: FxHashSet<usize> = members.iter().map(|x| x.site).collect();
    let more = pileup
        .sites
        .filter(|s| !listed.contains(s))
        .map(|s| {
            let w = tenths(m.whole(s));
            Whole {
                pos: pos(s),
                wa: unless(all_tenths(m.whole(s)), page_sum(&w)),
                w,
            }
        })
        .collect();
    Gene {
        name: query_symbol(name).to_string(),
        chr: pileup.chr.to_string(),
        lo,
        len: (hi - lo).max(1),
        exons,
        sites,
        more,
    }
}

/// `r`'s `[converted, all]` reads, to a tenth.
fn pair(r: Count) -> [f64; 2] {
    [tenth(r.on), tenth(r.total)]
}

/// Each cell module's `[converted, all]` reads, to a tenth.
fn tenths(reads: &[Count]) -> Vec<[f64; 2]> {
    reads.iter().copied().map(pair).collect()
}

/// All cell modules' `[converted, all]` reads, summed and then to a tenth.
fn all_tenths(reads: &[Count]) -> [f64; 2] {
    pair(reads.iter().copied().sum())
}

/// `pairs` added up in order from zero, as the page's script adds them.
pub(super) fn page_sum(pairs: &[[f64; 2]]) -> [f64; 2] {
    pairs
        .iter()
        .fold([0.0, 0.0], |a, p| [a[0] + p[0], a[1] + p[1]])
}

/// `exact` when the page, falling back to `page`, would not get it.
fn unless(exact: [f64; 2], page: [f64; 2]) -> Option<[f64; 2]> {
    (exact != page).then_some(exact)
}

/// The modules as a page, gene models drawn from `models`, written to `path`.
pub fn write(
    m: &Modules,
    models: &FxHashMap<Box<str>, GeneModel>,
    path: &Path,
) -> anyhow::Result<()> {
    save(&data(m, models), path)
}

/// The page of `data` (see [`data`]), streamed to `path`.
pub(super) fn save(data: &Data, path: &Path) -> anyhow::Result<()> {
    let title = format!("faba modules · {}", data.modality);
    crate::report::write(path, &title, FRAGMENT, data)
}
