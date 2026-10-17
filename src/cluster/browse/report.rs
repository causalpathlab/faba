//! `faba modules --html`, and `e` in the browser: the modules as one page,
//! the browser's levels and keys in it, its data the module tables (see
//! [`crate::report`]).

use super::Modules;
use crate::site_analysis::miami::genemodel::GeneModel;
use crate::site_analysis::pileup::{channel_names, query_symbol};
use rustc_hash::FxHashMap;
use serde::Serialize;
use std::path::Path;

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
    /// The cell modules' names and sizes, in the panes' column order, as
    /// the terminal shows them.
    cell_names: Vec<String>,
    cell_sizes: Vec<usize>,
    /// The modules' diagonal order for the heatmaps.
    diagonal: Vec<String>,
    modules: Vec<Module>,
}

/// A gene module, its genes in its order.
#[derive(Serialize)]
struct Module {
    name: String,
    size: usize,
    range: f64,
    genes: Vec<Gene>,
}

/// A gene over `len` bases from `lo` (its model's, else its sites' span with
/// a margin), its exons within it.
#[derive(Serialize)]
struct Gene {
    name: String,
    chr: String,
    lo: i64,
    len: i64,
    exons: Vec<(i64, i64)>,
    sites: Vec<Site>,
}

/// A site at `pos` within its gene (the page adds the gene's `lo` for its
/// place on the chromosome): each cell module's `[converted, all]` reads
/// (`c`), in column order, and its counts over every cell read, `[converted, unconverted,
/// cells]` (`s`), when the tables have them.
#[derive(Serialize)]
struct Site {
    pos: i64,
    c: Vec<[f64; 2]>,
    s: Option<[u64; 3]>,
}

/// To a tenth, which is all a page shows.
fn tenth(x: f64) -> f64 {
    (10.0 * x).round() / 10.0
}

/// To a thousandth, enough for a shift range's two decimals.
fn thousandth(x: f64) -> f64 {
    (1000.0 * x).round() / 1000.0
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
        cell_sizes: m.columns.iter().map(|&c| m.cell_sizes[c]).collect(),
        diagonal: m.diagonal.iter().map(|g| format!("M{}", g + 1)).collect(),
        modules,
    }
}

/// The sites `span` of gene module `g` (one gene's run of its members), on
/// that gene's model in `models` when it has one, their reads in column
/// order.
fn gene(
    m: &Modules,
    g: usize,
    span: std::ops::Range<usize>,
    models: &FxHashMap<Box<str>, GeneModel>,
) -> Gene {
    let members = &m.members[g][span];
    let first = &m.sites[members[0].site];
    let at = |i: usize| m.sites[members[i].site].pos;
    let (lo, hi, exons) = match models.get(&first.gene) {
        Some(model) => {
            let exons = model
                .exons
                .iter()
                .map(|&(s, e)| (s - model.lo, e - model.lo))
                .collect();
            (model.lo, model.hi - 1, exons)
        }
        None => {
            let (lo, hi) = (at(0), at(members.len() - 1));
            let pad = ((hi - lo) / 20).max(50);
            (lo - pad, hi + pad, Vec::new())
        }
    };
    Gene {
        name: query_symbol(&first.gene).to_string(),
        chr: first.chr.to_string(),
        lo,
        len: (hi - lo).max(1),
        exons,
        sites: members
            .iter()
            .map(|member| {
                let s = member.site;
                Site {
                    // At the model's ends when outside it, as the strip
                    // draws them, so no reads drop from the gene's tracks.
                    pos: m.sites[s].pos.clamp(lo, hi) - lo,
                    c: member
                        .reads
                        .iter()
                        .map(|r| [tenth(r.on), tenth(r.total)])
                        .collect(),
                    s: m.stats.as_ref().map(|t| [t.on[s], t.off[s], t.cells[s]]),
                }
            })
            .collect(),
    }
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
