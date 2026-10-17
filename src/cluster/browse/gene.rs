//! The fourth level: one gene's full pileup, read from the site matrices
//! and browsed with `faba pileup`'s browser. Its tracks are the levels
//! above it unsummed: all cells (the margin, on a scale of its own), then
//! one per cell module, each converted reads above unconverted. Being read
//! from the matrices, it holds every cell with a cell module, placed and
//! held-out reads too, not only the model's fitted reads.

use super::Modules;
use crate::cluster::read::site_row;
use crate::cluster::read::Input;
use crate::qc::matrix::open_matrix;
use crate::site_analysis::pileup::tui::{show_pileup, Exit, Settings, Track, View};
use crate::site_analysis::pileup::{channel_names, query_symbol, SharedModels};
use rayon::prelude::*;
use rustc_hash::FxHashMap;
use std::sync::{Arc, OnceLock};

/// One channel's reads by position, sorted.
type Channel = Vec<(i64, f64)>;

/// A gene's reads by group, the groups all cells and then each cell module
/// in the order the panes show them: per group, converted and unconverted.
#[derive(Debug, Default, PartialEq)]
pub(super) struct GeneReads {
    pub chr: Box<str>,
    pub groups: Vec<(Channel, Channel)>,
}

/// Per group and position, converted and unconverted reads.
type Tally = Vec<FxHashMap<i64, (f64, f64)>>;

/// A gene's rows in one input: its chromosome, and each row as `(row,
/// position, converted)`.
type GeneRows = (Box<str>, Vec<(usize, i64, bool)>);

/// One input's rows of the modules' modality by gene, and each column's
/// group (its cell module's place among the columns, from 1), none for a
/// cell without a cell module.
struct FileIndex {
    genes: FxHashMap<Box<str>, GeneRows>,
    groups: Vec<Option<usize>>,
}

/// Per batch, each cell's group by its barcode.
type Groups = FxHashMap<Box<str>, FxHashMap<Box<str>, usize>>;

/// The site matrices behind the modules, their rows indexed by gene on a
/// thread of their own as the browser opens, so a gene opened reads only
/// its rows.
pub(super) struct Matrices {
    inputs: Arc<Vec<Input>>,
    /// The index once built, or why it could not be, set by the thread
    /// building it and kept for every ask.
    index: Arc<OnceLock<Result<Vec<FileIndex>, String>>>,
}

impl Matrices {
    /// Index `inputs` for `m` in the background.
    pub(super) fn start(inputs: Vec<Input>, m: &Modules) -> Self {
        let inputs = Arc::new(inputs);
        let group_of_role: FxHashMap<usize, usize> = m
            .columns
            .iter()
            .enumerate()
            .map(|(k, &c)| (c, k + 1))
            .collect();
        let mut groups: Groups = FxHashMap::default();
        for ((batch, barcode), r) in &m.cells {
            if let Some(&k) = group_of_role.get(r) {
                groups
                    .entry(batch.clone())
                    .or_default()
                    .insert(barcode.clone(), k);
            }
        }
        let (files, modality) = (inputs.clone(), m.modality.clone());
        let index = Arc::new(OnceLock::new());
        let done = index.clone();
        std::thread::spawn(move || {
            // Set even when it panics, so no one waits for it for ever.
            let built = std::panic::catch_unwind(std::panic::AssertUnwindSafe(|| {
                files
                    .par_iter()
                    .map(|input| file_index(input, &modality, groups.get(&input.batch)))
                    .collect::<anyhow::Result<Vec<_>>>()
            }));
            let _ = done.set(match built {
                Ok(r) => r.map_err(|e| format!("{e:#}")),
                Err(_) => Err("indexing the matrices failed".to_string()),
            });
        });
        Self { inputs, index }
    }

    /// The index, waiting for it to be built the first time.
    fn index(&self) -> anyhow::Result<&[FileIndex]> {
        self.index
            .wait()
            .as_deref()
            .map_err(|e| anyhow::anyhow!("{e}"))
    }

    /// `gene`'s reads in every input, its cells grouped by their main cell
    /// role; cells without one are left out, so all cells is the sum of the
    /// cell modules.
    pub(super) fn read(&self, m: &Modules, gene: &str) -> anyhow::Result<GeneReads> {
        let index = self.index()?;
        let Some(chr) = index
            .iter()
            .find_map(|f| f.genes.get(gene).map(|g| g.0.clone()))
        else {
            anyhow::bail!("no input holds {}'s sites", query_symbol(gene));
        };
        let groups = 1 + m.cell_roles;
        let parts: Vec<anyhow::Result<Tally>> = self
            .inputs
            .par_iter()
            .zip(index)
            .map(|(input, file)| tally(input, file, gene, groups))
            .collect();
        Ok(GeneReads {
            chr,
            groups: sorted(merge(parts, groups)?),
        })
    }
}

/// The inputs' tallies, added up.
fn merge(parts: Vec<anyhow::Result<Tally>>, groups: usize) -> anyhow::Result<Tally> {
    let mut all: Tally = vec![FxHashMap::default(); groups];
    for part in parts {
        for (into, from) in all.iter_mut().zip(part?) {
            for (pos, (on, off)) in from {
                let e = into.entry(pos).or_default();
                e.0 += on;
                e.1 += off;
            }
        }
    }
    Ok(all)
}

/// Each group's two channels, by position.
fn sorted(all: Tally) -> Vec<(Channel, Channel)> {
    all.into_iter()
        .map(|t| {
            let mut at: Vec<(i64, (f64, f64))> = t.into_iter().collect();
            at.sort_unstable_by_key(|p| p.0);
            at.into_iter()
                .map(|(p, (on, off))| ((p, on), (p, off)))
                .unzip()
        })
        .collect()
}

/// `input`'s rows of `modality` by gene, and its columns' groups from its
/// batch's cells.
fn file_index(
    input: &Input,
    modality: &str,
    groups: Option<&FxHashMap<Box<str>, usize>>,
) -> anyhow::Result<FileIndex> {
    let data = open_matrix(&input.path)?;
    let mut genes: FxHashMap<Box<str>, GeneRows> = FxHashMap::default();
    for (i, name) in data.row_names()?.iter().enumerate() {
        let Some((g, chr, pos, converted)) = site_row(name, modality) else {
            continue;
        };
        genes
            .entry(g.into())
            .or_insert_with(|| (chr.into(), Vec::new()))
            .1
            .push((i, pos, converted));
    }
    let groups = data
        .column_names()?
        .iter()
        .map(|b| groups.and_then(|g| g.get(b)).copied())
        .collect();
    Ok(FileIndex { genes, groups })
}

/// `gene`'s reads in one input by group and position.
fn tally(input: &Input, file: &FileIndex, gene: &str, groups: usize) -> anyhow::Result<Tally> {
    let mut tally: Tally = vec![FxHashMap::default(); groups];
    let Some((_, rows)) = file.genes.get(gene) else {
        return Ok(tally);
    };
    let data = open_matrix(&input.path)?;
    let (_, _, triplets) = data.read_triplets_by_rows(rows.iter().map(|r| r.0).collect())?;
    for (row, col, value) in triplets {
        let Some(Some(k)) = file.groups.get(col as usize) else {
            continue;
        };
        let (_, pos, converted) = rows[row as usize];
        for t in [0, *k] {
            let e = tally[t].entry(pos).or_default();
            if converted {
                e.0 += value as f64;
            } else {
                e.1 += value as f64;
            }
        }
    }
    Ok(tally)
}

/// Browse `gene`'s full pileup until the user goes back (false) or quits
/// (true). Zooming out past the gene widens it in place; a search, which
/// the pileup cannot answer here, says to search from the modules.
/// `settings` carries how the last one was drawn.
pub(super) fn show(
    matrices: &Matrices,
    m: &Modules,
    gene: &str,
    models: Option<&SharedModels>,
    settings: &mut Settings,
) -> anyhow::Result<bool> {
    let symbol = query_symbol(gene);
    let what = format!("reading {symbol} from {} file(s)", matrices.inputs.len());
    let reads = crate::tui::busy(true, "modules", &what, || matrices.read(m, gene))?;
    let (on, off) = channel_names(&m.modality);
    let labels = track_labels(m, on, off);
    let pairs: Vec<(usize, usize)> = (0..reads.groups.len())
        .map(|k| (2 * k, 2 * k + 1))
        .collect();
    let mut extent = extent(&reads);
    let (mut cursor, mut status) = (None, None);
    let keys = [Box::<str>::from(gene)];
    loop {
        let view = View {
            title: &symbol,
            chr: &reads.chr,
            extent,
            on,
            off,
            keys: Some(&keys),
            genes: models.cloned(),
            cursor,
            switchable: false,
            pairs: &pairs,
        };
        match show_pileup(view, tracks(&reads, &labels), status.take(), settings)? {
            Exit::Quit => return Ok(true),
            Exit::Genes => return Ok(false),
            Exit::Locus { lo, hi, at } => (extent, cursor) = ((lo, hi), Some(at)),
            Exit::Search(query) => {
                status = Some(format!(
                    "{query}: search from the modules view (Esc, then /)"
                ))
            }
            Exit::Modality { at } => cursor = Some(at),
        }
    }
}

/// Each group's converted and unconverted tracks, under `labels`.
fn tracks<'a>(reads: &'a GeneReads, labels: &'a [String]) -> Vec<Track<'a>> {
    reads
        .groups
        .iter()
        .enumerate()
        .flat_map(|(k, (up, down))| {
            // All cells sums the cell modules: a scale of its own.
            let track = |label, reads| {
                let t = Track::single(label, "reads", reads, false);
                if k == 0 {
                    t.on_own_scale()
                } else {
                    t
                }
            };
            [track(&labels[2 * k], up), track(&labels[2 * k + 1], down)]
        })
        .collect()
}

/// The first and last position of all cells' converted reads, `(0, 0)` when
/// there are none.
fn extent(reads: &GeneReads) -> (i64, i64) {
    let positions = reads.groups[0].0.iter().map(|p| p.0);
    (
        positions.clone().min().unwrap_or(0),
        positions.max().unwrap_or(0),
    )
}

/// Each group's two tracks' labels: all cells, then each cell module with
/// its cells, converted and then unconverted.
fn track_labels(m: &Modules, on: &str, off: &str) -> Vec<String> {
    let sizes = super::data::role_sizes(m.cells.values().copied(), m.cell_roles);
    std::iter::once(format!("all {} cells", m.cells.len()))
        .chain(
            m.columns
                .iter()
                .map(|&c| format!("C{} {} cells", c + 1, sizes[c])),
        )
        .flat_map(|n| [format!("{n}: {on}"), format!("{n}: {off}")])
        .collect()
}
