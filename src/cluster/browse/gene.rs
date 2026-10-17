//! The fourth level: one gene's full pileup, from the posterior, browsed
//! with `faba pileup`'s browser. Its tracks are the levels above it
//! unsummed: all cells (the margin, on a scale of its own), then one per
//! cell module, each converted reads above unconverted. Each site's reads
//! are its posterior-mean reads over every gene role (`{m}_site_cells`),
//! split across the cell modules by the links' role draws, so the gene's
//! share in its main module (level 3) is a part of it. Placed cells are not
//! in the posterior tally: the tracks hold the fitted cells' reads.

use super::data::bounds;
use super::{Count, Modules};
use crate::site_analysis::pileup::tui::{show_pileup, Exit, Settings, Track, View};
use crate::site_analysis::pileup::{channel_names, query_symbol, SharedModels};

/// One channel's reads by position, sorted.
type Channel = Vec<(i64, f64)>;

/// A gene's reads by group, the groups all cells and then each cell module
/// in the order the panes show them: per group, converted and unconverted;
/// on the chromosome with most of its sites, the others' left out.
#[derive(Debug)]
pub(super) struct GeneReads<'a> {
    pub chr: &'a str,
    pub groups: Vec<(Channel, Channel)>,
    pub dropped: Vec<&'a str>,
}

impl<'a> GeneReads<'a> {
    /// `gene`'s posterior-mean reads at each of its sites, over every gene
    /// role: all cells the sum over the cell modules.
    pub(super) fn of(m: &'a Modules, gene: &str) -> anyhow::Result<Self> {
        let Some(pileup) = m.pileup(gene) else {
            anyhow::bail!("the modules hold no sites of {}", query_symbol(gene));
        };
        let mut groups = vec![(Channel::new(), Channel::new()); 1 + m.cell_roles];
        for s in pileup.sites {
            let (pos, whole) = (m.sites[s].pos, m.whole(s));
            let all: Count = whole.iter().copied().sum();
            for (group, r) in groups.iter_mut().zip(std::iter::once(&all).chain(whole)) {
                group.0.push((pos, r.on));
                group.1.push((pos, r.off()));
            }
        }
        Ok(Self {
            chr: pileup.chr,
            groups,
            dropped: pileup.dropped,
        })
    }
}

/// What the footer says first of `reads`: the chromosomes left out, and
/// that older tables make the pileup approximate; none when neither.
pub(super) fn note(m: &Modules, reads: &GeneReads) -> Option<String> {
    let mut parts = Vec::new();
    if !reads.dropped.is_empty() {
        parts.push(format!(
            "its sites on {} left out (shown: {})",
            reads.dropped.join(", "),
            reads.chr
        ));
    }
    parts.extend(m.approximate.pileup_note(&m.modality));
    (!parts.is_empty()).then(|| parts.join(" · "))
}

/// Browse `gene`'s full pileup until the user goes back (false) or quits
/// (true). Zooming out past the gene widens it in place; a search, which
/// the pileup cannot answer here, says to search from the modules.
/// `settings` carries how the last one was drawn; when not every cell
/// module's row fits, they start at the one in column `first`.
pub(super) fn show(
    m: &Modules,
    (gene, first): (&str, usize),
    models: Option<&SharedModels>,
    settings: &mut Settings,
) -> anyhow::Result<bool> {
    let reads = GeneReads::of(m, gene)?;
    let title = format!("{} · {POSTERIOR}", query_symbol(gene));
    let (on, off) = channel_names(&m.modality);
    let labels = track_labels(m, on, off);
    let keys = [Box::<str>::from(gene)];
    let mut place = Place {
        extent: extent(&reads),
        cursor: None,
        // Group 0 is all cells; each group has two tracks.
        track_from: 2 * (first + 1),
        status: note(m, &reads),
    };
    loop {
        let view = View {
            title: &title,
            chr: reads.chr,
            extent: place.extent,
            on,
            off,
            keys: Some(&keys),
            genes: models.cloned(),
            cursor: place.cursor,
            switchable: false,
            track_from: place.track_from,
        };
        let status = place.status.take();
        let (exit, track_from) = show_pileup(view, tracks(&reads, &labels), status, settings)?;
        if let Some(quit) = place.left(exit, track_from) {
            return Ok(quit);
        }
    }
}

/// Where a gene's browser opens again after leaving for what it cannot do
/// in place.
#[derive(Debug, Default, PartialEq)]
pub(super) struct Place {
    /// The span shown, widened by `-` past the gene.
    pub(super) extent: (i64, i64),
    /// Where the cursor was; the first site when none.
    pub(super) cursor: Option<i64>,
    /// The track scrolled to.
    pub(super) track_from: usize,
    /// What the footer says first.
    pub(super) status: Option<String>,
}

impl Place {
    /// Left by `exit`, scrolled to `track_from`: whether the user quit
    /// (`Some(true)`) or went back (`Some(false)`), else where to open again.
    pub(super) fn left(&mut self, exit: Exit, track_from: usize) -> Option<bool> {
        self.track_from = track_from;
        match exit {
            Exit::Quit => return Some(true),
            Exit::Genes => return Some(false),
            Exit::Locus { lo, hi, at } => (self.extent, self.cursor) = ((lo, hi), Some(at)),
            Exit::Search(query) => {
                self.status = Some(format!(
                    "{query}: search from the modules view (Esc, then /)"
                ))
            }
            Exit::Modality { at } => self.cursor = Some(at),
        }
        None
    }
}

/// Each group's converted and unconverted tracks, under `labels`, mirrored
/// in one row: converted above, in the accent.
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
            let up = track(&labels[2 * k], up).converted_above();
            [up, track(&labels[2 * k + 1], down)]
        })
        .collect()
}

/// The first and last position of all cells' converted reads, `(0, 0)` when
/// there are none.
fn extent(reads: &GeneReads) -> (i64, i64) {
    bounds(reads.groups[0].0.iter().map(|p| p.0))
}

/// What the full pileup's numbers are, in its title.
const POSTERIOR: &str = "posterior-mean reads, fitted cells";

/// Each group's two tracks' labels: all cells, then each cell module with
/// its posterior size, converted and then unconverted.
fn track_labels(m: &Modules, on: &str, off: &str) -> Vec<String> {
    let all: f64 = m.cell_sizes.iter().sum();
    std::iter::once(format!("all {} cells", m.cells_label(all)))
        .chain(
            m.columns
                .iter()
                .map(|&c| format!("C{} {} cells", c + 1, m.cells_label(m.cell_sizes[c]))),
        )
        .flat_map(|n| [format!("{n}: {on}"), format!("{n}: {off}")])
        .collect()
}
