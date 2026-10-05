//! The gene list and the selected gene's profile.

use super::*;
use crate::site_analysis::miami::bin::BinEdges;
use crate::site_analysis::pileup::{gene_matches, query_symbol};

/// The gene list of the focused view: its order, filter, the genes shown,
/// the selection, and each gene's kept sites under the current thresholds.
#[derive(Default)]
pub(super) struct GeneList {
    /// Per view, gene ids in list order: pinned first, then by kept sites.
    order: Vec<Vec<u32>>,
    /// Per view, each gene's place among the pinned genes, `usize::MAX` for
    /// the others.
    pin: Vec<Vec<usize>>,
    /// Per view, each gene's symbol in lowercase, for the filter and ties.
    lower: Vec<Vec<Box<str>>>,
    /// The filter, matched against symbols, case-insensitive.
    pub(super) find: String,
    /// The focused view's genes passing the filter, in list order.
    shown: Vec<u32>,
    /// Position of the selected gene in `shown`.
    pub(super) at: usize,
    /// The gene selected when the search began, to go back to on Esc.
    before_find: Option<u32>,
    /// The focused view's kept sites per gene id.
    kept: Vec<u32>,
}

impl GeneList {
    /// Every view's genes, with `pinned` (as `faba pileup --genes` matches
    /// them) to be listed first, in the order given.
    fn new(views: &[SiteView], pinned: &[Box<str>]) -> Self {
        let pinned: Vec<(&str, Box<str>)> =
            pinned.iter().map(|q| (&**q, query_symbol(q))).collect();
        let rank = |key: &str| pinned.iter().position(|(q, s)| gene_matches(q, s, key));
        let mut list = Self::default();
        for v in views {
            let Some(g) = &v.genes else {
                list.order.push(Vec::new());
                list.pin.push(Vec::new());
                list.lower.push(Vec::new());
                continue;
            };
            let n = g.keys.len();
            list.order.push((0..n as u32).collect());
            list.pin.push(
                g.keys
                    .iter()
                    .map(|k| rank(k).unwrap_or(usize::MAX))
                    .collect(),
            );
            list.lower
                .push((0..n).map(|i| g.symbol(i).to_lowercase().into()).collect());
        }
        list
    }

    /// Rebuild `shown` for `view` after the filter, the order or the view
    /// changed, keeping the selection on gene `keep` while it is shown.
    fn refilter(&mut self, view: usize, keep: Option<u32>) {
        let find = self.find.to_lowercase();
        let lower = &self.lower[view];
        self.shown = self.order[view]
            .iter()
            .copied()
            .filter(|&i| lower[i as usize].contains(&*find))
            .collect();
        self.at = keep
            .and_then(|g| self.shown.iter().position(|&s| s == g))
            .unwrap_or(self.at)
            .min(self.shown.len().saturating_sub(1));
    }

    /// Recount each gene's kept sites in view `v`, and order its genes by
    /// them: pinned first, then most kept, most putative, then by symbol.
    fn recount(&mut self, v: usize, view: &SiteView) {
        let Some(g) = &view.genes else {
            self.kept.clear();
            self.shown.clear();
            return;
        };
        self.kept = g
            .rows
            .iter()
            .map(|rows| {
                rows.iter()
                    .filter(|&&i| view.fails[i as usize] == 0)
                    .count() as u32
            })
            .collect();
        let (pin, lower, kept) = (&self.pin[v], &self.lower[v], &self.kept);
        self.order[v].sort_by(|&a, &b| {
            let key = |i: u32| {
                let i = i as usize;
                (
                    pin[i],
                    std::cmp::Reverse(kept[i]),
                    std::cmp::Reverse(g.rows[i].len()),
                )
            };
            key(a)
                .cmp(&key(b))
                .then_with(|| lower[a as usize].cmp(&lower[b as usize]))
        });
        let keep = self.selected().map(|g| g as u32);
        self.refilter(v, keep);
    }

    pub(super) fn shown(&self) -> &[u32] {
        &self.shown
    }

    pub(super) fn selected(&self) -> Option<usize> {
        self.shown.get(self.at).map(|&g| g as usize)
    }
}

/// One gene's sites along its span, all and kept, with its gene model.
pub(super) struct GeneProfile {
    pub(super) symbol: String,
    pub(super) model: Option<GeneModel>,
    pub(super) edges: BinEdges,
    pub(super) all: Vec<usize>,
    pub(super) kept: Vec<usize>,
    /// What `all` and `kept` add up.
    pub(super) unit: &'static str,
}

impl GeneProfile {
    /// Base pairs per bin.
    fn bin_bp(&self) -> i64 {
        (self.edges.span() as f64 / self.edges.num_bins.max(1) as f64).round() as i64
    }

    /// The gene, where it is, and what one bar holds.
    pub(super) fn title(&self) -> String {
        let (lo, hi) = (self.edges.min_pos, self.edges.max_pos);
        let place = match &self.model {
            Some(m) => format!(
                "{}:{lo}-{hi} ({})",
                m.chr,
                if m.forward { "+" } else { "-" }
            ),
            None => format!("{lo}-{hi}"),
        };
        let sum = |v: &[usize]| v.iter().sum::<usize>();
        format!(
            "{}  {place} · kept {} of {} · y: {} per {} bp",
            self.symbol,
            sum(&self.kept),
            sum(&self.all),
            self.unit,
            self.bin_bp()
        )
    }

    /// Whether bin `b` overlaps an exon; `None` without a gene model.
    pub(super) fn exonic(&self, b: usize) -> Option<bool> {
        let (a, z) = self.edges.col_range(b);
        let m = self.model.as_ref()?;
        Some(m.exons.iter().any(|&(s, e)| s < z.max(a + 1) && e > a))
    }
}

impl<'a> SitePicker<'a> {
    /// List `pinned` genes first; see [`GeneList::new`].
    pub(super) fn pin_genes(mut self, pinned: &[Box<str>]) -> Self {
        self.list = GeneList::new(&self.views, pinned);
        self.list.recount(self.modality, &self.views[self.modality]);
        self
    }

    /// The focused view's genes that pass the filter, in list order.
    pub(super) fn gene_list(&self) -> &[u32] {
        self.list.shown()
    }

    /// The selected gene's dense id.
    pub(super) fn gene(&self) -> Option<usize> {
        self.list.selected()
    }

    pub(super) fn step_gene(&mut self, delta: isize) {
        let n = self.list.shown.len() as isize;
        if n > 0 {
            self.list.at = (self.list.at as isize + delta).clamp(0, n - 1) as usize;
        }
    }

    /// Change the filter and show what passes it, from the top.
    pub(super) fn set_find(&mut self, edit: impl FnOnce(&mut String)) {
        edit(&mut self.list.find);
        self.list.at = 0;
        self.list.refilter(self.modality, None);
    }

    /// Start a gene search, remembering the gene selected now.
    pub(super) fn start_find(&mut self) {
        self.list.before_find = self.list.selected().map(|g| g as u32);
        self.mode = Mode::Find;
    }

    /// End the search on the whole list: with `focus`, on the gene under the
    /// cursor, in the gene panel; else back on the gene selected before.
    pub(super) fn end_find(&mut self, focus: bool) {
        let keep = if focus {
            self.list.selected().map(|g| g as u32)
        } else {
            self.list.before_find
        };
        if focus && keep.is_some() {
            self.panel = Panel::Genes;
        }
        self.list.find.clear();
        self.list
            .refilter(self.modality, keep.or(self.list.before_find));
        self.mode = Mode::Browse;
    }

    /// The gene list for a new focused view, or new thresholds.
    pub(super) fn refresh_genes(&mut self, view_changed: bool) {
        if view_changed {
            self.list.at = 0;
            self.list.shown.clear();
        }
        self.list.recount(self.modality, &self.views[self.modality]);
    }

    /// Kept and all sites of gene `g` in the focused view.
    pub(super) fn gene_kept(&self, g: usize) -> (usize, usize) {
        let all = self
            .view()
            .genes
            .as_ref()
            .map_or(0, |genes| genes.rows[g].len());
        (self.list.kept.get(g).map_or(0, |&k| k as usize), all)
    }

    /// The selected gene's sites along its span in `n` bins, all and kept.
    pub(super) fn gene_profile(&self, n: usize) -> Option<GeneProfile> {
        let g = self.gene()?;
        let genes = self.view().genes.as_ref()?;
        let rows = &genes.rows[g];
        let model = match &self.meta {
            Meta::Ready(a) => a.models.get(&genes.keys[g]).cloned(),
            _ => None,
        };
        let pos = |i: &u32| genes.pos.value(*i as usize);
        // The model's span, widened to every site: an annotation other than
        // the one the sites were called on may not cover them all.
        let (min, max) = (rows.iter().map(pos).min()?, rows.iter().map(pos).max()? + 1);
        let (lo, hi) = model
            .as_ref()
            .map_or((min, max), |m| (m.lo.min(min), m.hi.max(max)));
        let edges = BinEdges::new(lo, hi, n.max(1));
        let (mut all, mut kept) = (vec![0usize; edges.num_bins], vec![0usize; edges.num_bins]);
        let fails = &self.view().fails;
        for i in rows {
            let (b, w) = (edges.col_of(pos(i)), self.site_weight(*i as usize));
            all[b] += w;
            if fails[*i as usize] == 0 {
                kept[b] += w;
            }
        }
        Some(GeneProfile {
            symbol: genes.symbol(g).to_string(),
            model,
            edges,
            all,
            kept,
            unit: self.weight.unit(),
        })
    }
}
