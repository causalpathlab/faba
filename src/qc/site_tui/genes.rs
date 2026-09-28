//! The gene list and the selected gene's profile.

use super::*;

/// One gene's sites along its span, all and kept, with its exons.
pub(super) struct GeneProfile {
    pub(super) symbol: String,
    pub(super) chr: Option<String>,
    pub(super) lo: i64,
    pub(super) hi: i64,
    pub(super) forward: Option<bool>,
    pub(super) all: Vec<usize>,
    pub(super) kept: Vec<usize>,
    /// What `all` and `kept` add up.
    pub(super) unit: &'static str,
    /// Per bin, whether it overlaps an exon; `None` without a gene model.
    pub(super) exons: Option<Vec<bool>>,
}

impl GeneProfile {
    /// Base pairs per bin.
    pub(super) fn bin_bp(&self) -> i64 {
        ((self.hi - self.lo) as f64 / self.all.len().max(1) as f64).round() as i64
    }

    pub(super) fn title(&self) -> String {
        let strand = match self.forward {
            Some(true) => " (+)",
            Some(false) => " (-)",
            None => "",
        };
        let chr = self
            .chr
            .as_deref()
            .map_or(String::new(), |c| format!("{c}:"));
        format!("{}  {chr}{}-{}{strand}", self.symbol, self.lo, self.hi)
    }
}

impl<'a> SitePicker<'a> {
    /// List `pinned` (symbols or keys, case-insensitive) first, in the order
    /// given, then every other gene by number of putative sites.
    pub(super) fn pin_genes(mut self, pinned: &[Box<str>]) -> Self {
        let rank = |key: &str, symbol: &str| {
            pinned
                .iter()
                .position(|p| p.eq_ignore_ascii_case(key) || p.eq_ignore_ascii_case(symbol))
        };
        self.gene_order = self
            .views
            .iter()
            .map(|v| {
                let Some(g) = &v.genes else {
                    return Vec::new();
                };
                let mut order: Vec<u32> = (0..g.keys.len() as u32).collect();
                order.sort_by_key(|&i| {
                    let i = i as usize;
                    let pin = rank(&g.keys[i], g.symbol(i)).unwrap_or(usize::MAX);
                    (
                        pin,
                        std::cmp::Reverse(g.rows[i].len()),
                        g.symbol(i).to_string(),
                    )
                });
                order
            })
            .collect();
        self
    }

    /// The focused view's genes that pass the filter, in list order.
    pub(super) fn gene_list(&self) -> Vec<u32> {
        let Some(g) = &self.view().genes else {
            return Vec::new();
        };
        let find = self.gene_find.to_lowercase();
        self.gene_order[self.modality]
            .iter()
            .copied()
            .filter(|&i| find.is_empty() || g.symbol(i as usize).to_lowercase().contains(&find))
            .collect()
    }

    /// The selected gene's dense id.
    pub(super) fn gene(&self) -> Option<usize> {
        self.gene_list().get(self.gene_at).map(|&g| g as usize)
    }

    pub(super) fn step_gene(&mut self, delta: isize) {
        let n = self.gene_list().len() as isize;
        if n > 0 {
            self.gene_at = (self.gene_at as isize + delta).clamp(0, n - 1) as usize;
        }
    }

    /// Kept and all sites of gene `g` in the focused view.
    pub(super) fn gene_kept(&self, g: usize) -> (usize, usize) {
        let Some(genes) = &self.view().genes else {
            return (0, 0);
        };
        let rows = &genes.rows[g];
        let fails = &self.view().fails;
        let kept = rows.iter().filter(|&&i| fails[i as usize] == 0).count();
        (kept, rows.len())
    }

    /// The selected gene's sites along its span in `n` bins: all, kept,
    /// and which bins overlap an exon (`None` without a gene model).
    pub(super) fn gene_profile(&self, n: usize) -> Option<GeneProfile> {
        let g = self.gene()?;
        let genes = self.view().genes.as_ref()?;
        let rows = &genes.rows[g];
        let model = match &self.meta {
            Meta::Ready(a) => a.models.get(&genes.keys[g]),
            _ => None,
        };
        let pos = |i: &u32| genes.pos[*i as usize];
        // The model's span, widened to every site: an annotation other than
        // the one the sites were called on may not cover them all.
        let (min, max) = (rows.iter().map(pos).min()?, rows.iter().map(pos).max()? + 1);
        let (lo, hi) = model.map_or((min, max), |m| (m.lo.min(min), m.hi.max(max)));
        let n = n.max(1);
        let span = (hi - lo).max(1);
        let bin = |p: i64| (((p - lo).clamp(0, span - 1) * n as i64) / span) as usize;
        let (mut all, mut kept) = (vec![0usize; n], vec![0usize; n]);
        let fails = &self.view().fails;
        for i in rows {
            let b = bin(pos(i));
            let w = self.site_weight(*i as usize);
            all[b] += w;
            if fails[*i as usize] == 0 {
                kept[b] += w;
            }
        }
        let exons = model.map(|m| {
            (0..n as i64)
                .map(|b| {
                    let (a, z) = (lo + b * span / n as i64, lo + (b + 1) * span / n as i64);
                    m.exons.iter().any(|&(s, e)| s < z.max(a + 1) && e > a)
                })
                .collect()
        });
        Some(GeneProfile {
            symbol: genes.symbol(g).to_string(),
            chr: model.map(|m| m.chr.to_string()),
            lo,
            hi,
            forward: model.map(|m| m.forward),
            all,
            kept,
            unit: self.weight.unit(),
            exons,
        })
    }
}
