//! The Metropolis-Hastings swaps: two roles exchanged on all of a cell's
//! or a gene's links, each swap its own inverse.

use super::{Model, State, Urn};
use rand::rngs::SmallRng;
use rand::RngExt;

/// Buffers a sweep's moves reuse: the relabelling's log weights and
/// running products per role, or per block when drawn jointly (see [`Model::ln_link_blocks`]), and a swap's
/// moved links and the (site, block) entries they touch.
pub(super) struct Scratch {
    pub(super) w: Vec<f64>,
    pub(super) rising: Vec<(f64, i64)>,
    moved: Vec<usize>,
    touched: Vec<usize>,
}

impl Scratch {
    /// Empty buffers sized for `model`.
    pub(super) fn new(model: &Model) -> Self {
        let (aa, bb) = (model.cells.roles, model.settings.gene_roles);
        let roles = if model.joint_relabel() {
            aa * bb
        } else {
            aa.max(bb)
        };
        Self {
            w: vec![0.0; roles],
            rising: vec![(1.0, 0); roles],
            moved: Vec::new(),
            touched: Vec::new(),
        }
    }
}

impl Model<'_> {
    /// Cell `unit`'s links (`cell` true) or gene `unit`'s that hold role
    /// `r` or `r2`, into `out`.
    fn holding(
        &self,
        st: &State,
        (unit, cell): (usize, bool),
        (r, r2): (u16, u16),
        out: &mut Vec<usize>,
    ) {
        out.clear();
        let holds = |&l: &usize| {
            let x = if cell { st.a[l] } else { st.b[l] };
            x == r || x == r2
        };
        if cell {
            out.extend((self.links.cell_ptr[unit]..self.links.cell_ptr[unit + 1]).filter(holds));
        } else {
            out.extend(
                self.gene_link[self.gene_ptr[unit]..self.gene_ptr[unit + 1]]
                    .iter()
                    .copied()
                    .filter(holds),
            );
        }
    }

    /// The change in a cell's urn (`cell` true) or a gene's as roles `r`
    /// and `r2` exchange their links.
    fn unit_exchange(
        &self,
        urn: &Urn,
        st: &State,
        (unit, cell): (usize, bool),
        (r, r2): (u16, u16),
    ) -> f64 {
        let pair = (r as usize, r2 as usize);
        if cell {
            self.cells.exchange(urn, unit, pair)
        } else {
            self.gene_side.exchange(&st.genes, unit, pair)
        }
    }

    /// Roles `r` and `r2` exchanged on links `moved` (each holding one of
    /// them) of a cell (`cell` true) or a gene, their counts moved with
    /// them: its own inverse.
    fn flip(
        &self,
        urn: &mut Urn,
        st: &mut State,
        moved: &[usize],
        cell: bool,
        (r, r2): (u16, u16),
    ) {
        for &l in moved {
            self.put(urn, st, l, false);
            let x = if cell { &mut st.a[l] } else { &mut st.b[l] };
            *x = if *x == r { r2 } else { r };
        }
        for &l in moved {
            self.put(urn, st, l, true);
        }
    }

    /// Exchanges roles `r` and `r2` on cell `unit`'s links (`cell` true) or
    /// gene `unit`'s, in place (see [`Self::flip`]), and returns the change
    /// in log joint, the links moved left in `scratch`: the unit's urn
    /// (counted once, whatever holds its links) and this model's part (see
    /// [`Self::exchange_links`]); none when no link holds either role.
    pub(super) fn exchange(
        &self,
        urn: &mut Urn,
        st: &mut State,
        (unit, cell): (usize, bool),
        pair: (u16, u16),
        scratch: &mut Scratch,
    ) -> Option<f64> {
        self.holding(st, (unit, cell), pair, &mut scratch.moved);
        if scratch.moved.is_empty() {
            return None;
        }
        let unit_change = self.unit_exchange(urn, st, (unit, cell), pair);
        let (sites_before, sites_after, reads) = self.exchange_links(urn, st, cell, pair, scratch);
        Some(unit_change + (sites_after - sites_before) + reads)
    }

    /// This model's part of [`Self::exchange`]: roles `r` and `r2`
    /// exchanged, in place, on the links in `scratch`'s moved list (a cell's,
    /// `cell` true, or a gene's, each holding one of them), and the site
    /// terms of every block a moved link leaves or enters before and after,
    /// or at infinite strength (the reads factoring by link, each at its
    /// cell's level) each moved link's reads in the block it enters less
    /// those in the block it leaves.
    fn exchange_links(
        &self,
        urn: &mut Urn,
        st: &mut State,
        cell: bool,
        (r, r2): (u16, u16),
        scratch: &mut Scratch,
    ) -> (f64, f64, f64) {
        let (bb, kk) = (self.settings.gene_roles, self.blocks());
        let Scratch { moved, touched, .. } = scratch;
        // The block a link enters from block `k`.
        let other = |k: usize| {
            let (a, b) = ((k / bb) as u16, (k % bb) as u16);
            let swap = |x: u16| if x == r { r2 } else { r };
            if cell {
                swap(a) as usize * bb + b as usize
            } else {
                a as usize * bb + swap(b) as usize
            }
        };
        touched.clear();
        let mut reads = 0.0;
        for &l in moved.iter() {
            let from = self.block(st, l);
            let to = other(from);
            if self.finite() {
                for &s in self.links.of(l).0 {
                    touched.extend([s as usize * kk + from, s as usize * kk + to]);
                }
            } else {
                reads += self.ln_link(st, l, to) - self.ln_link(st, l, from);
            }
        }
        touched.sort_unstable();
        touched.dedup();
        let sites = |st: &State| {
            touched
                .iter()
                .map(|&at| self.site_term(st, at))
                .sum::<f64>()
        };
        let before = sites(st);
        self.flip(urn, st, moved, cell, (r, r2));
        (before, sites(st), reads)
    }

    /// A swap of two roles drawn uniformly on all of cell `unit`'s links
    /// (`cell` true) or gene `unit`'s: its own inverse, so the proposal is
    /// symmetric and the Metropolis-Hastings ratio the change in log joint.
    /// Rejected, the labels are flipped back, unscored.
    pub(super) fn swap(
        &self,
        urn: &mut Urn,
        st: &mut State,
        (unit, cell): (usize, bool),
        scratch: &mut Scratch,
        rng: &mut SmallRng,
    ) {
        let roles = if cell {
            self.cells.roles
        } else {
            self.settings.gene_roles
        };
        if roles < 2 {
            return;
        }
        let r = rng.random_range(0..roles);
        let r2 = (r + rng.random_range(1..roles)) % roles;
        let pair = (r as u16, r2 as u16);
        if let Some(delta) = self.exchange(urn, st, (unit, cell), pair, scratch) {
            if rng.random::<f64>().ln() >= delta {
                self.flip(urn, st, &scratch.moved, cell, pair);
            }
        }
    }
}
