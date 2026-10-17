//! The pinned genes: pinning, unpinning and moving them, and the site in
//! focus along the strip.

use super::{ModulesView, Pin};

impl ModulesView<'_> {
    /// The pin for the gene of the module's member `i`.
    pub(super) fn pin_of(&self, i: usize) -> Option<Pin> {
        let role = self.role()?;
        let gene = self.m.genes(role).iter().position(|s| s.contains(&i))?;
        Some(Pin { role, gene })
    }

    /// Whether member `i` of gene role `g` belongs to a pinned gene.
    pub(super) fn pinned(&self, g: usize, i: usize) -> bool {
        let spans = self.m.genes(g);
        self.pins
            .iter()
            .any(|p| p.role == g && spans[p.gene].contains(&i))
    }

    /// Focus pinned gene `k`, at its member `site` when given, else its first.
    pub(super) fn focus_pin(&mut self, k: usize, site: Option<usize>) {
        let Some(p) = self.pins.get(k) else { return };
        let span = self.m.genes(p.role)[p.gene].clone();
        self.pin = k;
        self.pin_site = site.filter(|s| span.contains(s)).unwrap_or(span.start);
    }

    /// Pin the gene of the module's member `i` at `at` (the end when
    /// past it), or focus it when it is pinned already.
    pub(super) fn pin_site_gene(&mut self, i: usize, at: usize) {
        let Some(p) = self.pin_of(i) else { return };
        let k = match self.pins.iter().position(|&q| q == p) {
            Some(k) => k,
            None => {
                let at = at.min(self.pins.len());
                self.pins.insert(at, p);
                at
            }
        };
        self.focus_pin(k, Some(i));
    }

    pub(super) fn unpin(&mut self, k: usize) {
        if k >= self.pins.len() {
            return;
        }
        self.pins.remove(k);
        let next = if k < self.pin { self.pin - 1 } else { self.pin };
        self.focus_pin(next.min(self.pins.len().saturating_sub(1)), None);
    }

    /// Move pinned gene `k` to place `to` among the pins.
    pub(super) fn move_pin(&mut self, k: usize, to: usize) {
        if k >= self.pins.len() || k == to {
            return;
        }
        let p = self.pins.remove(k);
        let to = to.min(self.pins.len());
        self.pins.insert(to, p);
        self.focus_pin(to, Some(self.pin_site));
    }

    /// The strip's site `by` sites along from the one in focus, crossing
    /// into the next or previous pinned gene at a gene's ends.
    pub(super) fn step_strip(&mut self, by: isize) {
        let Some(p) = self.pins.get(self.pin) else {
            return;
        };
        let span = self.m.genes(p.role)[p.gene].clone();
        let to = self.pin_site as isize + by;
        if to < span.start as isize && self.pin > 0 {
            let prev = self.pins[self.pin - 1];
            let end = self.m.genes(prev.role)[prev.gene].end - 1;
            self.focus_pin(self.pin - 1, Some(end));
        } else if to >= span.end as isize && self.pin + 1 < self.pins.len() {
            self.focus_pin(self.pin + 1, None);
        } else {
            let to = (to.max(span.start as isize) as usize).min(span.end - 1);
            self.focus_pin(self.pin, Some(to));
        }
    }
}
