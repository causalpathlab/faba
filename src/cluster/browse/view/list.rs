//! The sites list: the module in focus's members, sorted as `o` asks, and
//! the site in focus among them.

use super::{ModulesView, Sort};

impl ModulesView<'_> {
    /// The module's members in the sites list's order.
    #[cfg(test)]
    pub(in crate::cluster::browse) fn listing(&self) -> &[usize] {
        &self.list
    }

    /// Sort the module's members as the list asks.
    pub(super) fn relist(&mut self) {
        let Some(g) = self.role() else {
            self.list.clear();
            return;
        };
        let members = &self.m.members[g];
        let stats = self.m.stats.as_ref();
        // Each member's key; none last whichever way.
        let key = |i: usize| -> Option<f64> {
            let s = members[i].site;
            match self.sort {
                Sort::Module => Some(i as f64),
                Sort::Column(k) => members[i].reads[k].fraction(),
                Sort::Reads => Some(stats.map_or_else(
                    || members[i].reads.iter().map(|r| r.total).sum(),
                    |t| (t.on[s] + t.off[s]) as f64,
                )),
                Sort::Cells => stats.map(|t| t.cells[s] as f64),
            }
        };
        let keys: Vec<Option<f64>> = (0..members.len()).map(key).collect();
        let mut list: Vec<usize> = (0..members.len()).collect();
        list.sort_by(|&a, &b| match (keys[a], keys[b]) {
            (Some(x), Some(y)) => {
                let by = if self.descending {
                    y.total_cmp(&x)
                } else {
                    x.total_cmp(&y)
                };
                by.then(a.cmp(&b))
            }
            (Some(_), None) => std::cmp::Ordering::Less,
            (None, Some(_)) => std::cmp::Ordering::Greater,
            (None, None) => a.cmp(&b),
        });
        self.list = list;
    }

    pub(super) fn sort_name(&self) -> String {
        let name = match self.sort {
            Sort::Module => return "module order".into(),
            Sort::Column(k) => format!("C{} %", self.m.columns[k] + 1),
            Sort::Reads => "reads".into(),
            Sort::Cells => "cells".into(),
        };
        format!("{name} {}", if self.descending { "↓" } else { "↑" })
    }

    /// Sort the list anew, from its top, the site in focus kept.
    pub(super) fn resort(&mut self, sort: Sort, descending: bool) {
        (self.sort, self.descending, self.from) = (sort, descending, 0);
        self.relist();
        self.select_site(self.site);
    }

    /// The next sort: the module's order, each column, reads, cells.
    pub(super) fn next_sort(&mut self) {
        let sort = match self.sort {
            Sort::Module => Sort::Column(0),
            Sort::Column(k) if k + 1 < self.m.cell_roles => Sort::Column(k + 1),
            Sort::Column(_) => Sort::Reads,
            Sort::Reads => Sort::Cells,
            Sort::Cells => Sort::Module,
        };
        self.resort(sort, sort != Sort::Module);
        if let Sort::Column(k) = sort {
            self.show_column(k);
        }
    }

    pub(super) fn select_module(&mut self, g: usize) {
        if g != self.module && self.m.order.contains(&g) {
            (self.module, self.site, self.from) = (g, 0, 0);
            self.relist();
        }
    }

    /// The site in focus's place in the list.
    pub(super) fn listed_at(&self) -> Option<usize> {
        self.list.iter().position(|&i| i == self.site)
    }

    /// Focus member `site`, the list scrolled to keep it in view.
    pub(super) fn select_site(&mut self, site: usize) {
        self.site = site.min(self.sites().saturating_sub(1));
        let at = self.listed_at().unwrap_or(0);
        if at < self.from {
            self.from = at;
        } else if at >= self.from + self.rows {
            self.from = at + 1 - self.rows;
        }
    }

    /// Focus the site `by` places down the list from the one in focus.
    pub(super) fn step_listed(&mut self, by: isize) {
        let Some(at) = self.listed_at() else {
            return;
        };
        self.select_site(self.list[crate::tui::moved(at, by, self.list.len() - 1)]);
    }
}
