//! The `/` search: genes found by symbol or key, the view following the
//! match in focus, and back where it was on Esc.

use super::{ModulesView, Pane};
use crate::cluster::browse::Modules;
use crate::figure::{Edit, LineInput};
use crate::site_analysis::pileup::query_symbol;
use data_beans::interactive::ui::{DIM, HIGHLIGHT};
use ratatui::crossterm::event::{KeyCode, KeyEvent};
use ratatui::text::{Line, Span};

/// A gene the `/` search can find: its symbol and key, lower case, its
/// gene role and first member.
pub(super) struct Entry {
    symbol: String,
    key: String,
    role: usize,
    first: usize,
}

/// The `/` search: what is typed, the genes it matches as `(gene role,
/// first member)`, the one shown, and where the view was before.
pub(super) struct Find {
    input: LineInput,
    hits: Vec<(usize, usize)>,
    at: usize,
    back: (usize, usize, Pane),
}

/// Every gene the search can find, gene role by gene role.
pub(super) fn catalog(m: &Modules) -> Vec<Entry> {
    let mut catalog = Vec::new();
    for g in 0..m.gene_roles {
        for span in m.genes(g) {
            let key = m.gene_of(g, span.start);
            catalog.push(Entry {
                symbol: query_symbol(key).to_lowercase(),
                key: key.to_lowercase(),
                role: g,
                first: span.start,
            });
        }
    }
    catalog
}

impl ModulesView<'_> {
    /// `/`: start a search, the view's place kept for Esc.
    pub(super) fn open_find(&mut self) {
        let mut input = LineInput::new(64);
        input.open("");
        self.find = Some(Find {
            input,
            hits: Vec::new(),
            at: 0,
            back: (self.module, self.site, self.pane),
        });
    }

    /// The `/` search's matches for `typed`: symbols starting with it, then
    /// containing it, then genes whose key contains it; shorter first.
    fn matches(&self, typed: &str) -> Vec<(usize, usize)> {
        let q = typed.trim().to_lowercase();
        if q.is_empty() {
            return Vec::new();
        }
        let mut hits: Vec<(u8, &Entry)> = self
            .catalog
            .iter()
            .filter_map(|e| {
                let rank = if e.symbol.starts_with(&q) {
                    0
                } else if e.symbol.contains(&q) {
                    1
                } else if e.key.contains(&q) {
                    2
                } else {
                    return None;
                };
                Some((rank, e))
            })
            .collect();
        hits.sort_by(|(ra, a), (rb, b)| {
            (ra, a.symbol.len(), &a.symbol).cmp(&(rb, b.symbol.len(), &b.symbol))
        });
        hits.into_iter().map(|(_, e)| (e.role, e.first)).collect()
    }

    /// Show the search's match in focus: its module, at its first site.
    fn jump(&mut self) {
        let Some(&(g, site)) = self.find.as_ref().and_then(|f| f.hits.get(f.at)) else {
            return;
        };
        self.select_module(g);
        self.select_site(site);
        self.pane = Pane::Sites;
    }

    pub(super) fn find_key(&mut self, key: KeyEvent) {
        let Some(find) = self.find.as_mut() else {
            return;
        };
        match key.code {
            KeyCode::Down => find.at = (find.at + 1).min(find.hits.len().saturating_sub(1)),
            KeyCode::Up => find.at = find.at.saturating_sub(1),
            _ => match find.input.handle(key) {
                Edit::Submitted(_) => {
                    self.find = None;
                    return;
                }
                Edit::Cancelled => {
                    let (g, site, pane) = find.back;
                    self.find = None;
                    self.select_module(g);
                    self.select_site(site);
                    self.pane = pane;
                    return;
                }
                Edit::Typing => {
                    let typed = find.input.text().unwrap_or_default().to_string();
                    let hits = self.matches(&typed);
                    let find = self.find.as_mut().expect("searching");
                    (find.hits, find.at) = (hits, 0);
                }
            },
        }
        self.jump();
    }

    pub(super) fn find_line(&self, find: &Find) -> Line<'static> {
        let typed = find.input.text().unwrap_or_default();
        let mut spans = vec![
            Span::styled(" find gene: ", DIM),
            Span::styled(format!("{typed}▏"), HIGHLIGHT),
        ];
        for (k, &(g, site)) in find.hits.iter().enumerate().take(6) {
            let name = query_symbol(self.m.gene_of(g, site));
            let style = if k == find.at { HIGHLIGHT } else { DIM };
            spans.push(Span::styled(format!("  {name} → M{}", g + 1), style));
        }
        if find.hits.len() > 6 {
            let more = format!("  +{} more", find.hits.len() - 6);
            spans.push(Span::styled(more, DIM));
        }
        if find.hits.is_empty() && !typed.trim().is_empty() {
            spans.push(Span::styled("  no gene matches", DIM));
        }
        spans.push(Span::styled("   ↑/↓ match  Enter keep  Esc back", DIM));
        Line::from(spans)
    }
}
