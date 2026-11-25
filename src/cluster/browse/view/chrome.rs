//! What frames the grids and the strip: the panes' titles, the levels, the
//! site in focus's numbers and the footer's keys.

use super::{Full, Look, ModulesView, Pane, Sort, Spot};
use crate::site_analysis::pileup::query_symbol;
use crate::tui::hidden_marks;
use crate::tui::ui::{compact, help_line, panel, DIM, HIGHLIGHT};
use ratatui::layout::Rect;
use ratatui::text::{Line, Span};
use ratatui::Frame;

impl ModulesView<'_> {
    /// The site in focus: its gene, place and module, then each cell
    /// module's % and reads from the first shown, as many as fit in `width`;
    /// those left out either side counted.
    pub(super) fn stats_line(&self, width: u16) -> Line<'static> {
        let Some((g, i)) = self.focus() else {
            return Line::from(Span::styled(" no module has sites", DIM));
        };
        let Some(member) = self.m.members[g].get(i) else {
            return Line::default();
        };
        let s = &self.m.sites[member.site];
        let mut spans = vec![
            Span::styled(format!(" {}", query_symbol(&s.gene)), HIGHLIGHT),
            Span::styled(format!("  {}:{}", s.chr, s.pos), DIM),
            Span::styled(format!("   of M{}", g + 1), DIM),
            Span::raw("  "),
        ];
        if self.col_from > 0 {
            spans.push(Span::styled(
                format!("  {}", hidden_marks(self.col_from, 0)),
                DIM,
            ));
        }
        let from = self.col_from.min(self.shown.len());
        let mut used: usize = spans.iter().map(|s| s.width()).sum();
        let cols = self.shown.len();
        let shown = self.shown[from..].iter();
        for (k, (c, r)) in shown
            .map(|&p| (self.m.columns[p], &member.reads[p]))
            .enumerate()
        {
            let name = format!("  C{} ", c + 1);
            let reads = format!("{} of {}", r.pct(), compact(r.total));
            // Room kept for counting those after it, unless it is the last.
            let after = cols - from - k - 1;
            let mark = if after > 0 {
                hidden_marks(0, after).chars().count() + 2
            } else {
                0
            };
            let len = name.chars().count() + reads.chars().count();
            if used + len + mark > width as usize {
                spans.push(Span::styled(
                    format!("  {}", hidden_marks(0, after + 1)),
                    DIM,
                ));
                break;
            }
            used += len;
            spans.push(Span::styled(name, DIM));
            spans.push(Span::raw(reads));
        }
        Line::from(spans)
    }

    /// The footer: the search being typed, a status, or the keys that act
    /// where the focus is.
    pub(super) fn footer(&self, width: u16) -> Vec<Line<'static>> {
        if let Some(find) = &self.find {
            return vec![self.find_line(find)];
        }
        match (&self.status, self.footer_pairs()) {
            (Some(msg), _) => vec![crate::figure::status_line(msg)],
            (None, keys) => crate::tui::wrap_help(help_line(&keys.unwrap_or_default()), width),
        }
    }

    /// The footer's keys where the focus is, a full pane's own or the
    /// panes': none while a gene is typed or (but over the trace and GO)
    /// a status shows.
    pub(super) fn footer_pairs(&self) -> Option<Vec<(&'static str, &'static str)>> {
        let close = [("q", "quit"), ("?", "keys")];
        // The trace and GO panes say no status.
        match self.full {
            _ if self.find.is_some() => return None,
            Some(Full::Trace) => return Some([&[("t/Esc", "close")][..], &close].concat()),
            Some(Full::Go) => {
                let go = [("↑/↓", "gene module"), ("G/Esc", "close")];
                return Some([&go[..], &close].concat());
            }
            _ if self.status.is_some() => return None,
            Some(Full::Batches) => {
                let batches = [("b/Esc", "close"), ("B", "filter"), ("↑/↓", "scroll")];
                return Some([&batches[..], &close].concat());
            }
            None => {}
        }
        let grid = |look: Look| {
            let mut keys = vec![("v", "bars/%/z/Hinton")];
            if look == Look::Percent {
                keys.push(("y", "% scale"));
            }
            keys
        };
        let mut keys: Vec<(&str, &str)> = match self.pane {
            Pane::Modules => [vec![("↑/↓", "module")], grid(self.looks[0])].concat(),
            Pane::Sites => {
                let mut keys = vec![("↑/↓", "site"), ("Space", "pin gene")];
                let pinned = self.pin_of(self.site).filter(|p| self.pins.contains(p));
                if pinned.is_some() {
                    keys.push(("Del", "unpin"));
                }
                keys.extend([
                    ("Ctrl-drag", "pin onto genes"),
                    ("PgUp/PgDn", "page"),
                    if self.sort == Sort::Module {
                        ("o", "sort")
                    } else {
                        ("o/O", "sort/reverse")
                    },
                ]);
                [keys, grid(self.looks[1])].concat()
            }
            Pane::Genes => {
                let mut keys = Vec::new();
                if !self.pins.is_empty() {
                    keys.push(("←/→ n/N", "site"));
                }
                if self.pins.len() > 1 {
                    keys.push(("↑/↓", "pinned gene"));
                }
                if !self.pins.is_empty() {
                    keys.extend([("Del", "unpin"), ("Ctrl-drag", "move, off to unpin")]);
                }
                keys
            }
        };
        if !self.pins.is_empty() {
            keys.push(("X", "clear pins"));
        }
        // Into the genes pane with none pinned, there is no gene to open.
        match self.pane {
            Pane::Modules => keys.push(("Enter", "its sites")),
            Pane::Sites => keys.push(("Enter", "the pinned genes")),
            Pane::Genes if !self.pins.is_empty() => keys.push(("Enter", "open the gene")),
            Pane::Genes => {}
        }
        if self.pane != Pane::Modules {
            keys.push(("Esc", "back"));
        }
        if self.col_fit < self.shown.len() {
            keys.push(("[/]", "cell modules"));
        }
        keys.extend([("Tab", "pane"), ("/", "find gene")]);
        if self.report.is_some() {
            keys.push((crate::tui::REPORT_KEY, "export report"));
        }
        if self.m.batches.is_some() {
            keys.extend([("b", "batches"), ("B", "filter")]);
        }
        if self.go_ask.is_some() {
            keys.extend([("a", "GO column"), ("G", "GO terms")]);
        }
        keys.extend([("t", "trace")]);
        keys.extend(close);
        Some(keys)
    }

    /// Every key where the focus is, for the `?` pop-up: a full pane's
    /// own, or the panes'.
    pub(super) fn help_keys(&self) -> Vec<(&'static str, &'static str)> {
        let close = [("Esc", "close"), ("q", "quit"), ("?", "these keys")];
        match self.full {
            Some(Full::Trace) => return [&[("t", "close the trace")][..], &close].concat(),
            Some(Full::Batches) => {
                return [
                    &[
                        ("↑/↓ k/j", "scroll the cell modules"),
                        ("B", "batch filter"),
                        ("b", "close the batches"),
                    ][..],
                    &close,
                ]
                .concat()
            }
            Some(Full::Go) => {
                return [&[("↑/↓ k/j", "gene module"), ("G", "close GO")][..], &close].concat()
            }
            None => {}
        }
        let mut keys = vec![("Tab/Shift+Tab", "next/previous pane")];
        match self.pane {
            Pane::Modules => keys.extend([("↑/↓ k/j", "gene module"), ("Enter", "into its sites")]),
            Pane::Sites => {
                keys.extend([
                    ("↑/↓ k/j", "site"),
                    ("PgUp/PgDn", "a page of sites"),
                    ("Home/End", "first/last site"),
                    ("Space", "pin the site's gene"),
                ]);
                let pinned = self.pin_of(self.site).filter(|p| self.pins.contains(p));
                if pinned.is_some() {
                    keys.push(("Del/Backspace", "unpin it"));
                }
                keys.push(("o", "sort by the next column, reads, cells"));
                // In the module's order there is no other way round.
                if self.sort != Sort::Module {
                    keys.push(("O", "reverse the sort"));
                }
                keys.extend([("Enter", "into the pinned genes"), ("Esc", "back")]);
            }
            Pane::Genes => {
                if self.pins.len() > 1 {
                    keys.push(("↑/↓ k/j", "pinned gene"));
                }
                if !self.pins.is_empty() {
                    keys.extend([
                        ("←/→ h/l n/N", "previous/next site"),
                        ("Del/Backspace", "unpin the gene"),
                        ("Enter", "open the gene"),
                    ]);
                }
                keys.push(("Esc", "back"));
            }
        }
        if !self.pins.is_empty() {
            keys.push(("X", "clear every pin (asks)"));
        }
        // The genes pane has no grid to draw another way.
        if self.pane != Pane::Genes {
            keys.push(("v", "bars, %, z or Hinton (squares sized by |z|)"));
            if self.looks[self.pane.index()] == Look::Percent {
                keys.push(("y", "% scale"));
            }
        }
        if self.col_fit < self.shown.len() {
            keys.push(("[ ]", "scroll the cell modules"));
        }
        keys.extend([("/", "find a gene"), ("t", "trace")]);
        if self.report.is_some() {
            keys.push((crate::tui::REPORT_KEY, "export the report"));
        }
        // Without them, the keys say why there are none.
        keys.extend(match self.m.batches {
            Some(_) => &[("b", "batches"), ("B", "batch filter")][..],
            None => &[("b/B", "batches: the cells are of one batch")],
        });
        keys.extend(match self.go_ask {
            Some(_) => &[
                (
                    "a",
                    "GO column beside the gene modules: read GO, then show or hide it",
                ),
                ("G", "every GO term of the gene module in focus"),
            ][..],
            None => &[("a/G", "GO: none for this view")],
        });
        keys.extend([("q", "quit"), ("?", "these keys")]);
        keys
    }

    /// How grid `k` draws, for its title.
    fn look_name(&self, k: usize) -> String {
        match self.looks[k] {
            Look::Bars => "bars".into(),
            Look::Percent => format!("% heatmap, {}", self.ramp.name()),
            Look::Z => "z heatmap".into(),
            Look::Hinton => "Hinton, z".into(),
        }
    }

    /// The levels, the one in focus highlighted.
    pub(super) fn levels(&self) -> Line<'static> {
        let at = self.pane.index();
        let mut spans = vec![Span::raw(" ")];
        for (k, name) in ["modules", "sites", "genes", "pileup"].iter().enumerate() {
            if k > 0 {
                spans.push(Span::styled(" › ", DIM));
            }
            let style = if k == at { HIGHLIGHT } else { DIM };
            spans.push(Span::styled(format!("{} {name}", k + 1), style));
        }
        Line::from(spans)
    }

    /// The three panes' frames, titled, the one in focus highlighted, each
    /// a place to click; the areas inside them.
    pub(super) fn panels(
        &mut self,
        frame: &mut Frame,
        [modules, sites, genes]: [Rect; 3],
    ) -> [Rect; 3] {
        let module_name = self.role().map_or(String::new(), |g| format!("M{}", g + 1));
        let at = self.listed_at().map_or(0, |a| a + 1);
        let order = match self.looks[0] {
            Look::Bars => "by shift range",
            _ => "diagonal",
        };
        let pinned = match self.pins.len() {
            0 => "none pinned".to_string(),
            n => format!("{n} pinned"),
        };
        let blocks = [
            (
                modules,
                Pane::Modules,
                format!(" gene modules · {order} · {} ", self.look_name(0)),
            ),
            (
                sites,
                Pane::Sites,
                format!(
                    " sites · {module_name} · {at}/{} · {} · {} ",
                    self.sites(),
                    self.sort_name(),
                    self.look_name(1)
                ),
            ),
            (
                genes,
                Pane::Genes,
                format!(
                    " genes · {pinned} · posterior-mean reads in the main module, fitted cells "
                ),
            ),
        ];
        let mut inner = [Rect::default(); 3];
        for (k, (area, pane, title)) in blocks.into_iter().enumerate() {
            let block = panel(title, self.pane == pane);
            inner[k] = block.inner(area);
            self.hits.add(area, Spot::Pane(pane));
            frame.render_widget(block, area);
        }
        inner
    }
}
