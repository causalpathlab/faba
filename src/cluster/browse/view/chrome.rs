//! What frames the grids and the strip: the panes' titles, the levels, the
//! site in focus's numbers and the footer's keys.

use super::{Look, ModulesView, Pane, Spot};
use crate::site_analysis::pileup::query_symbol;
use crate::tui::hidden_marks;
use data_beans::interactive::ui::{compact, help_line, panel, DIM, HIGHLIGHT};
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
        let from = self.col_from.min(self.m.columns.len());
        let mut used: usize = spans.iter().map(|s| s.width()).sum();
        let cols = self.m.columns.len();
        for (k, (&c, r)) in self.m.columns[from..]
            .iter()
            .zip(&member.reads[from..])
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
        if let Some(msg) = &self.status {
            return vec![crate::figure::status_line(msg)];
        }
        let grid = |look: Look| {
            let mut keys = vec![("v", "bars/%/z")];
            if look == Look::Percent {
                keys.push(("y", "% scale"));
            }
            keys
        };
        let mut keys: Vec<(&str, &str)> = match self.pane {
            Pane::Modules => [vec![("↑/↓", "module")], grid(self.looks[0])].concat(),
            Pane::Sites => [
                vec![
                    ("↑/↓", "site"),
                    ("Space", "pin gene"),
                    ("Del", "unpin"),
                    ("Alt-drag", "pin onto genes"),
                    ("PgUp/PgDn", "page"),
                    ("o/O", "sort/reverse"),
                ],
                grid(self.looks[1]),
            ]
            .concat(),
            Pane::Genes => vec![
                ("←/→ n/N", "site"),
                ("↑/↓", "pinned gene"),
                ("Del", "unpin"),
                ("X", "clear"),
                ("Alt-drag", "move, off to unpin"),
            ],
        };
        keys.push((
            "Enter",
            match self.pane {
                Pane::Modules => "its sites",
                Pane::Sites => "the pinned genes",
                Pane::Genes => "open the gene",
            },
        ));
        if self.pane != Pane::Modules {
            keys.push(("Esc", "back"));
        }
        if self.col_fit < self.m.cell_roles {
            keys.push(("[/]", "cell modules"));
        }
        keys.extend([("Tab", "pane"), ("/", "find gene")]);
        if self.report.is_some() {
            keys.push(("e", "export report"));
        }
        keys.push(("q", "quit"));
        crate::tui::wrap_help(help_line(&keys), width)
    }

    /// How grid `k` draws, for its title.
    fn look_name(&self, k: usize) -> String {
        match self.looks[k] {
            Look::Bars => "bars".into(),
            Look::Percent => format!("% heatmap, {}", self.ramp.name()),
            Look::Z => "z heatmap".into(),
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
