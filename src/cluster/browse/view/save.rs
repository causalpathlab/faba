//! `s`: the screen saved as a figure, drawn by the figure the command
//! line's `--screen` draws (see [`crate::cluster::browse::figure::screen`])
//! from the view's own state, its size the terminal's; the state saved
//! beside it as JSON, for `--state` to draw it again.

use super::{Look, ModulesView};
use crate::cluster::browse::figure::screen::{self, Beside, Preset, Region, State};
use crate::cluster::browse::figure::GridArgs;
use crate::figure::style::Sheet;

/// The least size a saved screen takes, in points.
const LEAST: (f64, f64) = (432.0, 288.0);

impl ModulesView<'_> {
    /// The state the screen shows now: the grid's order and look (the
    /// blocks' bars, else the diagonal Hinton), the cell modules in the
    /// panes' order, the gene module in focus, the GO column or the sites
    /// beside the grid, and the pinned genes.
    pub fn screen_state(&self) -> State {
        let m = self.m;
        let gene = |p: &super::Pin| m.gene_of(p.role, m.genes(p.role)[p.gene].start).to_string();
        let go = self.go_column && matches!(self.go_table, Some(Ok(_)));
        State {
            preset: match self.looks[0] {
                Look::Bars => Preset::Blocks,
                _ => Preset::Diagonal,
            },
            cell_order: m.columns.iter().map(|c| format!("C{}", c + 1)).collect(),
            module: self.role().map(|g| format!("M{}", g + 1)),
            side: if go { Beside::Go } else { Beside::Sites },
            go_terms: self.go_terms,
            pins: self.pins.iter().map(gene).collect(),
            genes_yscale: self.genes_yscale,
        }
    }

    /// The screen as a figure, the size of the terminal it was last drawn
    /// on (no smaller than [`LEAST`]).
    pub fn screen_svg(&self) -> Result<String, String> {
        self.screen_on(self.screen_sheet())
    }

    /// The sheet a saved screen takes: the terminal's size, no smaller
    /// than [`LEAST`], in the figures' own look.
    pub(super) fn screen_sheet(&self) -> Sheet {
        let page = self.controls.page(self.area, LEAST).unwrap_or(LEAST);
        Sheet {
            page: (Some(page.0), Some(page.1)),
            ..Sheet::OWN
        }
    }

    /// The screen as a figure on `sheet`.
    pub(super) fn screen_on(&self, sheet: Sheet) -> Result<String, String> {
        let state = self.screen_state();
        let placed = state.place(self.m).map_err(|e| format!("{e:#}"))?;
        let go = match &self.go_table {
            Some(Ok(t)) => Some(t),
            _ => None,
        };
        let ask = screen::Ask {
            go,
            models: &self.models,
            grid: &GridArgs::default(),
            region: Region::All,
            genes: Default::default(),
            title: self.figure_title.as_deref().unwrap_or(&self.title),
            sheet,
        };
        screen::draw(self.m, &placed, &ask)
    }

    /// Save the screen under `prefix`, and its state beside it as
    /// `{prefix}.state.json`.
    pub(super) fn save_screen(&mut self, prefix: &str) {
        match self.screen_svg() {
            Ok(svg) => self.controls.save(&svg, prefix),
            Err(e) => return self.status = Some(format!("not saved: {e}")),
        }
        let stem = prefix.strip_suffix(".pdf").unwrap_or(prefix);
        let json = serde_json::to_string_pretty(&self.screen_state()).unwrap_or_default();
        if let Err(e) = std::fs::write(format!("{stem}.state.json"), json) {
            self.status = Some(format!("the state not saved: {e}"));
        }
    }
}
