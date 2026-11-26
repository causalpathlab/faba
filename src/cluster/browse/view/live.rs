//! The browser as a live view (see [`crate::tui::live`]): its state is
//! the screen's [`State`] (`preset`, `cell_order`, `module`, `side`,
//! `go_terms`, `pins`, `genes_yscale`, as `--state` reads it), and its figure the screen
//! `s` saves, drawn by [`crate::cluster::browse::figure::screen::draw`].

use super::{Look, ModulesView};
use crate::cluster::browse::figure::screen::{Beside, Preset, State};
use crate::figure::size::Figure;
use crate::figure::style::Sheet;
use crate::tui::live::Live;
use anyhow::Context;
use serde_json::Value;

impl Live for ModulesView<'_> {
    fn view_name(&self) -> &'static str {
        "cluster_view"
    }

    fn view_state(&self) -> Value {
        serde_json::to_value(self.screen_state()).unwrap_or_default()
    }

    fn set_view_state(&mut self, state: &Value) -> anyhow::Result<()> {
        let given = state.as_object().context("the state is a JSON object")?;
        // The keys it names over the screen's own; others passed over.
        let mut merged = self.view_state();
        for (k, v) in given {
            if merged.get(k).is_some() {
                merged[k] = v.clone();
            }
        }
        let want: State = serde_json::from_value(merged).context("the state")?;
        let placed = want.place(self.m)?;
        anyhow::ensure!(
            placed.columns == self.m.columns,
            "the view keeps its cell modules in the order it opened with; cell_order is set on `--screen`"
        );
        match want.side {
            Beside::Go => {
                anyhow::ensure!(
                    self.read_go(),
                    "no GO for this view: open it with GO terms to set side go"
                );
                self.go_column = true;
            }
            Beside::Sites => self.go_column = false,
            Beside::None => {
                anyhow::bail!("side none is drawn by `--screen` only: the view shows sites or go")
            }
        }
        self.looks[0] = match want.preset {
            Preset::Blocks => Look::Bars,
            Preset::Diagonal => Look::Hinton,
        };
        self.look_set = true;
        if let Some(g) = placed.module {
            self.select_module(g);
        }
        self.go_terms = want.go_terms;
        self.genes_yscale = want.genes_yscale;
        if self.pins != placed.pins {
            self.pins = placed.pins;
            self.focus_pin(0, None);
        }
        Ok(())
    }

    fn own_sheet(&self) -> (Sheet, Figure) {
        (self.screen_sheet(), Figure::Modules)
    }

    fn saved_figure(&self, sheet: Sheet) -> anyhow::Result<String> {
        self.screen_on(sheet).map_err(anyhow::Error::msg)
    }
}
