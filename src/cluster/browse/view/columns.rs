//! The cell modules the grids and the strip show: how wide each is, how
//! many fit, and which one comes first, the same in every pane.

use super::grid::{CELL, COMPACT, COUNTS, MODULE_EXTRA, MODULE_LABEL, SITE_LABEL};
use super::ModulesView;
use std::ops::Range;

impl ModulesView<'_> {
    /// How many cell modules fit in `room` columns from column `from`.
    fn visible_from(&self, from: usize, room: u16) -> usize {
        let fit = (room / self.cell_w).max(1) as usize;
        fit.min(self.m.cell_roles - from.min(self.m.cell_roles))
    }

    /// The cell modules a grid with `room` columns for its cells shows:
    /// as many as fit from the first shown, at least one.
    pub(super) fn visible(&self, room: u16) -> Range<usize> {
        let from = self.col_from.min(self.m.cell_roles.saturating_sub(1));
        from..from + self.visible_from(from, room)
    }

    /// Shift the cell modules shown `by` places, no further than the
    /// last ones filling the narrowest grid.
    pub(super) fn scroll_columns(&mut self, by: isize) {
        let last = self.m.cell_roles.saturating_sub(self.col_fit);
        self.col_from = crate::tui::moved(self.col_from, by, last);
    }

    /// Scroll so that column `k` is shown.
    pub(super) fn show_column(&mut self, k: usize) {
        if k < self.col_from {
            self.col_from = k;
        } else if k >= self.col_from + self.col_fit {
            self.col_from = k + 1 - self.col_fit;
        }
    }

    /// The widths the modules and sites panes need, borders included, with
    /// every cell module `cell_w` columns wide.
    fn needs(&self, cell_w: u16) -> (u16, u16) {
        let cells = self.m.cell_roles as u16 * cell_w;
        (
            MODULE_LABEL + cells + MODULE_EXTRA + 2,
            SITE_LABEL + cells + COUNTS + 4,
        )
    }

    /// The widths of the modules and sites panes: what each needs, or when
    /// both do not fit, `width` split in proportion to it, the modules
    /// pane kept wide enough for one cell module when there is room.
    pub(super) fn split(&self, width: u16) -> u16 {
        let (modules, sites) = self.needs(self.cell_w);
        if modules + sites <= width {
            return modules;
        }
        let share = (width as u32 * modules as u32 / (modules + sites) as u32) as u16;
        let least = MODULE_LABEL + self.cell_w + MODULE_EXTRA + 2;
        share.max(least.min(width / 2))
    }

    /// For grids `width` wide together: compact cells when not every cell
    /// module fits at full width, and how many the narrower grid has room
    /// for; every pane starts at the same one.
    pub(super) fn fit_columns(&mut self, width: u16) {
        let (modules, sites) = self.needs(CELL);
        self.cell_w = if modules + sites <= width {
            CELL
        } else {
            COMPACT
        };
        let modules = self.split(width);
        let room = |pane: u16, kept: u16| pane.saturating_sub(2 + kept);
        let fits = [
            room(modules, MODULE_LABEL + MODULE_EXTRA),
            room(width - modules, SITE_LABEL + COUNTS + 1),
        ]
        .map(|r| self.visible_from(0, r));
        self.col_fit = fits.into_iter().min().unwrap_or(1);
        self.scroll_columns(0);
    }
}
