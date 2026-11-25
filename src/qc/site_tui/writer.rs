//! Writing the fileset: the writer the run hands the view, and the
//! figures saved with it.

use super::*;
use crate::figure::style::Sheet;

/// Writes the fileset once the thresholds are confirmed: `start` spawns the
/// writer with them, and `progress` is how far it has got.
pub struct Writer<'a> {
    /// The output directory, for the confirmation and the progress.
    pub output: &'a str,
    pub progress: &'a Progress,
    pub start: Option<StartWrite<'a>>,
    /// The page and look the figures are laid out for (see
    /// [`crate::figure::SizeArgs`] and [`crate::figure::style::StyleArgs`]).
    pub sheet: Sheet,
    /// The site panels drawn with the view's figures, as `faba site-plot`
    /// draws them (see [`crate::qc::site_plot::Panels`]); none when `None`.
    pub panels: Option<crate::qc::site_plot::Plan>,
    /// Each batch's cell QC, for the cells figures and the report.
    pub cells: CellsShared,
}

/// Starts the write with the applied thresholds and the view's figures.
pub type StartWrite<'a> = Box<dyn FnOnce(SiteFilterArgs, Vec<Figure>) + 'a>;

/// A figure of the view to save with the fileset: a file stem and its SVG.
pub struct Figure {
    pub stem: String,
    pub svg: String,
}

/// Progress that is already over, for a writer with nothing to write.
static NOTHING_TO_WRITE: Progress = Progress::finished();

impl Writer<'_> {
    /// A writer that writes nothing: applying only decides.
    pub(super) fn idle() -> Self {
        Writer {
            output: "",
            progress: &NOTHING_TO_WRITE,
            start: None,
            sheet: Sheet::default(),
            panels: None,
            cells: CellsShared::default(),
        }
    }
}

impl SitePicker<'_> {
    /// Refused when `sheet`'s page is short of what a knob's figure needs
    /// with the boxes the view opens with: the gene's when a table has
    /// genes, the metagene's when an annotation is named (see
    /// [`knob_room`]). Checked before the view opens, as the width is.
    pub(super) fn page_problem(&self, sheet: &Sheet) -> anyhow::Result<()> {
        let genes = self.views.iter().any(|v| v.genes.is_some());
        let boxes = usize::from(genes) + usize::from(self.annotation.gff.is_some());
        knob_room(sheet, knob_least(sheet, Some(boxes)), "a qc figure")
    }

    /// Apply the confirmed thresholds: start the write, which [`Screen::tick`]
    /// then follows to the end.
    pub(super) fn apply(&mut self) {
        self.mode = Mode::Writing;
    }

    /// Draw the figures and start the writer, once the writing pop-up is up.
    /// A figure the page is still too small for (an annotation read in the
    /// view adds the metagene) is said, and the view stays open, unwritten.
    pub(super) fn start_writing(&mut self) {
        let Some(start) = self.writer.start.take() else {
            return;
        };
        let plan = self.writer.panels.clone();
        match self.written_figures(plan.as_ref(), self.writer.sheet) {
            Ok(figures) => start(self.filter.clone(), figures),
            Err(e) => {
                self.writer.start = Some(start);
                self.mode = Mode::Browse;
                self.status = Some(format!("not written: {e}"));
            }
        }
    }

    /// The figures saved with the fileset, laid out for `sheet`: the
    /// view's (see [`Self::figures`]), then the site panels `plan` asks
    /// for, if any (see [`Self::site_panels`]).
    pub(super) fn written_figures(
        &mut self,
        plan: Option<&crate::qc::site_plot::Plan>,
        sheet: Sheet,
    ) -> anyhow::Result<Vec<Figure>> {
        let mut out = self.figures(sheet)?;
        if let Some(plan) = plan {
            out.extend(self.site_panels(plan, sheet)?);
        }
        Ok(out)
    }

    /// The site panels `plan` asks for of every modality at the thresholds
    /// now, as `faba site-plot` draws them, under `panels/`. A quantity
    /// asked for that a modality lacks is left out of that modality's.
    pub(super) fn site_panels(
        &self,
        plan: &crate::qc::site_plot::Plan,
        sheet: Sheet,
    ) -> anyhow::Result<Vec<Figure>> {
        let mut out = Vec::new();
        for v in &self.views {
            let panels =
                crate::qc::site_plot::Panels::of_view(&v.sites(), &self.filter, plan, sheet)?;
            out.extend(panels.figures().into_iter().map(|f| Figure {
                stem: format!("panels/{}", f.stem),
                svg: f.svg,
            }));
        }
        Ok(out)
    }

    /// The view's figure for every modality and every knob, under the
    /// applied thresholds and the scales chosen per knob, and each
    /// modality's pairs (not alone, `--panel-only`: the site panels are
    /// those), each laid out for `sheet` (see [`crate::figure::SizeArgs`]);
    /// refused when one is too small. The modality on screen keeps its selected gene; the others
    /// show their top gene, the gene filter aside. The view is left as it
    /// was.
    pub(super) fn figures(&mut self, sheet: Sheet) -> anyhow::Result<Vec<Figure>> {
        let (modality, focus, at) = (self.modality, self.focus, self.list.at);
        let find = std::mem::take(&mut self.list.find);
        let out = self.each_figure(sheet, (modality, at), &find);
        self.restore((modality, focus, at), find);
        let mut out = out?;
        out.extend(cell_figures(&self.cells, sheet));
        Ok(out)
    }

    /// [`Self::figures`]' drawings, the view stepped through each modality
    /// and knob (`modality` the one on screen, its gene at `at`, found by
    /// `find`).
    fn each_figure(
        &mut self,
        sheet: Sheet,
        (modality, at): (usize, usize),
        find: &str,
    ) -> anyhow::Result<Vec<Figure>> {
        let mut out = Vec::new();
        for m in 0..self.views.len() {
            self.modality = m;
            if m == modality {
                find.clone_into(&mut self.list.find);
            }
            self.refresh_genes(true);
            self.list.find.clear();
            if m == modality {
                self.list.at = at;
            }
            // Alone, each knob's histogram only, which needs no column of
            // the view's own; and the pairs are the site panels, a figure
            // each.
            if sheet.look.panel_only {
                out.extend(self.knob_panels_of(m, sheet)?);
                continue;
            }
            for f in 0..self.views[m].criteria.len() {
                self.focus = f;
                self.rebuild_column();
                out.push(Figure {
                    stem: super::draw::knob_stem(self.view().table, self.criterion()),
                    svg: self.figure(sheet)?,
                });
            }
            out.push(Figure {
                stem: super::pairs_view::pairs_stem(&self.view().table.modality),
                svg: self.pairs_figure(m, sheet),
            });
        }
        Ok(out)
    }

    /// The modality, knob and gene the view had before
    /// [`Self::figures`] stepped through them, put back.
    fn restore(&mut self, (modality, focus, at): (usize, usize, usize), find: String) {
        self.modality = modality;
        self.list.find = find;
        self.refresh_genes(true);
        self.list.at = at;
        self.focus = focus;
        self.rebuild_column();
    }
}
