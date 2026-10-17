//! Writing the fileset: the writer the run hands the view, and the
//! figures saved with it.

use super::*;

/// Writes the fileset once the thresholds are confirmed: `start` spawns the
/// writer with them, and `progress` is how far it has got.
pub struct Writer<'a> {
    /// The output directory, for the confirmation and the progress.
    pub output: &'a str,
    pub progress: &'a Progress,
    pub start: Option<StartWrite<'a>>,
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
        }
    }
}

impl SitePicker<'_> {
    /// Apply the confirmed thresholds: start the write, which [`Screen::tick`]
    /// then follows to the end.
    pub(super) fn apply(&mut self) {
        self.mode = Mode::Writing;
    }

    /// Draw the figures and start the writer, once the writing pop-up is up.
    pub(super) fn start_writing(&mut self) {
        if let Some(start) = self.writer.start.take() {
            let figures = self.figures();
            start(self.filter.clone(), figures);
        }
    }

    /// The view's figure for every modality and every knob, under the
    /// applied thresholds and the scales chosen per knob. The modality on
    /// screen keeps its selected gene; the others show their top gene, the
    /// gene filter aside.
    pub(super) fn figures(&mut self) -> Vec<Figure> {
        let (modality, focus, at) = (self.modality, self.focus, self.list.at);
        let find = std::mem::take(&mut self.list.find);
        let mut out = Vec::new();
        for m in 0..self.views.len() {
            self.modality = m;
            if m == modality {
                self.list.find.clone_from(&find);
            }
            self.refresh_genes(true);
            self.list.find.clear();
            if m == modality {
                self.list.at = at;
            }
            for f in 0..self.views[m].criteria.len() {
                self.focus = f;
                self.rebuild_column();
                let knob = self.criterion().flag().trim_start_matches("--site-");
                out.push(Figure {
                    stem: format!("{}_{knob}", self.view().table.modality),
                    svg: self.figure(),
                });
            }
        }
        self.modality = modality;
        self.list.find = find;
        self.refresh_genes(true);
        self.list.at = at;
        self.focus = focus;
        self.rebuild_column();
        out
    }
}
