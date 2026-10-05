//! The annotation read on a thread, and what the gene and metagene bars add up.

use super::*;

/// Total metagene bins across 5'UTR, CDS and 3'UTR: enough for the shape,
/// and at two columns a bin they fill a panel on a 150-column terminal.
pub(super) const META_BINS: usize = 48;

/// What the gene and metagene bars add up: the read coverage at the sites
/// (the default), their converted reads, or the sites themselves.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Weight {
    Coverage,
    Converted,
    Sites,
}

impl Weight {
    /// The y axis, as the titles name it.
    pub(super) fn unit(self) -> &'static str {
        match self {
            Weight::Coverage => "read coverage",
            Weight::Converted => "converted reads",
            Weight::Sites => "sites",
        }
    }

    /// The one `c` switches to.
    pub(super) fn next(self) -> Self {
        match self {
            Weight::Coverage => Weight::Converted,
            Weight::Converted => Weight::Sites,
            Weight::Sites => Weight::Coverage,
        }
    }
}

/// What the annotation gives the picker: per view, in view order, the
/// metagene layout (`None` where no site is on a coding transcript); and the
/// gene models of the genes the tables name, by `{gene_id}_{symbol}`.
pub(super) struct Annotation {
    pub(super) views: Vec<Option<MetaLayout>>,
    /// Why there is no metagene at all, when the annotation is to blame.
    pub(super) no_metagene: Option<String>,
    pub(super) models: FxHashMap<Box<str>, GeneModel>,
}

/// The annotation's state.
pub(super) enum Meta {
    /// Nothing to draw; the reason, for the panels.
    Unavailable(String),
    /// The annotation being read and the sites placed, on a thread.
    Pending(std::thread::JoinHandle<Result<Annotation, String>>),
    Ready(Annotation),
}

/// What an annotation is read for, kept so another can be read in its place.
#[derive(Default)]
pub(super) struct AnnotationSource {
    /// The annotation the metagene and gene models come from.
    pub(super) gff: Option<Box<str>>,
    /// Where the annotation browser opens when `gff` has no directory.
    pub(super) dir: std::path::PathBuf,
    /// Every table's sites, in view order.
    pub(super) batches: Vec<RecordBatch>,
    /// The genes whose models to keep.
    pub(super) keys: FxHashSet<Box<str>>,
}

/// Whether the annotation browser lists `name`: a GTF or GFF, gzipped or not.
pub(super) fn is_annotation(name: &str) -> bool {
    let name = name.to_ascii_lowercase();
    let name = name.strip_suffix(".gz").unwrap_or(&name);
    [".gtf", ".gff", ".gff3"]
        .iter()
        .any(|ext| name.ends_with(ext))
}

impl Meta {
    /// Read `source`'s annotation once, on a thread: place every table's
    /// sites on the metagene, and keep the models of its genes.
    pub(super) fn start(source: &AnnotationSource) -> Self {
        let Some(gff) = source.gff.as_deref().map(str::to_string) else {
            return Meta::Unavailable(
                "no annotation: pass --gff, or keep the run record next to the sites".into(),
            );
        };
        // A run record names the annotation where the run was, which may be
        // another machine.
        if !std::path::Path::new(&gff).is_file() {
            return Meta::Unavailable(format!("{gff}: not found"));
        }
        let (batches, keys) = (source.batches.clone(), source.keys.clone());
        Meta::Pending(std::thread::spawn(move || {
            let fail = |e: anyhow::Error| format!("{gff}: {e:#}");
            // Gene models use gene and exon lines; the metagene's transcripts
            // exon, CDS and stop codon lines. Nothing else is parsed.
            let features = ["gene", "Gene", "exon", "CDS", "cds", "stop_codon"];
            let records = read_records_of(&gff, &features).map_err(fail)?;
            let meta = MetaModels::from_records(&records);
            let no_metagene = meta.is_empty().then(|| {
                format!("{gff}: no coding transcript (exon and CDS lines) to build a metagene")
            });
            let views = batches
                .iter()
                .map(|b| {
                    let sites = genomic_sites(b).map_err(|e| format!("{e:#}"))?;
                    Ok(meta.layout(&sites, META_BINS))
                })
                .collect::<Result<Vec<_>, String>>()?;
            let models = gene_models_from_records(&records, |k| keys.contains(k))
                .map_err(fail)?
                .into_iter()
                .map(|m| (m.key.clone(), m))
                .collect();
            Ok(Annotation {
                views,
                no_metagene,
                models,
            })
        }))
    }

    #[cfg(test)]
    pub(super) fn ready(layouts: Vec<Option<MetaLayout>>, models: Vec<GeneModel>) -> Self {
        Meta::Ready(Annotation {
            views: layouts,
            no_metagene: None,
            models: models.into_iter().map(|m| (m.key.clone(), m)).collect(),
        })
    }

    /// Take the thread's result if it has arrived; true when it just did.
    pub(super) fn poll(&mut self) -> bool {
        let Meta::Pending(handle) = self else {
            return false;
        };
        if !handle.is_finished() {
            return false;
        }
        let Meta::Pending(handle) = std::mem::replace(self, Meta::Unavailable(String::new()))
        else {
            unreachable!()
        };
        *self = match handle.join() {
            Ok(Ok(annotation)) => Meta::Ready(annotation),
            Ok(Err(e)) => Meta::Unavailable(e),
            Err(_) => Meta::Unavailable("reading the annotation panicked".into()),
        };
        true
    }
}

/// A view's metagene: all sites, the kept ones, and the bins per region.
/// A site on k isoforms adds 1/k to each, so bins are fractional.
pub(super) struct MetaCounts<'m> {
    pub(super) all: &'m [f64],
    pub(super) kept: &'m [f64],
    pub(super) regions: [usize; 3],
    pub(super) unassigned: usize,
}

impl<'a> SitePicker<'a> {
    pub(super) fn with_refreshed_meta(mut self) -> Self {
        self.refresh_meta();
        self
    }

    #[cfg(test)]
    pub(super) fn set_meta(&mut self, meta: Meta) {
        self.meta = meta;
        self.refresh_meta();
    }

    fn view_meta(&self) -> Option<&MetaLayout> {
        let Meta::Ready(a) = &self.meta else {
            return None;
        };
        a.views.get(self.modality)?.as_ref()
    }

    /// Site `i`'s bar weight in the focused view.
    pub(super) fn site_weight(&self, i: usize) -> usize {
        let t = &self.view().table;
        match self.weight {
            Weight::Coverage => t.coverage[i] as usize,
            Weight::Converted => t.converted[i] as usize,
            Weight::Sites => 1,
        }
    }

    /// Recount the focused view's metagene bars; redraw only if they changed.
    pub(super) fn refresh_meta(&mut self) {
        let fails = &self.views[self.modality].fails;
        let bars = self
            .view_meta()
            .map(|m| {
                let all = m.counts(|i| self.site_weight(i) as f64);
                let kept = m.counts(|i| {
                    if fails[i] == 0 {
                        self.site_weight(i) as f64
                    } else {
                        0.0
                    }
                });
                (all, kept)
            })
            .unwrap_or_default();
        if bars != self.meta_bars {
            self.meta_bars = bars;
            self.meta_plot.invalidate();
        }
    }

    pub(super) fn switch_weight(&mut self) {
        self.weight = self.weight.next();
        self.refresh_meta();
    }

    /// The focused view's metagene under the current thresholds, `None`
    /// until the layouts arrive or when no site is on a coding transcript.
    pub(super) fn meta_counts(&self) -> Option<MetaCounts<'_>> {
        let m = self.view_meta()?;
        Some(MetaCounts {
            all: &self.meta_bars.0,
            kept: &self.meta_bars.1,
            regions: m.region_bins(),
            unassigned: m.unassigned,
        })
    }

    /// Open the annotation browser where the current annotation is, else
    /// where the sites are.
    pub(super) fn open_gff_browser(&mut self) {
        let dir = self
            .annotation
            .gff
            .as_deref()
            .and_then(|g| std::path::Path::new(g).parent())
            .filter(|d| d.is_dir())
            .map_or_else(|| self.annotation.dir.clone(), |d| d.to_path_buf());
        let mut browser = Browser::new(dir.clone());
        with_gff_listing(|l| browser.open(dir, l));
        self.mode = Mode::Gff(browser);
    }

    /// A key in the annotation browser: Enter on a file reads it.
    pub(super) fn gff_key(&mut self, key: KeyEvent) {
        let Mode::Gff(browser) = &mut self.mode else {
            return;
        };
        match with_gff_listing(|l| browser.key(key, l)) {
            Nav::Moved => {}
            Nav::Picked(path) => {
                self.mode = Mode::Browse;
                self.load_gff(path.to_string_lossy().into());
            }
            Nav::Ignored => {
                if matches!(key.code, KeyCode::Esc | KeyCode::Char('q')) {
                    self.mode = Mode::Browse;
                }
            }
        }
    }

    /// Read `gff` in place of the current annotation.
    pub(super) fn load_gff(&mut self, gff: Box<str>) {
        self.annotation.gff = Some(gff);
        self.meta = Meta::start(&self.annotation);
        self.refresh_meta();
    }

    /// The annotation the view uses now, for the run record.
    pub(super) fn gff(&self) -> Option<&str> {
        self.annotation.gff.as_deref()
    }

    /// Why there is no metagene to draw.
    pub(super) fn meta_status(&self) -> String {
        match &self.meta {
            Meta::Unavailable(why) => why.clone(),
            Meta::Pending(_) => "reading gene models ...".into(),
            Meta::Ready(a) => a
                .no_metagene
                .clone()
                .unwrap_or_else(|| "no site on a coding transcript".into()),
        }
    }
}

/// `f` with the annotation browser's listing: directories and annotation
/// files, none tagged.
fn with_gff_listing<R>(f: impl FnOnce(&mut Listing) -> R) -> R {
    f(&mut Listing {
        keep: &is_annotation,
        tag: &mut |_| false,
    })
}
