//! The annotation read on a thread, and what the gene and metagene bars add up.

use super::*;

/// Total metagene bins across 5'UTR, CDS and 3'UTR: enough for the shape,
/// and at two columns a bin they fill a panel on a 150-column terminal.
pub(super) const META_BINS: usize = 48;

/// What the gene and metagene bars add up: sites, or their converted reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub(super) enum Weight {
    Sites,
    Converted,
}

impl Weight {
    /// The y axis, as the titles name it.
    pub(super) fn unit(self) -> &'static str {
        match self {
            Weight::Sites => "sites",
            Weight::Converted => "converted reads",
        }
    }

    pub(super) fn other(self) -> Self {
        match self {
            Weight::Sites => Weight::Converted,
            Weight::Converted => Weight::Sites,
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

impl Meta {
    /// Read `gff` once, on a thread: place every table's sites (`batches`,
    /// in view order) on the metagene, and keep the models of the genes in
    /// `keys`.
    pub(super) fn start(
        gff: Option<&str>,
        batches: Vec<RecordBatch>,
        keys: FxHashSet<Box<str>>,
    ) -> Self {
        let Some(gff) = gff.map(str::to_string) else {
            return Meta::Unavailable(
                "no annotation: pass --gff, or keep the run record next to the sites".into(),
            );
        };
        Meta::Pending(std::thread::spawn(move || {
            let fail = |e: anyhow::Error| format!("{gff}: {e:#}");
            let records = read_gff_record_vec(&gff).map_err(fail)?;
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
pub(super) struct MetaCounts<'m> {
    pub(super) all: &'m [usize],
    pub(super) kept: &'m [usize],
    pub(super) regions: [usize; 3],
    pub(super) unassigned: usize,
}

impl<'a> SitePicker<'a> {
    pub(super) fn with_meta(mut self, meta: Meta) -> Self {
        self.set_meta(meta);
        self
    }

    pub(super) fn set_meta(&mut self, meta: Meta) {
        self.meta = meta;
        self.refresh_meta();
    }

    pub(super) fn view_meta(&self) -> Option<&MetaLayout> {
        let Meta::Ready(a) = &self.meta else {
            return None;
        };
        a.views.get(self.modality)?.as_ref()
    }

    /// Site `i`'s bar weight in the focused view.
    pub(super) fn site_weight(&self, i: usize) -> usize {
        match self.weight {
            Weight::Sites => 1,
            Weight::Converted => self.view().table.converted[i] as usize,
        }
    }

    /// Recount the focused view's metagene over all sites, after the view,
    /// the annotation or the weight changed; then the kept sites.
    pub(super) fn refresh_meta(&mut self) {
        self.meta_all = self
            .view_meta()
            .map(|m| m.counts(|i| self.site_weight(i)))
            .unwrap_or_default();
        self.meta_plot.invalidate();
        self.update_meta();
    }

    /// Recount the focused view's kept sites; redraw only if they changed.
    pub(super) fn update_meta(&mut self) {
        let fails = &self.views[self.modality].fails;
        let kept = self
            .view_meta()
            .map(|m| {
                m.counts(|i| {
                    if fails[i] == 0 {
                        self.site_weight(i)
                    } else {
                        0
                    }
                })
            })
            .unwrap_or_default();
        if kept != self.meta_kept {
            self.meta_kept = kept;
            self.meta_plot.invalidate();
        }
    }

    pub(super) fn switch_weight(&mut self) {
        self.weight = self.weight.other();
        self.refresh_meta();
    }

    /// The focused view's metagene under the current thresholds, `None`
    /// until the layouts arrive or when no site is on a coding transcript.
    pub(super) fn meta_counts(&self) -> Option<MetaCounts<'_>> {
        let m = self.view_meta()?;
        Some(MetaCounts {
            all: &self.meta_all,
            kept: &self.meta_kept,
            regions: m.region_bins(),
            unassigned: m.unassigned,
        })
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
