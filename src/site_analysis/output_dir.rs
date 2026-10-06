//! What `pileup` and `metagene` need from a faba output directory (from
//! `faba run`, a producer, or `faba qc`), so they can be pointed at the
//! directory, or at a run record in it, rather than at its files.

use std::path::Path;

use crate::qc::layout::{scan_input_dir, SITE_MODALITIES};

/// One modality's sites in an output directory.
#[derive(Debug)]
pub struct OutputSites {
    pub modality: Box<str>,
    /// `(batch, path)` of each `{batch}_{modality}_site` matrix, by batch.
    pub matrices: Vec<(Box<str>, Box<str>)>,
    /// `{modality}_sites.parquet`, when there is one.
    pub site_table: Option<Box<str>>,
}

/// `modality`'s sites in `dir`, or (without one) those of the first of
/// [`SITE_MODALITIES`] the directory has.
pub fn output_sites(dir: &Path, modality: Option<&str>) -> anyhow::Result<OutputSites> {
    let layout = scan_input_dir(&dir.to_string_lossy())?;
    let of = |m: &str| {
        let mut matrices: Vec<(Box<str>, Box<str>)> = layout
            .matrices
            .iter()
            .filter(|f| f.site_modality() == Some(m))
            .map(|f| (f.batch.clone(), f.path.clone()))
            .collect();
        matrices.sort();
        OutputSites {
            modality: m.into(),
            matrices,
            site_table: layout.site_tables.get(m).cloned(),
        }
    };
    let has = |s: &OutputSites| !s.matrices.is_empty() || s.site_table.is_some();
    if let Some(m) = modality {
        let sites = of(m);
        anyhow::ensure!(
            has(&sites),
            "{} has no {m} site matrices or site table",
            dir.display()
        );
        return Ok(sites);
    }
    SITE_MODALITIES
        .iter()
        .map(|m| of(m))
        .find(has)
        .ok_or_else(|| {
            anyhow::anyhow!(
                "{} has no site matrices or site tables ({})",
                dir.display(),
                SITE_MODALITIES.join(", ")
            )
        })
}

#[cfg(test)]
#[path = "tests/output_dir.rs"]
mod tests;
