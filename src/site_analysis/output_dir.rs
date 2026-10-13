//! What `pileup`, `metagene` and `pwm` need from a faba output directory (from
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

/// The output written without `-o` in batch: `{sites}{suffix}` in the
/// current directory after the first site table, `_combined` for several.
pub fn default_output(site_files: &[Box<str>], suffix: &str) -> Box<str> {
    let name = crate::qc::layout::file_name(&site_files[0]);
    let stem = name.strip_suffix(".parquet").unwrap_or(&name);
    let together = if site_files.len() > 1 {
        "_combined"
    } else {
        ""
    };
    format!("{stem}{together}{suffix}").into()
}

/// The view's title: the site tables' file names.
pub fn title(site_files: &[Box<str>]) -> String {
    let names: Vec<String> = site_files
        .iter()
        .map(|f| crate::qc::layout::file_name(f).to_string())
        .collect();
    names.join(" + ")
}

/// The site tables to read: each of `given` (a site table, or an output
/// directory or run record holding one of `modality`), or with none given,
/// those chosen in a browser for `command`. `None` when that is cancelled.
pub fn site_tables(
    command: &'static str,
    batch: bool,
    given: Vec<Box<str>>,
    modality: Option<&str>,
) -> anyhow::Result<Option<Vec<Box<str>>>> {
    let mut given = given;
    if given.is_empty() {
        let ask = crate::site_analysis::input_picker::ask_inputs(
            command,
            batch,
            "site tables",
            |n| crate::qc::layout::site_table_modality(n).is_some(),
            false,
        );
        let Some(chosen) = ask? else {
            return Ok(None);
        };
        given = chosen.paths;
    }
    given
        .iter()
        .map(|g| site_table_of(g, modality))
        .collect::<anyhow::Result<_>>()
        .map(Some)
}

/// `given`, or the site table of the output directory it names.
fn site_table_of(given: &str, modality: Option<&str>) -> anyhow::Result<Box<str>> {
    let Some(dir) = crate::run_record::output_dir_of(given) else {
        return Ok(given.into());
    };
    let sites = output_sites(&dir, modality)?;
    let table = sites.site_table.ok_or_else(|| {
        anyhow::anyhow!("{} has no {}_sites.parquet", dir.display(), sites.modality)
    })?;
    log::info!("{}: reading {table}", dir.display());
    Ok(table)
}

#[cfg(test)]
#[path = "tests/output_dir.rs"]
mod tests;
