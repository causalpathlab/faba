use rust_htslib::bam::{self, Read};
use std::hash::Hash;
use std::path::Path;
use std::thread;

/// BAM file sample name
///
#[derive(Debug, PartialEq, Eq, Hash, Clone)]
#[allow(dead_code)]
pub enum BamSample {
    Combined,
    Barcode(Box<str>),
}

/// Display sample names
///
#[allow(dead_code)]
impl std::fmt::Display for BamSample {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            BamSample::Combined => write!(f, "."),
            BamSample::Barcode(barcode) => write!(f, "{}", barcode),
        }
    }
}

/// Check random access BAM index
///
#[allow(dead_code)]
pub fn check_bam_index(
    bam_file_name: &str,
    idx_file_name: Option<&str>,
) -> anyhow::Result<Box<str>> {
    // log::info!("Checking BAM index");

    let idx_file = match idx_file_name {
        Some(x) => String::from(x),
        None => format!("{}.bai", bam_file_name),
    };

    if Path::new(&idx_file).exists() {
        return Ok(idx_file.into_boxed_str());
    }

    let ncore = thread::available_parallelism()
        .expect("failed to figure out number of cores")
        .get();

    // log::info!(
    //     "Creating a new index file {} using {} cores",
    //     &idx_file,
    //     &ncore
    // );

    // need to build an index for this bam file
    bam::index::build(
        bam_file_name,
        Some(&idx_file),
        bam::index::Type::Bai,
        ncore as u32,
    )?;

    Ok(idx_file.into_boxed_str())
}
