use extendr_api::prelude::*;

// mod aggregate;
// mod depth;
mod sift;
mod util;

use crate::util::bam::check_bam_index;

/// Sift through BAM records to identify potential variant sites.  The
/// resulting sites may not be necessarily true hits, but they can
/// form a good starting candidate pool.  Later, routines in
/// `aggregate` could revisit them all, regardless of significance
/// levels, and collect sufficient statistics for further tests.  So,
/// we just need to return a vector of (chr: Box<str>, lb: i64, ub:
/// i64)
///
/// @export
#[extendr]
fn compare_case_control_bam(fg_bam: &str, bg_bam: &str, block_size: Option<usize>) -> List {
    // sift::compare::search_case_control()

    println!("Establishing BAM File Sifters...");

    let fg_bai = (fg_bam.to_string() + ".bai").into_boxed_str();
    let bg_bai = (bg_bam.to_string() + ".bai").into_boxed_str();

    let _ = check_bam_index(fg_bam, Some(&fg_bai))
        .expect(&format!("failed to generate index for: {}", fg_bam));

    let _ = check_bam_index(bg_bam, Some(&bg_bai))
        .expect(&format!("failed to generate index for: {}", bg_bam));

    let mut data_fg = sift::sifter::BamSifter::from_file(fg_bam, &fg_bai, block_size);
    let mut data_bg = sift::sifter::BamSifter::from_file(bg_bam, &bg_bai, block_size);

    println!("Searching for variable positions");

    data_fg.sweep_variable_positions();
    data_bg.sweep_variable_positions();

    // update variable positions to each other
    data_bg.add_missed_positions(&data_fg);
    data_fg.add_missed_positions(&data_bg);

    println!("Collecting sufficient statistics");

    data_fg.populate_statistics();
    data_bg.populate_statistics();

    todo!("need to report");

    list!()
}

// Macro to generate exports.
// This ensures exported functions are registered with R.
// See corresponding C code in `entrypoint.c`.
extendr_module! {
    mod faba;
    fn compare_case_control_bam;
}
