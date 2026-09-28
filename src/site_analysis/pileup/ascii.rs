//! The ASCII pileup and its TSV.

use super::*;

pub(super) struct BinnedPileup {
    pub(super) gene: Box<str>,
    pub(super) chr: Box<str>,
    pub(super) bins: Vec<f64>,
    /// Distinct site coordinates (genomic order) backing the bins; used to
    /// place `+` axis markers and the right-side location list.
    pub(super) sites: Vec<i64>,
    pub(super) min_pos: i64,
    pub(super) max_pos: i64,
    pub(super) num_sites: usize,
    pub(super) track_label: Box<str>,
    pub(super) signal_name: &'static str,
}

pub(super) fn bin_positions_with_extent(
    positions: &[(i64, f64)],
    num_bins: usize,
    min_pos: i64,
    max_pos: i64,
    log_transform: bool,
) -> Vec<f64> {
    BinEdges::new(min_pos, max_pos, num_bins).bin(positions, log_transform)
}

/// Distinct coordinates from a position/value list (already sorted by
/// position), for the axis markers and right-side location list.
pub(super) fn distinct_positions(positions: &[(i64, f64)]) -> Vec<i64> {
    let mut out: Vec<i64> = positions.iter().map(|(p, _)| *p).collect();
    out.dedup();
    out
}

/// Group digits in threes for readability: `26781984` -> `26,781,984`.
pub(crate) fn fmt_thousands(n: i64) -> String {
    let digits = n.unsigned_abs().to_string();
    let bytes = digits.as_bytes();
    let mut out = String::with_capacity(digits.len() + digits.len() / 3);
    if n < 0 {
        out.push('-');
    }
    for (i, b) in bytes.iter().enumerate() {
        if i > 0 && (bytes.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(*b as char);
    }
    out
}

/// Map a coordinate to its bin column under the same binning rule as
/// [`bin_positions_with_extent`].
pub(super) fn pos_to_col(pos: i64, min_pos: i64, max_pos: i64, num_bins: usize) -> usize {
    BinEdges::new(min_pos, max_pos, num_bins).col_of(pos)
}

/// Right-side legend listing each site location top-to-bottom (genomic
/// order, mirroring the `+` marks left-to-right). First line is a title;
/// the list is capped to the rows available, with an overflow note.
pub(super) fn build_site_legend(pileup: &BinnedPileup, height: usize) -> Vec<String> {
    if pileup.sites.is_empty() {
        return Vec::new();
    }
    // Mixture rows have no chromosome — the "positions" are component ordinals.
    let kind = if pileup.chr.is_empty() || pileup.chr.as_ref() == "component" {
        "components".to_string()
    } else {
        format!("sites @ {}", pileup.chr)
    };
    let mut out = vec![format!("{} ({}):", kind, pileup.sites.len())];

    let capacity = height.saturating_sub(1); // rows left under the title
    let n = pileup.sites.len();
    if n <= capacity {
        out.extend(pileup.sites.iter().map(|&p| fmt_thousands(p)));
    } else {
        let shown = capacity.saturating_sub(1);
        out.extend(pileup.sites.iter().take(shown).map(|&p| fmt_thousands(p)));
        out.push(format!("... (+{} more)", n - shown));
    }
    out
}

pub(super) fn print_vertical_histogram(pileup: &BinnedPileup, height: usize) {
    let header = format!(
        "  {}  {}:{}-{}  [{}] signal: {}  sites: {}",
        pileup.gene,
        pileup.chr,
        fmt_thousands(pileup.min_pos),
        fmt_thousands(pileup.max_pos),
        pileup.track_label,
        pileup.signal_name,
        pileup.num_sites
    );

    let max_val = pileup.bins.iter().cloned().fold(0.0f64, f64::max);
    if max_val <= 0.0 {
        eprintln!("{}", header);
        eprintln!("  (no signal)");
        return;
    }

    eprintln!();
    eprintln!("{}", header);
    eprintln!();

    let max_label = format!("{:.1}", max_val);
    let label_width = max_label.len().max(4);

    // Locations listed down the right of the plot (one per row), so the
    // x-axis only needs `+` markers rather than crowded text.
    let legend = build_site_legend(pileup, height);

    for (i, row) in (1..=height).rev().enumerate() {
        let threshold = max_val * row as f64 / height as f64;

        let label = if row == height {
            format!("{:>w$.1}", max_val, w = label_width)
        } else if row == height / 2 {
            format!("{:>w$.1}", max_val / 2.0, w = label_width)
        } else if row == 1 {
            format!("{:>w$.1}", max_val / height as f64, w = label_width)
        } else {
            " ".repeat(label_width)
        };

        let mut bar = String::with_capacity(pileup.bins.len());
        for &val in &pileup.bins {
            if val >= threshold {
                bar.push('#');
            } else {
                bar.push(' ');
            }
        }

        let mut line = format!("{label} |{bar}");
        if let Some(entry) = legend.get(i) {
            line.push_str("   ");
            line.push_str(entry);
        }
        eprintln!("{}", line.trim_end());
    }

    // Axis line: `+` at every column that holds a site, `-` elsewhere.
    let mut axis = vec!['-'; pileup.bins.len()];
    for &pos in &pileup.sites {
        axis[pos_to_col(pos, pileup.min_pos, pileup.max_pos, pileup.bins.len())] = '+';
    }
    let axis: String = axis.into_iter().collect();
    eprintln!("{} +{}", " ".repeat(label_width), axis);
    eprintln!();
}

pub(super) fn write_pileup_tsv(tracks: &[&BinnedPileup], output: &str) -> anyhow::Result<()> {
    let mut writer = legume_numeric::matrix::common_io::open_buf_writer(output)?;

    for pileup in tracks {
        writeln!(
            writer,
            "#track={}\tgene={}\tchr={}\tmin_pos={}\tmax_pos={}\tsignal={}\tnum_sites={}",
            pileup.track_label,
            pileup.gene,
            pileup.chr,
            pileup.min_pos,
            pileup.max_pos,
            pileup.signal_name,
            pileup.num_sites
        )?;
        writeln!(writer, "bin\tgenomic_start\tgenomic_stop\tvalue")?;

        let span = (pileup.max_pos - pileup.min_pos).max(1);
        let num_bins = pileup.bins.len();

        for (i, &val) in pileup.bins.iter().enumerate() {
            let bin_start = pileup.min_pos + (i as i64 * span / num_bins as i64);
            let bin_stop = pileup.min_pos + ((i + 1) as i64 * span / num_bins as i64);
            writeln!(writer, "{}\t{}\t{}\t{:.4}", i, bin_start, bin_stop, val)?;
        }
    }

    writer.flush()?;
    Ok(())
}
