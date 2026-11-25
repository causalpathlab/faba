//! Site tables for the tests of the site views and the site panels: built
//! from rows, their derived columns as `read_site_table` makes them.

use super::sites::SiteTable;
use arrow::datatypes::Schema;
use arrow::record_batch::RecordBatch;
use std::sync::Arc;

/// A site's p-value, coverage, converted reads, and the control's two.
pub type Row = (f32, u64, u64, u64, u64);

/// A `modality` table of `rows` over `n_genes` genes, site `i` keyed
/// `chr1:{key(i)}` in gene `i % n_genes`.
pub fn site_table(
    modality: &str,
    n_genes: usize,
    rows: &[Row],
    key: impl Fn(usize) -> usize,
) -> SiteTable {
    let mut t = SiteTable {
        modality: modality.into(),
        batch: RecordBatch::new_empty(Arc::new(Schema::empty())),
        key: Vec::new(),
        gene_id: Vec::new(),
        n_genes,
        pv: Vec::new(),
        coverage: Vec::new(),
        converted: Vec::new(),
        control_coverage: Vec::new(),
        edit_ratio: Vec::new(),
        fold: Vec::new(),
        raw_log_odds: Vec::new(),
    };
    for (i, &(pv, cov, conv, ctl_cov, ctl_conv)) in rows.iter().enumerate() {
        let rate = |a: u64, n: u64| if n == 0 { 0.0 } else { a as f32 / n as f32 };
        let (w, m) = (rate(conv, cov), rate(ctl_conv, ctl_cov));
        t.key.push(format!("chr1:{}", key(i)).into());
        t.gene_id.push((i % n_genes) as u32);
        t.pv.push(pv);
        t.coverage.push(cov);
        t.converted.push(conv);
        t.control_coverage.push(ctl_cov);
        t.edit_ratio.push(w);
        t.fold.push(if m <= 0.0 { f32::INFINITY } else { w / m });
        t.raw_log_odds.push(faba::hypothesis_tests::log_odds_ratio(
            conv,
            cov - conv,
            ctl_conv,
            ctl_cov - ctl_conv,
        ));
    }
    t
}

/// Cells per site for `n` sites, 0 to 59.
pub fn cells(n: usize) -> Vec<usize> {
    (0..n).map(|i| i * 31 % 60).collect()
}
