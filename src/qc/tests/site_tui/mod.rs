//! Tests of the site-threshold picker. The fixtures are here; the tests
//! are by what they cover.

use super::*;
use crate::site_analysis::miami::genemodel::GeneModel;
use crate::site_analysis::site_io::GenomicSite;
use arrow::datatypes::Schema;
use arrow::record_batch::RecordBatch;
use data_beans::aux::feature_rows::{ATOI, M6A};
use ratatui::backend::TestBackend;
use ratatui::crossterm::event::KeyModifiers;
use ratatui::Terminal;
use std::sync::Arc;

mod confirm;
mod genes;
mod metagene;
mod report;
mod thresholds;

/// A deterministic site table: `n` sites over seven genes, with p-values
/// down to 0, clean and noisy controls, and a spread of depths.
fn table(modality: &str, n: usize) -> SiteTable {
    let mut t = SiteTable {
        modality: modality.into(),
        batch: RecordBatch::new_empty(Arc::new(Schema::empty())),
        key: Vec::new(),
        gene_id: Vec::new(),
        n_genes: 7,
        pv: Vec::new(),
        coverage: Vec::new(),
        converted: Vec::new(),
        control_coverage: Vec::new(),
        edit_ratio: Vec::new(),
        fold: Vec::new(),
        raw_log_odds: Vec::new(),
    };
    for i in 0..n {
        let coverage = 1 + (i * 37 % 200) as u64;
        let converted = (i * 13 % 50) as u64 % (coverage + 1);
        let (control_cov, control_conv) = if i % 5 == 0 {
            (10, 0)
        } else {
            (20, (i % 7) as u64)
        };
        let rate_w = converted as f32 / coverage as f32;
        let rate_m = control_conv as f32 / control_cov as f32;
        t.key.push(format!("chr1:{i}").into());
        t.gene_id.push((i % 7) as u32);
        t.pv.push(if i % 97 == 0 {
            0.0
        } else {
            10f32.powf(-((i * 7919 % 1000) as f32) / 100.0)
        });
        t.coverage.push(coverage);
        t.converted.push(converted);
        t.control_coverage.push(control_cov);
        t.edit_ratio.push(rate_w);
        t.fold.push(if rate_m <= 0.0 {
            f32::INFINITY
        } else {
            rate_w / rate_m
        });
        t.raw_log_odds.push(faba::hypothesis_tests::log_odds_ratio(
            converted,
            coverage - converted,
            control_conv,
            control_cov - control_conv,
        ));
    }
    t
}

fn cells(n: usize) -> Vec<usize> {
    (0..n).map(|i| i * 31 % 60).collect()
}

/// A picker over one table, applying on Enter.
fn picker<'a>(
    t: &'a SiteTable,
    n_cells: Option<Vec<usize>>,
    start: SiteFilterArgs,
) -> SitePicker<'a> {
    let view = SiteView::new(t, n_cells, &start);
    SitePicker::new("x", vec![view], start, Meta::Unavailable("none".into()))
}

fn press(p: &mut SitePicker, code: KeyCode) {
    crate::tui::feed(p, KeyEvent::new(code, KeyModifiers::NONE));
}

fn apply_key(p: &mut SitePicker) {
    crate::tui::feed(p, KeyEvent::new(KeyCode::Char('r'), KeyModifiers::CONTROL));
}

fn kept_linear(t: &SiteTable, c: Option<&[usize]>, f: &SiteFilterArgs) -> usize {
    f.reasons(t, c).iter().filter(|r| r.is_none()).count()
}

fn focus_on(p: &mut SitePicker, c: Criterion) {
    while p.criterion() != c {
        press(p, KeyCode::Down);
    }
}

/// `f` with only `c` taken from `from`, every other threshold off.
fn only(c: Criterion, from: &SiteFilterArgs) -> SiteFilterArgs {
    let mut f = SiteFilterArgs::permissive();
    c.set(&mut f, c.get(from));
    f
}

/// A layout for `n` sites spread over one forward two-exon transcript on
/// chr1 (exons 1000..1199 and 1500..1699, CDS 1100..1599), one site off it
/// in every ten.
fn meta_layout(n: usize) -> MetaLayout {
    use genomic_data::gff::{FeatureType, GeneId, GeneSymbol, GeneType, GffRecord, TranscriptId};
    use genomic_data::sam::Strand;
    let rec = |feature_type, start, stop| GffRecord {
        seqname: "chr1".into(),
        feature_type,
        start,
        stop,
        strand: Strand::Forward,
        gene_id: GeneId::Ensembl("GENE1".into()),
        gene_name: GeneSymbol::Symbol("GENE1".into()),
        gene_type: GeneType::CodingGene,
        transcript_id: TranscriptId::Ensembl("T1".into()),
    };
    let records = [
        rec(FeatureType::Exon, 1000, 1199),
        rec(FeatureType::Exon, 1500, 1699),
        rec(FeatureType::CDS, 1100, 1199),
        rec(FeatureType::CDS, 1500, 1599),
        rec(FeatureType::StopCodon, 1600, 1602),
    ];
    let exonic: Vec<i64> = (1000..1200).chain(1500..1700).collect();
    let sites: Vec<GenomicSite> = (0..n)
        .map(|i| GenomicSite {
            chr: "chr1".into(),
            position: if i % 10 == 9 {
                5000
            } else {
                exonic[i * 7 % exonic.len()] - 1
            },
            strand: Strand::Forward,
        })
        .collect();
    MetaModels::from_records(&records)
        .layout(&sites, META_BINS)
        .expect("sites on the transcript")
}

/// `t` with a batch carrying the gene and position columns the gene view
/// reads: gene `GENE{g}` for the fixture's dense id `g`, sites spread along
/// `[1000, 3000)`.
fn with_genes(mut t: SiteTable) -> SiteTable {
    use arrow::array::{ArrayRef, Int64Array, StringArray};
    use arrow::datatypes::{DataType, Field};
    let gene: Vec<String> = t
        .gene_id
        .iter()
        .map(|g| format!("ENSG{g}_GENE{g}"))
        .collect();
    let pos: Vec<i64> = (0..t.len() as i64).map(|i| 1000 + i * 37 % 2000).collect();
    let schema = Schema::new(vec![
        Field::new("gene", DataType::Utf8, false),
        Field::new("primary_pos", DataType::Int64, false),
    ]);
    let columns: Vec<ArrayRef> = vec![
        Arc::new(StringArray::from(gene)),
        Arc::new(Int64Array::from(pos)),
    ];
    t.batch = RecordBatch::try_new(Arc::new(schema), columns).unwrap();
    t
}
