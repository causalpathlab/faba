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
    p.handle_key(KeyEvent::new(code, KeyModifiers::NONE));
}

fn shift_enter(p: &mut SitePicker) {
    p.handle_key(KeyEvent::new(KeyCode::Enter, KeyModifiers::SHIFT));
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

#[test]
fn every_flag_has_a_criterion() {
    use clap::Args;
    let n = SiteFilterArgs::augment_args(clap::Command::new("x"))
        .get_arguments()
        .count();
    assert_eq!(n, Criterion::ALL.len());
    for (i, c) in Criterion::ALL.iter().enumerate() {
        assert_eq!(*c as usize, i, "ALL is in declaration order");
    }
}

#[test]
fn counts_agree_with_the_qc_rule() {
    let t = table(M6A, 1500);
    let c = cells(t.len());
    let mut p = picker(&t, Some(c.clone()), SiteFilterArgs::default_values());
    for crit in p.view().criteria.clone() {
        focus_on(&mut p, crit);
        for _ in 0..4 {
            press(&mut p, KeyCode::Right);
            assert_eq!(p.tally.kept, kept_linear(&t, Some(&c), &p.filter));
            for &k in &p.view().criteria {
                let dropped = t.len() - kept_linear(&t, Some(&c), &only(k, &p.filter));
                assert_eq!(p.tally.fail[k as usize], dropped, "fail {k:?}");
            }
        }
    }
}

#[test]
fn subset_is_the_sites_passing_every_other_knob() {
    let t = table(M6A, 800);
    let p = picker(&t, None, SiteFilterArgs::default_values());
    let mut others = p.filter.clone();
    Criterion::MaxPv.set(&mut others, Criterion::MaxPv.permissive());
    let want = kept_linear(&t, None, &others);
    assert_eq!(p.tally.subset.iter().sum::<usize>(), want);
    assert_eq!(p.column.hist.counts.iter().sum::<usize>(), t.len());
}

#[test]
fn a_stop_keeps_its_own_bin() {
    let t = table(M6A, 1200);
    let mut p = picker(&t, None, SiteFilterArgs::permissive());
    for crit in p.view().criteria.clone() {
        focus_on(&mut p, crit);
        for &(k, raw) in &p.column.stops.clone() {
            let mut f = SiteFilterArgs::permissive();
            crit.set(&mut f, raw);
            // Every site in bin `k` or on its kept side survives.
            for i in 0..t.len() {
                let b = p.column.hist.kmin + p.column.slot[i] as i32;
                let kept_side = if crit.keeps_high() { b >= k } else { b <= k };
                if kept_side {
                    assert!(f.reason(&t, i, None).is_none(), "{crit:?} bin {k} site {i}");
                }
            }
        }
    }
}

#[test]
fn stepping_walks_on_and_back_off() {
    let t = table(M6A, 1000);
    let mut p = picker(&t, None, SiteFilterArgs::permissive());
    let c = Criterion::MaxPv;
    assert!(c.is_off(&p.filter));
    press(&mut p, KeyCode::Right);
    assert!(!c.is_off(&p.filter), "right from off turns the knob on");
    let mut last = p.tally.kept;
    for _ in 0..10 {
        press(&mut p, KeyCode::Right);
        assert!(p.tally.kept <= last, "tightening never keeps more");
        last = p.tally.kept;
    }
    for _ in 0..200 {
        press(&mut p, KeyCode::Left);
    }
    assert!(c.is_off(&p.filter), "walking out the kept end turns it off");
    assert_eq!(p.tally.kept, t.len());

    // The upper edit-ratio bound keeps the low end: off is at the right.
    focus_on(&mut p, Criterion::MaxEditRatio);
    press(&mut p, KeyCode::Left);
    assert!(!Criterion::MaxEditRatio.is_off(&p.filter));
    for _ in 0..200 {
        press(&mut p, KeyCode::Right);
    }
    assert!(Criterion::MaxEditRatio.is_off(&p.filter));
}

#[test]
fn right_tightens_from_the_default_on_every_scale() {
    let t = table(M6A, 3000);
    for scale_presses in 0..3 {
        let mut p = picker(&t, None, SiteFilterArgs::default_values());
        for _ in 0..scale_presses {
            press(&mut p, KeyCode::Char('x'));
        }
        let mut last = p.filter.site_max_pv;
        for step in 0..5 {
            press(&mut p, KeyCode::Right);
            let now = p.filter.site_max_pv;
            assert!(now < last, "{:?} step {step}: {last} -> {now}", p.scale());
            last = now;
        }
    }
}

#[test]
fn keys_type_reset_off_and_decide() {
    let t = table(M6A, 300);
    let start = SiteFilterArgs::default_values();
    let mut p = picker(&t, Some(cells(t.len())), start.clone());
    focus_on(&mut p, Criterion::MinCoverage);
    for ch in "25".chars() {
        press(&mut p, KeyCode::Char(ch));
    }
    press(&mut p, KeyCode::Enter);
    assert_eq!(p.filter.site_min_coverage, 25);
    press(&mut p, KeyCode::Char('+'));
    assert_eq!(p.filter.site_min_coverage, 26);
    press(&mut p, KeyCode::Char('o'));
    assert_eq!(p.filter.site_min_coverage, 0);
    press(&mut p, KeyCode::Char('r'));
    assert_eq!(p.filter.site_min_coverage, start.site_min_coverage);
    assert!(!p.done());

    press(&mut p, KeyCode::Char('2'));
    press(&mut p, KeyCode::Esc);
    assert_eq!(p.filter.site_min_coverage, start.site_min_coverage);
    // Plain Enter does nothing, in the view or the confirmation.
    press(&mut p, KeyCode::Enter);
    assert!(!p.done() && matches!(p.mode, Mode::Browse));
    // Shift+Enter asks first; Esc goes back, Shift+Enter twice applies.
    shift_enter(&mut p);
    assert!(!p.done() && matches!(p.mode, Mode::Confirm));
    press(&mut p, KeyCode::Enter);
    assert!(!p.done() && matches!(p.mode, Mode::Confirm));
    press(&mut p, KeyCode::Esc);
    assert!(!p.done() && matches!(p.mode, Mode::Browse));
    shift_enter(&mut p);
    shift_enter(&mut p);
    // With nothing to write, the next tick ends the session.
    assert!(matches!(p.mode, Mode::Writing));
    p.tick();
    assert!(p.done());
    let Some(Picked::Apply(got)) = p.decision.clone() else {
        panic!("Shift+Enter, Shift+Enter applies");
    };
    assert_eq!(qc_flags(&got), qc_flags(&start));

    // With nothing changed, q and p leave at once.
    let mut q = picker(&t, None, start.clone());
    press(&mut q, KeyCode::Char('q'));
    assert!(matches!(q.decision, Some(Picked::Cancelled)));

    let mut r = picker(&t, None, start.clone());
    press(&mut r, KeyCode::Char('p'));
    assert!(matches!(r.decision, Some(Picked::PrintOnly(_))));

    // Esc never leaves the view: it is how every pop-up is closed.
    let mut e = picker(&t, None, start);
    press(&mut e, KeyCode::Esc);
    press(&mut e, KeyCode::Esc);
    assert!(!e.done() && matches!(e.mode, Mode::Browse));

    // With a threshold changed, q and p ask; Esc and other keys go back or
    // wait, and only the same key again leaves.
    press(&mut e, KeyCode::Char('+'));
    press(&mut e, KeyCode::Char('q'));
    assert!(matches!(e.mode, Mode::Leave(Picked::Cancelled)) && !e.done());
    press(&mut e, KeyCode::Char('p'));
    assert!(matches!(e.mode, Mode::Leave(_)) && !e.done());
    press(&mut e, KeyCode::Esc);
    assert!(matches!(e.mode, Mode::Browse) && !e.done());
    press(&mut e, KeyCode::Char('p'));
    assert!(matches!(e.mode, Mode::Leave(Picked::PrintOnly(_))));
    press(&mut e, KeyCode::Char('n'));
    assert!(matches!(e.mode, Mode::Browse) && !e.done());
    press(&mut e, KeyCode::Char('q'));
    press(&mut e, KeyCode::Char('q'));
    assert!(matches!(e.decision, Some(Picked::Cancelled)));
}

#[test]
fn modalities_share_thresholds_but_not_knobs() {
    let m6a = table(M6A, 400);
    let atoi = table(ATOI, 400);
    let start = SiteFilterArgs::default_values();
    let views = vec![
        SiteView::new(&m6a, None, &start),
        SiteView::new(&atoi, None, &start),
    ];
    let mut p = SitePicker::new("x", views, start, Meta::Unavailable("none".into()));
    assert!(p.view().criteria.contains(&Criterion::MinFold));
    assert!(!p.view().criteria.contains(&Criterion::MinCells));
    focus_on(&mut p, Criterion::MinFold);
    press(&mut p, KeyCode::Char('m'));
    assert!(!p.view().criteria.contains(&Criterion::MinFold));
    assert!(!p.view().criteria.contains(&Criterion::MinLogOdds));
    assert_eq!(p.tally.kept, kept_linear(&atoi, None, &p.filter));

    focus_on(&mut p, Criterion::MinCoverage);
    press(&mut p, KeyCode::Right);
    let cov = p.filter.site_min_coverage;
    assert_eq!(p.tally.kept, kept_linear(&atoi, None, &p.filter));
    press(&mut p, KeyCode::Char('M'));
    assert_eq!(p.criterion(), Criterion::MinCoverage);
    assert_eq!(p.filter.site_min_coverage, cov);
    assert_eq!(p.tally.kept, kept_linear(&m6a, None, &p.filter));
}

#[test]
fn flags_round_trip_through_the_cli() {
    let mut f = SiteFilterArgs::permissive();
    Criterion::MaxPv.set(&mut f, 0.0123);
    Criterion::MinCoverage.set(&mut f, 7.0);
    let flags = qc_flags(&f);
    let parsed = SiteFilterArgs::parse_flags(flags.split(' '));
    assert_eq!(qc_flags(&parsed), flags);
    assert!(flags.contains("--site-min-log-odds=-inf"));
}

#[test]
fn short_numbers_fit_the_table() {
    let c = Criterion::MaxPv;
    assert_eq!(c.fmt_short(0.05), "0.05");
    assert_eq!(c.fmt_short(0.000257887), "2.58e-4");
    assert_eq!(c.fmt_short(0.0123456), "0.0123");
    assert_eq!(c.fmt_short(1.0), "1");
    assert_eq!(Criterion::MinFold.fmt_short(12.3456), "12.3");
    assert_eq!(Criterion::MinCoverage.fmt_short(1418521.0), "1418521");
    for v in [1e-30, 0.00031, 0.5, 3.0, 45.123, 99999.0, 1e7] {
        assert!(c.fmt_short(v).len() <= 9, "{v}");
    }
}

#[test]
fn renders_every_knob_and_scale() {
    let t = table(M6A, 600);
    let start = SiteFilterArgs::default_values();
    let view = SiteView::new(&t, Some(cells(t.len())), &start);
    let mut p = SitePicker::new("input", vec![view], start, Meta::Unavailable("none".into()));
    let mut term = Terminal::new(TestBackend::new(120, 30)).unwrap();
    for _ in 0..p.view().criteria.len() {
        for _ in 0..3 {
            term.draw(|f| p.render(f)).unwrap();
            press(&mut p, KeyCode::Char('x'));
            press(&mut p, KeyCode::Char('y'));
        }
        press(&mut p, KeyCode::Down);
    }
    let text: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    assert!(text.contains("kept"));
    assert!(text.contains("print flags"));
    assert!(text.contains("save"));

    // A tiny terminal must not panic either.
    let mut small = Terminal::new(TestBackend::new(30, 8)).unwrap();
    small.draw(|f| p.render(f)).unwrap();
}

#[test]
fn saves_the_view_as_pdf_and_png() {
    let t = table(M6A, 600);
    let mut p = picker(&t, Some(cells(t.len())), SiteFilterArgs::default_values());
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("view1");
    press(&mut p, KeyCode::Char('s'));
    for _ in 0..20 {
        press(&mut p, KeyCode::Backspace);
    }
    for ch in prefix.to_str().unwrap().chars() {
        press(&mut p, KeyCode::Char(ch));
    }
    press(&mut p, KeyCode::Enter);
    assert!(dir.path().join("view1.pdf").exists());
    assert!(dir.path().join("view1.png").exists());
    assert!(!p.done(), "saving does not leave the view");
    let svg = p.figure();
    assert!(svg.contains("kept") && svg.contains(figure::ACCENT));
}

#[test]
fn draws_through_the_image_path() {
    let t = table(M6A, 400);
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    p.controls = Controls::new("x").with_picker(ratatui_image::picker::Picker::halfblocks());
    let mut term = Terminal::new(TestBackend::new(120, 30)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
    press(&mut p, KeyCode::Right);
    term.draw(|f| p.render(f)).unwrap();
    press(&mut p, KeyCode::Char('i'));
    assert!(p.controls.images().is_none());
    term.draw(|f| p.render(f)).unwrap();
}

#[test]
fn edit_ratio_rows_are_named_by_bound() {
    let t = table(M6A, 200);
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    let mut term = Terminal::new(TestBackend::new(120, 30)).unwrap();
    focus_on(&mut p, Criterion::MaxEditRatio);
    term.draw(|f| p.render(f)).unwrap();
    let text: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    assert!(text.contains("min edit ratio") && text.contains("max edit ratio"));
    assert!(text.contains("high is variant-like"));
}

#[test]
fn only_is_what_turning_the_threshold_off_keeps() {
    let t = table(M6A, 1500);
    let c = cells(t.len());
    let mut p = picker(&t, Some(c.clone()), SiteFilterArgs::default_values());
    let kept = kept_linear(&t, Some(&c), &p.filter);
    for k in p.view().criteria.clone() {
        focus_on(&mut p, k);
        let mut off = p.filter.clone();
        k.set(&mut off, k.permissive());
        let regained = kept_linear(&t, Some(&c), &off) - kept;
        assert_eq!(p.tally.only, regained, "{k:?}");
    }
    let mut term = Terminal::new(TestBackend::new(150, 40)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
    let text: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    assert!(text.contains("filtered out: sites failing this threshold."));
    assert!(text.contains("filtered out only by this "));
    assert!(text.contains("█ all sites "));
    assert!(text.contains("pass the other thresholds "));
    assert!(!text.contains("rescued"));
    assert!(!text.contains("reason"));
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

#[test]
fn the_metagene_follows_the_thresholds_on_a_fixed_axis() {
    let t = table(M6A, 1500);
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    p.weight = Weight::Sites;
    assert!(p.meta_counts().is_none());
    p.set_meta(Meta::ready(vec![Some(meta_layout(t.len()))], Vec::new()));
    let placed = |i: usize| i % 10 != 9;

    let before = p.meta_counts().unwrap();
    let (all, kept_before) = (before.all.to_vec(), before.kept.iter().sum::<f64>());
    focus_on(&mut p, Criterion::MinCoverage);
    for _ in 0..5 {
        press(&mut p, KeyCode::Right);
    }
    let after = p.meta_counts().unwrap();
    // The axis and the all-site profile do not move; the kept profile is
    // the kept sites, and it shrinks as the threshold tightens.
    assert_eq!(after.all, all);
    let kept_placed = (0..t.len())
        .filter(|&i| p.view().fails[i] == 0 && placed(i))
        .count();
    assert_eq!(after.kept.iter().sum::<f64>(), kept_placed as f64);
    assert!(after.kept.iter().sum::<f64>() < kept_before);
    assert_eq!(
        after.unassigned,
        (0..t.len()).filter(|&i| !placed(i)).count()
    );
}

#[test]
fn the_metagene_panel_says_why_it_is_empty_and_draws_when_ready() {
    let t = table(M6A, 400);
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    p.set_meta(Meta::Unavailable("no annotation: pass --gff".into()));
    let screen = |p: &mut SitePicker| {
        let mut term = Terminal::new(TestBackend::new(140, 44)).unwrap();
        term.draw(|f| p.render(f)).unwrap();
        term.backend()
            .buffer()
            .content()
            .iter()
            .map(|c| c.symbol())
            .collect::<String>()
    };
    let text = screen(&mut p);
    assert!(text.contains("metagene"));
    assert!(text.contains("no annotation: pass --gff"));

    p.set_meta(Meta::ready(vec![Some(meta_layout(t.len()))], Vec::new()));
    let text = screen(&mut p);
    assert!(text.contains("kept / all"));
    assert!(text.contains("CDS"));
    // The export carries the metagene too; read coverage by default.
    assert!(p
        .figure()
        .contains("m6a metagene · y: read coverage per bin"));
    assert!(text.contains("metagene · y: read coverage per bin"));
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

#[test]
fn genes_list_pinned_first_then_by_kept_sites_and_filter_by_symbol() {
    let t = with_genes(table(M6A, 100));
    let p = picker(&t, None, SiteFilterArgs::default_values());
    let kept = |p: &SitePicker| -> Vec<usize> {
        p.gene_list()
            .iter()
            .map(|&g| p.gene_kept(g as usize).0)
            .collect()
    };
    let sorted = |v: &[usize]| v.windows(2).all(|w| w[0] >= w[1]);
    assert!(sorted(&kept(&p)), "{:?}", kept(&p));

    // The order follows the thresholds, and the selection stays on its gene.
    let mut p = p;
    press(&mut p, KeyCode::Char(']'));
    let selected = p.gene();
    focus_on(&mut p, Criterion::MinCoverage);
    for _ in 0..5 {
        press(&mut p, KeyCode::Right);
    }
    assert!(sorted(&kept(&p)), "{:?}", kept(&p));
    assert_eq!(p.gene(), selected);

    let mut p = p.pin_genes(&["gene5".into(), "ENSG2_GENE2".into()]);
    let first: Vec<&str> = {
        let genes = p.view().genes.as_ref().unwrap();
        p.gene_list()[..2]
            .iter()
            .map(|&g| genes.symbol(g as usize))
            .collect()
    };
    assert_eq!(first, ["GENE5", "GENE2"]);

    // Typing narrows the list; Enter focuses the match and leaves the
    // search on the whole list, the gene still selected.
    press(&mut p, KeyCode::Char('/'));
    for ch in "ne3".chars() {
        press(&mut p, KeyCode::Char(ch));
    }
    assert_eq!(p.gene_list().len(), 1);
    press(&mut p, KeyCode::Enter);
    assert!(matches!(p.mode, Mode::Browse));
    assert_eq!(p.gene_list().len(), 7);
    let symbol = |p: &SitePicker| {
        let g = p.gene().unwrap();
        p.view().genes.as_ref().unwrap().symbol(g).to_string()
    };
    assert_eq!(symbol(&p), "GENE3");
    assert_eq!(p.panel, Panel::Genes);

    // The arrows move through the matches; letters still type.
    press(&mut p, KeyCode::Char('/'));
    for ch in "gene".chars() {
        press(&mut p, KeyCode::Char(ch));
    }
    assert_eq!(p.list.find, "gene");
    let first = symbol(&p);
    press(&mut p, KeyCode::Down);
    press(&mut p, KeyCode::Down);
    let third = symbol(&p);
    assert_ne!(third, first);
    press(&mut p, KeyCode::Enter);
    assert_eq!(symbol(&p), third);

    // Esc goes back to the gene selected before the search.
    press(&mut p, KeyCode::Char('/'));
    for ch in "ne5".chars() {
        press(&mut p, KeyCode::Char(ch));
    }
    press(&mut p, KeyCode::Esc);
    assert_eq!(p.gene_list().len(), 7);
    assert_eq!(symbol(&p), third);
}

#[test]
fn the_gene_profile_follows_the_thresholds() {
    let t = with_genes(table(M6A, 1500));
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    p.weight = Weight::Sites;
    press(&mut p, KeyCode::Char(']'));
    let g = p.gene().unwrap();
    let before = p.gene_profile(40).unwrap();
    let (kept, all) = p.gene_kept(g);
    assert_eq!(before.all.iter().sum::<usize>(), all);
    assert_eq!(before.kept.iter().sum::<usize>(), kept);
    // No gene model yet: the span is the sites' own, and no exon track.
    assert!(before.model.is_none() && before.exonic(0).is_none());

    focus_on(&mut p, Criterion::MinCoverage);
    for _ in 0..5 {
        press(&mut p, KeyCode::Right);
    }
    let after = p.gene_profile(40).unwrap();
    assert_eq!(after.all, before.all);
    assert_eq!(after.kept.iter().sum::<usize>(), p.gene_kept(g).0);
    assert!(after.kept.iter().sum::<usize>() < kept);
}

#[test]
fn the_gene_panel_draws_its_model_once_the_annotation_arrives() {
    let t = with_genes(table(M6A, 400));
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    let g = p.gene().unwrap();
    let key = p.view().genes.as_ref().unwrap().keys[g].clone();
    let model = GeneModel {
        chr: "chr1".into(),
        lo: 900,
        hi: 3100,
        forward: true,
        exons: vec![(900, 1500), (2500, 3100)],
        symbol: key.split_once('_').unwrap().1.into(),
        key,
    };
    p.set_meta(Meta::ready(vec![None], vec![model]));
    let profile = p.gene_profile(44).unwrap();
    let exonic = |b| profile.exonic(b).unwrap();
    assert!(exonic(0) && exonic(43) && !exonic(22));

    let mut term = Terminal::new(TestBackend::new(150, 44)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
    let text: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    let sum = |v: &[usize]| v.iter().sum::<usize>();
    let kept = format!("kept {} of {}", sum(&profile.kept), sum(&profile.all));
    assert!(text.contains(&format!(
        "chr1:900-3100 (+) · {kept} · y: read coverage per "
    )));
    assert!(text.contains("genes: kept / all sites"));
    assert!(text.contains("[ ] move  / find"));
    assert!(text.contains("█ kept read coverage   c: converted reads"));
    assert!(text.contains("▬"));
    assert!(p.figure().contains("chr1:900-3100 (+)"));
}

#[test]
fn c_cycles_the_bars_from_read_coverage_to_converted_reads_and_sites() {
    let t = with_genes(table(M6A, 1500));
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    p.set_meta(Meta::ready(vec![Some(meta_layout(t.len()))], Vec::new()));
    let g = p.gene().unwrap();
    let rows = p.view().genes.as_ref().unwrap().rows[g].clone();
    let placed = |i: &usize| i % 10 != 9;
    // Each measure, per gene and on the metagene: the gene's total, its kept
    // total, and the metagene's total over the placed sites.
    let check = |p: &SitePicker, unit: &str, w: &dyn Fn(usize) -> u64| {
        let profile = p.gene_profile(40).unwrap();
        assert_eq!(profile.unit, unit);
        let want: u64 = rows.iter().map(|&i| w(i as usize)).sum();
        assert_eq!(profile.all.iter().sum::<usize>() as u64, want, "{unit}");
        let kept_want: u64 = rows
            .iter()
            .filter(|&&i| p.view().fails[i as usize] == 0)
            .map(|&i| w(i as usize))
            .sum();
        assert_eq!(profile.kept.iter().sum::<usize>() as u64, kept_want);
        let placed_want: u64 = (0..t.len()).filter(placed).map(w).sum();
        let meta: f64 = p.meta_counts().unwrap().all.iter().sum();
        assert_eq!(meta, placed_want as f64, "{unit}");
    };
    let screen = |p: &mut SitePicker| {
        let mut term = Terminal::new(TestBackend::new(150, 44)).unwrap();
        term.draw(|f| p.render(f)).unwrap();
        term.backend()
            .buffer()
            .content()
            .iter()
            .map(|c| c.symbol())
            .collect::<String>()
    };

    check(&p, "read coverage", &|i| t.coverage[i]);
    let text = screen(&mut p);
    assert!(text.contains("metagene · y: read coverage per bin"));
    assert!(text.contains("c: converted reads"));

    press(&mut p, KeyCode::Char('c'));
    check(&p, "converted reads", &|i| t.converted[i]);
    let text = screen(&mut p);
    assert!(text.contains("metagene · y: converted reads per bin"));
    assert!(text.contains("█ kept converted reads"));
    assert!(text.contains("c: sites"));

    press(&mut p, KeyCode::Char('c'));
    check(&p, "sites", &|_| 1);

    press(&mut p, KeyCode::Char('c'));
    check(&p, "read coverage", &|i| t.coverage[i]);
}

#[test]
fn an_annotation_arriving_on_its_thread_fills_the_metagene() {
    let t = with_genes(table(M6A, 400));
    let n = t.len();
    let pending = Meta::Pending(std::thread::spawn(move || {
        Ok(Annotation {
            views: vec![Some(meta_layout(n))],
            no_metagene: None,
            models: Default::default(),
        })
    }));
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    p.set_meta(pending);
    assert!(p.meta_counts().is_none());
    while !p.tick() {
        std::thread::sleep(std::time::Duration::from_millis(5));
    }
    let m = p.meta_counts().unwrap();
    assert_eq!(m.all.len(), m.regions.iter().sum::<usize>());
    assert_eq!(m.kept.len(), m.all.len());
    let mut term = Terminal::new(TestBackend::new(150, 44)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
}

#[test]
fn sites_outside_the_gene_model_widen_its_span() {
    let t = with_genes(table(M6A, 400));
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    let g = p.gene().unwrap();
    let key = p.view().genes.as_ref().unwrap().keys[g].clone();
    // A model narrower than the sites, which span [1000, 3000).
    let model = GeneModel {
        chr: "chr1".into(),
        lo: 1500,
        hi: 2500,
        forward: true,
        exons: vec![(1500, 2500)],
        symbol: "GENE".into(),
        key,
    };
    p.set_meta(Meta::ready(vec![None], vec![model]));
    let genes = p.view().genes.as_ref().unwrap();
    let pos: Vec<i64> = genes.rows[g]
        .iter()
        .map(|&i| genes.pos.value(i as usize))
        .collect();
    let (min, max) = (*pos.iter().min().unwrap(), *pos.iter().max().unwrap());
    assert!(
        min < 1500 && max >= 2500,
        "the fixture's sites overhang the model"
    );
    let profile = p.gene_profile(40).unwrap();
    let edges = &profile.edges;
    assert_eq!((edges.min_pos, edges.max_pos), (min, max + 1));
    assert!(profile.exonic(0) == Some(false) && profile.exonic(39) == Some(false));
}

#[test]
fn an_annotation_without_coding_transcripts_is_named_in_the_panel() {
    let t = with_genes(table(M6A, 100));
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    p.set_meta(Meta::Ready(Annotation {
        views: vec![None],
        no_metagene: Some("x.gff: no coding transcript".into()),
        models: Default::default(),
    }));
    assert_eq!(p.meta_status(), "x.gff: no coding transcript");
}

#[test]
fn the_confirmation_recaps_output_changes_and_every_modality() {
    let (m6a, atoi) = (table(M6A, 600), table(ATOI, 300));
    let start = SiteFilterArgs::default_values();
    let views = vec![
        SiteView::new(&m6a, None, &start),
        SiteView::new(&atoi, None, &start),
    ];
    let mut p = SitePicker::new("x", views, start, Meta::Unavailable("none".into()));
    p.writer.output = "out_qc";
    focus_on(&mut p, Criterion::MinCoverage);
    for _ in 0..3 {
        press(&mut p, KeyCode::Right);
    }
    shift_enter(&mut p);
    assert!(matches!(p.mode, Mode::Confirm));
    let text: String = p
        .confirm_lines()
        .iter()
        .map(|l| l.to_string() + "\n")
        .collect();
    assert!(text.contains("out_qc"));
    let coverage: Vec<&str> = text.lines().filter(|l| l.contains("coverage")).collect();
    assert!(
        coverage[0].starts_with("* ") && coverage[0].contains("(was ≥ 3)"),
        "{coverage:?}"
    );
    assert!(text.lines().any(|l| l.trim_start().starts_with(M6A)));
    assert!(text.lines().any(|l| l.trim_start().starts_with(ATOI)));

    let mut term = Terminal::new(TestBackend::new(150, 44)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
    let screen: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    assert!(screen.contains("apply these thresholds?"));
    assert!(screen.contains("apply and write"));
    press(&mut p, KeyCode::Char('n'));
    assert!(matches!(p.mode, Mode::Browse) && !p.done());
}

#[test]
fn tab_moves_between_panels_and_the_arrows_follow() {
    let t = with_genes(table(M6A, 300));
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    assert_eq!(p.panel, Panel::Thresholds);
    let (focus, gene) = (p.focus, p.list.at);
    press(&mut p, KeyCode::Down);
    assert_eq!((p.focus, p.list.at), (focus + 1, gene));

    press(&mut p, KeyCode::Tab);
    assert_eq!(p.panel, Panel::Genes);
    press(&mut p, KeyCode::Down);
    press(&mut p, KeyCode::Down);
    assert_eq!((p.focus, p.list.at), (focus + 1, gene + 2));
    press(&mut p, KeyCode::Up);
    assert_eq!(p.list.at, gene + 1);

    press(&mut p, KeyCode::BackTab);
    assert_eq!(p.panel, Panel::Thresholds);
    press(&mut p, KeyCode::Up);
    assert_eq!((p.focus, p.list.at), (focus, gene + 1));
}

#[test]
fn confirming_writes_in_the_view_and_shows_progress() {
    let progress = Progress::default();
    let started = std::cell::RefCell::new(None);
    let t = table(M6A, 300);
    let start = SiteFilterArgs::default_values();
    let mut p = picker(&t, None, start.clone());
    p.writer = Writer {
        output: "out_qc",
        progress: &progress,
        start: Some(Box::new(|f, figures| {
            *started.borrow_mut() = Some((f, figures))
        })),
    };

    // Only Shift+Enter applies: `A`, `G`, `y` and plain Enter do not.
    for code in [KeyCode::Char('A'), KeyCode::Char('G')] {
        press(&mut p, code);
        assert!(matches!(p.mode, Mode::Browse));
    }
    shift_enter(&mut p);
    assert!(matches!(p.mode, Mode::Confirm));
    for code in [
        KeyCode::Char('y'),
        KeyCode::Char('A'),
        KeyCode::Char('G'),
        KeyCode::Enter,
    ] {
        press(&mut p, code);
        assert!(matches!(p.mode, Mode::Confirm));
    }
    shift_enter(&mut p);
    assert!(matches!(p.mode, Mode::Writing) && !p.done());
    // The figures are drawn, and the writer started, on the next tick.
    assert!(started.borrow().is_none());
    assert!(p.tick() && !p.done());
    let (f, figures) = started.borrow_mut().take().expect("the writer starts");
    assert_eq!(qc_flags(&f), qc_flags(&start));
    // A figure per knob, named after its flag, and the view left as it was.
    assert_eq!(figures.len(), p.view().criteria.len());
    assert!(figures.iter().any(|f| f.stem == format!("{M6A}_max-pv")));
    assert!(figures.iter().all(|f| f.svg.starts_with("<svg")));
    assert_eq!((p.modality, p.focus), (0, 0));

    // Keys and Ctrl-C wait while the fileset is written.
    press(&mut p, KeyCode::Char('q'));
    p.interrupt();
    assert!(matches!(p.mode, Mode::Writing) && !p.done());

    progress.plan(4);
    progress.next("first");
    progress.next("m6a_sites.parquet");
    assert!(p.tick() && !p.done());
    let mut term = Terminal::new(TestBackend::new(150, 44)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
    let screen: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    assert!(screen.contains("writing"));
    assert!(screen.contains("1 / 4"));
    assert!(screen.contains("m6a_sites.parquet"));

    progress.finish();
    p.tick();
    assert!(p.done());
    assert!(matches!(p.decision, Some(Picked::Apply(_))));
}

#[test]
fn figures_cover_every_modality_and_leave_the_view_as_it_was() {
    let m6a = with_genes(table(M6A, 400));
    let atoi = with_genes(table(ATOI, 300));
    let start = SiteFilterArgs::default_values();
    let views = vec![
        SiteView::new(&m6a, None, &start),
        SiteView::new(&atoi, None, &start),
    ];
    let mut p = SitePicker::new("x", views, start, Meta::Unavailable("none".into()));
    press(&mut p, KeyCode::Char('m'));
    press(&mut p, KeyCode::Down);
    // A filter on the view, matching one gene of the modality on screen.
    p.set_find(|f| f.push_str("ne3"));
    assert_eq!(p.gene_list().len(), 1);
    let (modality, focus, gene) = (p.modality, p.focus, p.gene());

    let figures = p.figures();
    // The filter is put back for the view, and kept out of the figures:
    // every one, of either modality, has its gene panel.
    assert_eq!(p.list.find, "ne3");
    assert!(figures.iter().all(|f| f.svg.contains(" bp")));
    let n: usize = p.views.iter().map(|v| v.criteria.len()).sum();
    assert_eq!(figures.len(), n);
    for v in &p.views {
        let m = &*v.table.modality;
        assert!(figures.iter().any(|f| f.stem.starts_with(&format!("{m}_"))));
    }
    assert_eq!((p.modality, p.focus, p.gene()), (modality, focus, gene));
}

#[test]
fn g_browses_for_an_annotation_and_reads_the_one_picked() {
    let tmp = tempfile::tempdir().unwrap();
    let root = tmp.path().to_path_buf();
    std::fs::create_dir_all(root.join("ref")).unwrap();
    for f in ["ref/genes.gtf.gz", "ref/notes.txt"] {
        std::fs::write(root.join(f), b"").unwrap();
    }
    let t = table(M6A, 50);
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    // A run record names an annotation that is not on this machine.
    p.annotation = AnnotationSource {
        gff: Some("/elsewhere/gencode.gtf.gz".into()),
        dir: root.clone(),
        ..Default::default()
    };
    p.set_meta(Meta::start(&p.annotation));
    assert_eq!(p.meta_status(), "/elsewhere/gencode.gtf.gz: not found");
    let mut term = Terminal::new(TestBackend::new(150, 44)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
    let text: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    assert!(text.contains("choose an annotation (GTF/GFF)"));

    // The browser opens where the sites are, since the named one is absent,
    // and lists directories and annotations only.
    press(&mut p, KeyCode::Char('g'));
    let Mode::Gff(b) = &p.mode else {
        panic!("no browser")
    };
    assert_eq!(b.cwd, root);
    let at = b.entries.iter().position(|e| e.name == "ref").unwrap();
    term.draw(|f| p.render(f)).unwrap();
    let Mode::Gff(b) = &mut p.mode else {
        unreachable!()
    };
    b.at = at;
    press(&mut p, KeyCode::Enter);
    let Mode::Gff(b) = &p.mode else {
        panic!("left the browser")
    };
    let names: Vec<&str> = b.entries.iter().map(|e| e.name.as_str()).collect();
    assert_eq!(names, ["..", "genes.gtf.gz"]);

    // Picking the file reads it, and the view reports it for the record.
    press(&mut p, KeyCode::Down);
    press(&mut p, KeyCode::Enter);
    assert!(matches!(p.mode, Mode::Browse));
    let picked = root.join("ref/genes.gtf.gz");
    assert_eq!(p.gff(), Some(picked.to_str().unwrap()));
    assert!(!p.meta_status().ends_with("not found"));

    // Esc leaves the browser without reading anything.
    press(&mut p, KeyCode::Char('g'));
    press(&mut p, KeyCode::Esc);
    assert!(matches!(p.mode, Mode::Browse));
    assert_eq!(p.gff(), Some(picked.to_str().unwrap()));
}

#[test]
fn a_pick_that_fails_to_read_keeps_the_annotation_given_for_the_record() {
    let t = table(M6A, 50);
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    let tmp = tempfile::tempdir().unwrap();
    let given = tmp.path().join("given.gtf");
    std::fs::write(&given, b"").unwrap();
    let given: Box<str> = given.to_str().unwrap().into();
    p.annotation = AnnotationSource {
        gff: Some(given.clone()),
        given: Some(given.clone()),
        ..Default::default()
    };
    p.load_gff("/nowhere/picked.gtf".into());
    assert_eq!(p.meta_status(), "/nowhere/picked.gtf: not found");
    assert_eq!(p.gff(), Some(&*given));
}
