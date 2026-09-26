use super::*;
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
    SitePicker::new("x", vec![view], start, Purpose::Apply)
}

fn press(p: &mut SitePicker, code: KeyCode) {
    p.handle_key(KeyEvent::new(code, KeyModifiers::NONE));
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
                assert_eq!(p.tally.alone[k as usize], dropped, "alone {k:?}");
            }
            // `first` refines the written reason: summed per reason it matches.
            let reasons = p.filter.reasons(&t, Some(&c));
            for &k in &p.view().criteria {
                let r = k.drop_reason();
                let want = reasons.iter().filter(|x| **x == Some(r)).count();
                let got: usize = Criterion::ALL
                    .iter()
                    .filter(|x| x.drop_reason() == r)
                    .map(|x| p.tally.first[*x as usize])
                    .sum();
                assert_eq!(got, want, "first {r:?}");
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
    press(&mut p, KeyCode::Enter);
    assert!(p.done());
    let got = p.decision.clone().flatten().expect("Enter confirms");
    assert_eq!(qc_flags(&got), qc_flags(&start));

    let mut q = picker(&t, None, start);
    press(&mut q, KeyCode::Char('q'));
    assert!(q.done());
    assert!(q.decision.clone().flatten().is_none());
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
    let mut p = SitePicker::new("x", views, start, Purpose::Apply);
    assert!(p.view().criteria.contains(&Criterion::MinFold));
    assert!(!p.view().criteria.contains(&Criterion::MinCells));
    focus_on(&mut p, Criterion::MinFold);
    press(&mut p, KeyCode::Tab);
    assert!(!p.view().criteria.contains(&Criterion::MinFold));
    assert!(!p.view().criteria.contains(&Criterion::MinLogOdds));
    assert_eq!(p.tally.kept, kept_linear(&atoi, None, &p.filter));

    focus_on(&mut p, Criterion::MinCoverage);
    press(&mut p, KeyCode::Right);
    let cov = p.filter.site_min_coverage;
    assert_eq!(p.tally.kept, kept_linear(&atoi, None, &p.filter));
    press(&mut p, KeyCode::BackTab);
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
    let mut p = SitePicker::new("input", vec![view], start, Purpose::Explore);
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

    // A tiny terminal must not panic either.
    let mut small = Terminal::new(TestBackend::new(30, 8)).unwrap();
    small.draw(|f| p.render(f)).unwrap();
}
