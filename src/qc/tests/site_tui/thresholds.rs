//! The thresholds: knobs, stepping, keys, flags, drawing and saving.

use super::*;

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
    // Shift+Enter and Alt+Enter are not the apply key, nor plain Enter.
    for m in [KeyModifiers::SHIFT, KeyModifiers::ALT] {
        crate::tui::feed(&mut p, KeyEvent::new(KeyCode::Enter, m));
        assert!(!p.done() && matches!(p.mode, Mode::Browse) && p.enter_hint);
    }
    // Ctrl+J is no key: nothing is typed into the gene find.
    press(&mut p, KeyCode::Char('/'));
    crate::tui::feed(
        &mut p,
        KeyEvent::new(KeyCode::Char('j'), KeyModifiers::CONTROL),
    );
    assert!(matches!(p.mode, Mode::Find) && p.list.find.is_empty());
    press(&mut p, KeyCode::Esc);
    // Ctrl+S is not `s`: no save prompt opens.
    crate::tui::feed(
        &mut p,
        KeyEvent::new(KeyCode::Char('s'), KeyModifiers::CONTROL),
    );
    assert!(p.controls.footer().is_none());
    // Ctrl+R asks first; Esc goes back, Ctrl+R twice applies.
    apply_key(&mut p);
    assert!(!p.done() && matches!(p.mode, Mode::Confirm));
    press(&mut p, KeyCode::Enter);
    assert!(!p.done() && matches!(p.mode, Mode::Confirm));
    press(&mut p, KeyCode::Esc);
    assert!(!p.done() && matches!(p.mode, Mode::Browse));
    apply_key(&mut p);
    apply_key(&mut p);
    // With nothing to write, the next tick ends the session.
    assert!(matches!(p.mode, Mode::Writing));
    p.tick();
    assert!(p.done());
    let Some(Picked::Apply(got)) = p.decision.clone() else {
        panic!("Ctrl+R, Ctrl+R applies");
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
    assert!(!dir.path().join("view1.png").exists(), "PDF only");
    assert!(!p.done(), "saving does not leave the view");
    let svg = p.figure();
    assert!(svg.contains("kept") && svg.contains(figure::ACCENT));
}

#[test]
fn a_saved_view_keeps_its_shape_on_screen() {
    let t = table(M6A, 600);
    let mut p = picker(&t, Some(cells(t.len())), SiteFilterArgs::default_values());
    let (cw, ch) = figure::CELL;
    for (w, h) in [(160, 40), (100, 50)] {
        let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
        term.draw(|f| p.render(f)).unwrap();
        let (fw, fh) = figure::svg_size(&p.figure());
        let on_screen = w as f64 * cw / (h as f64 * ch);
        assert!((fw / fh - on_screen).abs() < 1e-6, "{w}x{h}: {fw}x{fh}");
        assert!(fw >= 720.0, "never narrower than its table and plots: {fw}");
    }
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
