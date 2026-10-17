//! The metagene panel and the bars' weight.

use super::*;

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
fn an_annotation_without_coding_transcripts_is_named_in_the_panel() {
    let t = with_genes(table(M6A, 100));
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    p.set_meta(Meta::Ready(
        Annotation {
            views: vec![None],
            no_metagene: Some("x.gff: no coding transcript".into()),
            models: Default::default(),
        }
        .into(),
    ));
    assert_eq!(p.meta_status(), "x.gff: no coding transcript");
}
