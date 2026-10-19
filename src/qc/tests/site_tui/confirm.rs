//! Confirming and writing, figures, picking an annotation, the footer
//! and Ctrl+C.

use super::*;

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
    apply_key(&mut p);
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

    // Only Ctrl+R applies: `A`, `G`, `y` and plain Enter do not.
    for code in [KeyCode::Char('A'), KeyCode::Char('G')] {
        press(&mut p, code);
        assert!(matches!(p.mode, Mode::Browse));
    }
    apply_key(&mut p);
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
    apply_key(&mut p);
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
    press(&mut p, KeyCode::Char('a'));
    let Mode::Gff(b) = &p.mode else {
        panic!("no browser")
    };
    assert_eq!(b.cwd, root);
    let _ = b;
    term.draw(|f| p.render(f)).unwrap();
    let Mode::Gff(b) = &mut p.mode else {
        unreachable!()
    };
    assert!(b.list.select("ref"));
    press(&mut p, KeyCode::Enter);
    let Mode::Gff(b) = &p.mode else {
        panic!("left the browser")
    };
    let names: Vec<&str> = b.list.shown().map(|e| e.name.as_str()).collect();
    assert_eq!(names, ["..", "genes.gtf.gz"]);

    // Picking the file reads it, and the view reports it for the record.
    press(&mut p, KeyCode::Down);
    press(&mut p, KeyCode::Enter);
    assert!(matches!(p.mode, Mode::Browse));
    let picked = root.join("ref/genes.gtf.gz");
    assert_eq!(p.gff(), Some(picked.to_str().unwrap()));
    assert!(!p.meta_status().ends_with("not found"));

    // Esc leaves the browser without reading anything.
    press(&mut p, KeyCode::Char('a'));
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

#[test]
fn every_footer_key_shows_at_80_columns() {
    let t = table(M6A, 600);
    let start = SiteFilterArgs::default_values();
    let view = SiteView::new(&t, Some(cells(t.len())), &start);
    let mut p = SitePicker::new("input", vec![view], start, Meta::Unavailable("none".into()));
    for panel in [KeyCode::Tab, KeyCode::Tab] {
        let mut term = Terminal::new(TestBackend::new(80, 30)).unwrap();
        term.draw(|f| p.render(f)).unwrap();
        let screen: String = term
            .backend()
            .buffer()
            .content()
            .iter()
            .map(|c| c.symbol())
            .collect();
        let missing = crate::tui::missing_keys(&p.footer_line(), &screen);
        assert!(missing.is_empty(), "{missing:?}");
        press(&mut p, panel);
    }
}

#[test]
fn ctrl_c_asks_first_as_q_does() {
    let t = table(M6A, 600);
    let start = SiteFilterArgs::default_values();
    let mut p = picker(&t, None, start);
    press(&mut p, KeyCode::Char('+'));
    p.interrupt();
    assert!(
        matches!(p.mode, Mode::Leave(Picked::Cancelled)) && !p.done(),
        "asks"
    );
    p.interrupt();
    assert!(p.done(), "the second leaves");
}

#[test]
fn a_kill_leaves_changed_thresholds_without_asking() {
    let t = table(M6A, 300);
    let changed = || {
        let mut p = picker(&t, None, SiteFilterArgs::default_values());
        focus_on(&mut p, Criterion::MinCoverage);
        press(&mut p, KeyCode::Right);
        assert_ne!(qc_flags(&p.filter), qc_flags(&p.initial), "changed");
        p
    };
    // Ctrl-C asks first, the thresholds having changed ...
    let mut asked = changed();
    asked.interrupt();
    assert!(matches!(asked.mode, Mode::Leave(_)) && !asked.done());
    // ... a kill does not: the change is dropped and nothing written.
    let mut killed = changed();
    View::quit(&mut killed);
    assert!(killed.done());
    assert!(matches!(killed.decision, Some(Picked::Cancelled)));
}

#[test]
fn a_kill_lets_a_write_under_way_finish() {
    let progress = Progress::default();
    let t = table(M6A, 300);
    let start = SiteFilterArgs::default_values();
    let mut p = picker(&t, None, start);
    p.writer = Writer {
        output: "out_qc",
        progress: &progress,
        start: Some(Box::new(|_, _| {})),
    };
    apply_key(&mut p);
    apply_key(&mut p);
    assert!(p.tick() && matches!(p.mode, Mode::Writing));
    // The kill waits for the write: the view is done once it is, with
    // the thresholds it wrote.
    View::quit(&mut p);
    assert!(matches!(p.mode, Mode::Writing) && !p.done());
    progress.finish();
    p.tick();
    assert!(p.done());
    assert!(matches!(p.decision, Some(Picked::Apply(_))));
}
