//! The gene list and the gene panel: order, filter, profile, model, span.

use super::*;

#[test]
fn clicks_pick_the_panel_its_row_and_the_wheel_moves() {
    use ratatui::crossterm::event::{MouseButton, MouseEvent, MouseEventKind as M};
    let t = with_genes(table(M6A, 100));
    let mut p = picker(&t, None, SiteFilterArgs::default_values());
    let mut term = Terminal::new(TestBackend::new(120, 40)).unwrap();
    let mut at = |p: &mut SitePicker, kind, hit: Hit| {
        term.draw(|f| p.render(f)).unwrap();
        let (column, row) = p
            .hits
            .spot(&hit)
            .unwrap_or_else(|| panic!("{hit:?} not drawn"));
        p.mouse(MouseEvent {
            kind,
            column,
            row,
            modifiers: KeyModifiers::NONE,
        });
    };
    let click = M::Down(MouseButton::Left);
    at(&mut p, click, Hit::Genes(Some(2)));
    assert_eq!((p.panel, p.list.at), (Panel::Genes, 2));
    at(&mut p, M::ScrollDown, Hit::Genes(None));
    assert_eq!(p.list.at, 3);
    at(&mut p, click, Hit::Thresholds(Some(1)));
    assert_eq!((p.panel, p.focus), (Panel::Thresholds, 1));
    at(&mut p, click, Hit::Modality(0));
    assert_eq!(p.modality, 0);
    // The apply button asks first, as Ctrl+R; back goes back.
    at(&mut p, click, Hit::Apply);
    assert!(matches!(p.mode, Mode::Confirm));
    at(&mut p, click, Hit::Back);
    assert!(matches!(p.mode, Mode::Browse));
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
