use super::*;
use ratatui::backend::TestBackend;
use ratatui::crossterm::event::KeyModifiers;
use ratatui::Terminal;

fn rows() -> Vec<&'static str> {
    vec![
        "ENSG1_GENE1/m6a/chr1:100/methylated",
        "ENSG1_GENE1/m6a/chr1:100/unmethylated",
        "ENSG1_GENE1/m6a/chr1:180/methylated",
        "ENSG2_GENE2/m6a/chr2:50/methylated",
        "ENSG3_GENE3/m6a/chr1:900/methylated",
        "ENSG3_GENE3/m6a/chr1:950/methylated",
        "ENSG3_GENE3/m6a/chr1:990/methylated",
        "ENSG4_GENE4/count/spliced",
    ]
}

fn press(p: &mut GenePicker, code: KeyCode) {
    p.handle_key(KeyEvent::new(code, KeyModifiers::NONE));
}

#[test]
fn catalog_counts_converted_sites_per_gene() {
    let c = catalog_from_rows(rows());
    let names: Vec<&str> = c.iter().map(|e| e.gene.as_ref()).collect();
    assert_eq!(names, ["ENSG3_GENE3", "ENSG1_GENE1", "ENSG2_GENE2"]);
    assert_eq!((c[1].lo, c[1].hi, c[1].sites), (100, 180, 2));
    assert_eq!(c[2].chr.as_ref(), "chr2");
}

#[test]
fn typing_filters_and_enter_picks() {
    let c = catalog_from_rows(rows());
    let mut p = GenePicker::new(&c);
    for ch in "gene1".chars() {
        press(&mut p, KeyCode::Char(ch));
    }
    assert_eq!(p.shown.len(), 1);
    press(&mut p, KeyCode::Enter);
    assert_eq!(p.decision, Some(Choice::Gene(1)));

    let mut p = GenePicker::new(&c);
    press(&mut p, KeyCode::Char('x'));
    assert!(p.shown.is_empty());
    press(&mut p, KeyCode::Enter);
    assert_eq!(p.decision, None, "nothing to pick");
    press(&mut p, KeyCode::Esc);
    assert_eq!(p.shown.len(), 3, "Esc clears the filter first");
    press(&mut p, KeyCode::Down);
    press(&mut p, KeyCode::Down);
    press(&mut p, KeyCode::Down);
    assert_eq!(p.selected, 2);
    press(&mut p, KeyCode::Esc);
    assert_eq!(p.decision, Some(Choice::Quit), "then quits");
}

#[test]
fn renders_the_list() {
    let c = catalog_from_rows(rows());
    let mut p = GenePicker::new(&c);
    let mut term = Terminal::new(TestBackend::new(80, 10)).unwrap();
    term.draw(|f| p.render(f)).unwrap();
    let text: String = term
        .backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect();
    assert!(text.contains("ENSG3_GENE3") && text.contains("chr1:900-990"));
    assert!(text.contains("3 of 3 genes"));
    let mut tiny = Terminal::new(TestBackend::new(10, 3)).unwrap();
    tiny.draw(|f| p.render(f)).unwrap();
}

#[test]
fn a_typed_locus_opens_as_a_locus() {
    let c = catalog_from_rows(rows());
    let mut p = GenePicker::new(&c);
    for ch in "chr1:100-200".chars() {
        press(&mut p, KeyCode::Char(ch));
    }
    press(&mut p, KeyCode::Enter);
    assert_eq!(p.decision, Some(Choice::Locus("chr1:100-200".into())));
    p.set_filter("gene3");
    assert_eq!(p.filter(), "gene3");
    assert_eq!(p.shown.len(), 1);
}
