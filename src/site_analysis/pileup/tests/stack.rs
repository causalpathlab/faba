use super::*;
use crate::site_analysis::pileup::tui::Track;
use data_beans::interactive::ui::Scale;
use ratatui::backend::TestBackend;
use ratatui::crossterm::event::KeyModifiers;
use ratatui::Terminal;

/// Sites every 97 bp over a 100 kb gene body, `scale` reads each.
fn positions(scale: f64) -> Vec<(i64, f64)> {
    (0..1000).map(|i| (1_000_000 + i * 97, scale)).collect()
}

const EXTENT: (i64, i64) = (1_000_000, 1_100_000);

fn stack<'a>(a: &'a [(i64, f64)], b: &'a [(i64, f64)]) -> StackView<'a> {
    let one = |name: &str, m: &'a [(i64, f64)]| {
        PileupView::new(
            name,
            "chr1",
            vec![Track::single("matrix", "sum", m, false)],
            EXTENT,
        )
    };
    StackView::new(vec![one("GENE1", a), one("GENE2", b)])
}

fn press(s: &mut StackView, code: KeyCode) {
    crate::tui::feed(s, KeyEvent::new(code, KeyModifiers::NONE));
}

fn screen(s: &mut StackView, w: u16, h: u16) -> String {
    let mut term = Terminal::new(TestBackend::new(w, h)).unwrap();
    term.draw(|f| s.render(f)).unwrap();
    term.backend()
        .buffer()
        .content()
        .iter()
        .map(|c| c.symbol())
        .collect()
}

#[test]
fn genes_stack_top_to_bottom_with_every_key_shown() {
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    let text = screen(&mut s, 80, 40);
    let (one, two) = (text.find("GENE1").unwrap(), text.find("GENE2").unwrap());
    assert!(one < two, "the first marked on top");
    assert!(text.contains("▸ keys"), "says which gene has the keys");
    let missing = crate::tui::missing_keys(&s.footer_line(), &text);
    assert!(missing.is_empty(), "{missing:?}");
}

#[test]
fn keys_go_to_one_gene_and_how_bars_draw_to_all() {
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    screen(&mut s, 80, 40);
    press(&mut s, KeyCode::Char('y'));
    assert!(s.views.iter().all(|v| v.y_scale == s.views[0].y_scale));
    assert_ne!(
        s.views[1].y_scale,
        Scale::Linear,
        "the scale applies to both"
    );

    // Tab moves the keys on; `Y` locks one y axis for both.
    press(&mut s, KeyCode::Tab);
    assert_eq!(s.focus, 1);
    let before = s.views[0].cursor;
    press(&mut s, KeyCode::Right);
    assert_eq!(
        s.views[0].cursor, before,
        "only the gene with the keys moves"
    );
    press(&mut s, KeyCode::Char('Y'));
    let lock = s.views[1].y_lock;
    assert!(
        lock.is_some() && s.views[0].y_lock == lock,
        "one axis for both"
    );

    press(&mut s, KeyCode::Char('/'));
    assert!(
        s.status.is_some() && s.exit.is_none(),
        "searching needs one gene"
    );
    press(&mut s, KeyCode::Esc);
    assert!(matches!(s.exit, Some(Exit::Genes)));
}

#[test]
fn the_figure_is_the_genes_one_under_another() {
    let (a, b) = (positions(1.0), positions(3.0));
    let s = stack(&a, &b);
    let svg = s.figure();
    assert!(svg.contains("GENE1") && svg.contains("GENE2"));
    let one = s.views[0].figure();
    let height = |svg: &str| -> f64 {
        let at = svg.find(" height=\"").unwrap() + 9;
        svg[at..].split('"').next().unwrap().parse().unwrap()
    };
    let two = s.views[1].figure();
    assert_eq!(height(&svg), height(&one) + height(&two));
    // It renders, as a save would.
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("both");
    crate::figure::save(&svg, out.to_str().unwrap()).unwrap();
    assert!(dir.path().join("both.pdf").exists());
}

#[test]
fn y_locks_every_gene_to_the_tallest() {
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    screen(&mut s, 100, 40);
    let tallest = s.views[1].bins().shared;
    // From the shorter gene, the lock is still the taller one's height.
    press(&mut s, KeyCode::Char('Y'));
    assert_eq!(s.views[0].bins().shared, tallest);
    assert_eq!(s.views[1].bins().shared, tallest);
    press(&mut s, KeyCode::Char('Y'));
    assert!(
        s.views.iter().all(|v| v.y_lock.is_none()),
        "Y again unlocks all"
    );
}

#[test]
fn every_gene_redraws_when_the_scale_or_bars_change() {
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    screen(&mut s, 100, 40);
    // Drawn once: no plot waits to be redrawn.
    for v in &mut s.views {
        v.plots.iter_mut().for_each(|p| *p = Default::default());
    }
    press(&mut s, KeyCode::Tab);
    press(&mut s, KeyCode::Char('y'));
    let other = &s.views[0];
    assert_ne!(other.y_scale, Scale::Linear);
    assert!(
        other
            .plots
            .iter()
            .all(crate::figure::term::PlotImage::is_stale),
        "the gene without the keys redraws too"
    );
}
