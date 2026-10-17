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
    let missing = crate::tui::missing_keys(&s.footer_line(s.views[0].rows_hidden()), &text);
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

#[test]
fn bar_cycles_the_layout_and_the_figure_follows() {
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    let width = |svg: &str| -> f64 {
        let at = svg.find(" width=\"").unwrap() + 8;
        svg[at..].split('"').next().unwrap().parse().unwrap()
    };
    let one = width(&s.views[0].figure());
    assert_eq!(width(&s.figure()), one, "one under another");
    press(&mut s, KeyCode::Char('|'));
    assert_eq!(s.arrange, Arrange::Row);
    assert_eq!(width(&s.figure()), 2.0 * one, "side by side");
    // Side by side on screen, too: GENE2 to the right of GENE1, same rows.
    let text = screen(&mut s, 160, 30);
    let row_of = |i: usize| i / 160;
    let (g1, g2) = (text.find("GENE1").unwrap(), text.find("GENE2").unwrap());
    assert_eq!(row_of(g1), row_of(g2), "{text}");
    press(&mut s, KeyCode::Char('|'));
    assert_eq!(s.arrange, Arrange::Grid);
    press(&mut s, KeyCode::Char('|'));
    assert_eq!(s.arrange, Arrange::Column);
}

#[test]
fn a_click_gives_the_gene_the_keys_and_moves_its_cursor() {
    use ratatui::crossterm::event::{MouseButton, MouseEvent, MouseEventKind};
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    screen(&mut s, 100, 40);
    let lower = s.views[1].drawn;
    let x = s.views[1].chart_x + 40;
    let click = MouseEvent {
        kind: MouseEventKind::Down(MouseButton::Left),
        column: x,
        row: lower.y + 3,
        modifiers: KeyModifiers::NONE,
    };
    let before = s.views[1].cursor;
    assert!(s.mouse(click));
    // Let go where it was pressed: a click.
    s.mouse(MouseEvent {
        kind: MouseEventKind::Up(MouseButton::Left),
        ..click
    });
    assert_eq!(s.focus, 1, "the clicked gene has the keys");
    assert_ne!(s.views[1].cursor, before, "its cursor moved to the bar");
    let (lo, hi) = s.views[1].bar_range(40);
    assert!((lo..=hi).contains(&s.views[1].cursor));
}

#[test]
fn shifted_zoom_keys_zoom_every_gene() {
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    screen(&mut s, 100, 40);
    let whole: Vec<_> = s.views.iter().map(|v| v.window).collect();
    press(&mut s, KeyCode::Char('='));
    assert_ne!(
        s.views[0].window, whole[0],
        "= zooms the gene with the keys"
    );
    assert_eq!(s.views[1].window, whole[1], "and only it");
    press(&mut s, KeyCode::Char('+'));
    assert_ne!(s.views[1].window, whole[1], "+ zooms every gene");
    press(&mut s, KeyCode::Char(')'));
    let now: Vec<_> = s.views.iter().map(|v| v.window).collect();
    assert_eq!(now, whole, ") shows every gene whole");
    press(&mut s, KeyCode::Char('+'));
    press(&mut s, KeyCode::Char('_'));
    let now: Vec<_> = s.views.iter().map(|v| v.window).collect();
    assert_eq!(now, whole, "_ zooms every gene back out");
    assert!(s.exit.is_none());
}

#[test]
fn a_drag_or_the_wheel_pans_the_gene_under_the_pointer() {
    use ratatui::crossterm::event::{MouseButton, MouseEvent, MouseEventKind};
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    screen(&mut s, 100, 40);
    // Zoomed in, so there is room to pan.
    press(&mut s, KeyCode::Char('+'));
    press(&mut s, KeyCode::Char('+'));
    screen(&mut s, 100, 40);
    let at = |kind, column, row| MouseEvent {
        kind,
        column,
        row,
        modifiers: KeyModifiers::NONE,
    };
    let v = &s.views[1];
    let (x, y, before) = (v.chart_x + 50, v.drawn.y + 3, v.window);
    // As far as the pointer: the window's span a column.
    let ten_columns = 10 * (before.1 - before.0) / v.columns as i64;
    s.mouse(at(MouseEventKind::Down(MouseButton::Left), x, y));
    s.mouse(at(MouseEventKind::Up(MouseButton::Left), x - 10, y));
    assert_eq!(s.focus, 1);
    let after = s.views[1].window;
    assert_eq!(
        after.0 - before.0,
        ten_columns,
        "dragged left: the view moves right"
    );
    assert_eq!(after.1 - after.0, before.1 - before.0, "the same width");

    // The wheel pans the gene under it, keys or not.
    let top = s.views[0].window;
    let y0 = s.views[0].drawn.y + 3;
    s.mouse(at(MouseEventKind::ScrollDown, x, y0));
    assert!(s.views[0].window.0 > top.0, "a notch to the right");
    s.mouse(at(MouseEventKind::ScrollLeft, x, y0));
    assert_eq!(s.views[0].window, top, "and back");
}

#[test]
fn a_saved_gene_keeps_its_shape_on_screen() {
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    let size = crate::figure::svg_size;
    screen(&mut s, 120, 40);
    let tall = size(&s.views[0].figure());
    press(&mut s, KeyCode::Char('|'));
    screen(&mut s, 120, 40);
    let wide = size(&s.views[0].figure());
    // One under another: each gene is the full width; side by side, half.
    assert_eq!(tall.0, 120.0 * crate::figure::CELL.0);
    assert_eq!(wide.0, 60.0 * crate::figure::CELL.0);
    assert!(wide.1 > tall.1, "side by side, each gene is taller");
}

#[test]
fn a_press_outside_every_gene_lets_go_of_a_grab() {
    use ratatui::crossterm::event::{MouseButton, MouseEvent, MouseEventKind};
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    press(&mut s, KeyCode::Char('+'));
    screen(&mut s, 100, 40);
    let at = |kind, column, row| MouseEvent {
        kind,
        column,
        row,
        modifiers: KeyModifiers::NONE,
    };
    let v = &s.views[0];
    let (x, y, before) = (v.chart_x + 30, v.drawn.y + 3, v.window);
    // Grabbed, then its release lost (let go outside the terminal).
    s.mouse(at(MouseEventKind::Down(MouseButton::Left), x, y));
    // A press and release on the footer, outside every gene.
    s.mouse(at(MouseEventKind::Down(MouseButton::Left), x + 20, 39));
    s.mouse(at(MouseEventKind::Up(MouseButton::Left), x + 20, 39));
    assert_eq!(s.views[0].window, before, "no stale grab pans the gene");
}

/// A margin and `n` groups, each converted above unconverted in one row.
fn group_tracks(m: &[(i64, f64)], n: usize) -> Vec<Track<'_>> {
    const LABELS: [&str; 2] = ["on", "off"];
    (0..=n)
        .flat_map(|_| {
            let up = Track::single(LABELS[0], "reads", m, false).converted_above();
            [up, Track::single(LABELS[1], "reads", m, false)]
        })
        .collect()
}

/// Two genes of a margin and 12 groups each, both opened at track `from`.
fn group_stack(m: &[(i64, f64)], from: usize) -> StackView<'_> {
    let one = |name: &str| {
        let mut v = PileupView::new(name, "chr1", group_tracks(m, 12), EXTENT);
        v.track_from = from;
        v
    };
    StackView::new(vec![one("GENE1"), one("GENE2")])
}

/// The first track of the first scrolled row each gene shows.
fn first_scrolled(s: &StackView) -> Vec<usize> {
    s.views
        .iter()
        .map(|v| match v.shown_rows(&v.rows())[1] {
            super::super::Row::Track(k) | super::super::Row::Mirror(k, _) => k,
            r => panic!("{r:?}"),
        })
        .collect()
}

#[test]
fn split_or_mirrored_every_gene_stays_at_its_group() {
    let m = positions(1.0);
    let mut s = group_stack(&m, 8);
    screen(&mut s, 200, 120);
    assert_eq!(first_scrolled(&s), vec![8, 8]);
    press(&mut s, KeyCode::Char('v'));
    screen(&mut s, 200, 120);
    assert_eq!(first_scrolled(&s), vec![8, 8], "split, both genes");
    press(&mut s, KeyCode::Char('v'));
    screen(&mut s, 200, 120);
    assert_eq!(first_scrolled(&s), vec![8, 8], "mirrored again");
}

#[test]
fn the_first_frame_names_the_scroll_keys() {
    let m = positions(1.0);
    let mut s = group_stack(&m, 0);
    let first = screen(&mut s, 200, 60);
    assert!(first.contains("[/]"), "{first}");
    assert_eq!(screen(&mut s, 200, 60), first, "as every frame after");
}
