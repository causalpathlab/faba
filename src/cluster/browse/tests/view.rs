//! The browser on a test terminal: keys, the mouse and what is drawn.

use super::{cell_of, heads, loaded, mouse, press, screen, typed};
use crate::cluster::browse::view::{Look, Pane, Pin, Ramp, Sort};
use crate::cluster::browse::ModulesView;
use data_beans::interactive::ui::Screen;
use ratatui::crossterm::event::{KeyCode, KeyModifiers, MouseButton, MouseEventKind};

#[test]
fn keys_walk_the_levels_and_back() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    assert_eq!((v.pane, v.module), (Pane::Modules, 1));
    press(&mut v, KeyCode::Down);
    assert_eq!(v.module, 0);
    press(&mut v, KeyCode::Down);
    assert_eq!(v.module, 0, "the last module stays");
    press(&mut v, KeyCode::Enter);
    assert_eq!(v.pane, Pane::Sites);
    press(&mut v, KeyCode::End);
    assert_eq!(v.site, 2);
    press(&mut v, KeyCode::Enter);
    assert_eq!(v.pane, Pane::Genes);
    press(&mut v, KeyCode::Enter);
    assert!(!v.done(), "nothing pinned, nothing to open");
    press(&mut v, KeyCode::Esc);
    assert_eq!(v.pane, Pane::Sites);
    press(&mut v, KeyCode::Esc);
    press(&mut v, KeyCode::Esc);
    assert!(!v.done(), "Esc at the top does not quit");
    press(&mut v, KeyCode::Up);
    assert_eq!(
        (v.module, v.site),
        (1, 0),
        "a new module starts at its first site"
    );
    press(&mut v, KeyCode::Char('q'));
    assert!(v.done());
}

#[test]
fn space_pins_or_selects_delete_unpins_and_x_clears() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    press(&mut v, KeyCode::Down);
    press(&mut v, KeyCode::Enter);
    // GENE1's first site, then GENE2's.
    press(&mut v, KeyCode::Char(' '));
    press(&mut v, KeyCode::End);
    press(&mut v, KeyCode::Char(' '));
    let gene1 = Pin { role: 0, gene: 0 };
    let gene2 = Pin { role: 0, gene: 1 };
    assert_eq!(v.pins, vec![gene1, gene2]);
    // Space on a pinned gene selects it; it never unpins.
    press(&mut v, KeyCode::Home);
    press(&mut v, KeyCode::Char(' '));
    assert_eq!((v.pins.len(), v.pin), (2, 0));
    // Pins come from any module, in the order pinned.
    press(&mut v, KeyCode::Esc);
    press(&mut v, KeyCode::Up);
    press(&mut v, KeyCode::Enter);
    press(&mut v, KeyCode::Char(' '));
    assert_eq!(v.pins[2], Pin { role: 1, gene: 0 });
    // In the strip: ←/→ walk the sites across the genes, Enter opens one.
    press(&mut v, KeyCode::Enter);
    assert_eq!(v.pane, Pane::Genes);
    press(&mut v, KeyCode::Up);
    press(&mut v, KeyCode::Up);
    assert_eq!((v.pin, v.pin_site), (0, 0));
    press(&mut v, KeyCode::Right);
    press(&mut v, KeyCode::Right);
    assert_eq!((v.pin, v.pin_site), (1, 2), "into GENE2");
    press(&mut v, KeyCode::Enter);
    assert_eq!(v.take_open().as_deref(), Some("GENE2"));
    // Delete unpins the gene in focus.
    press(&mut v, KeyCode::Delete);
    assert_eq!(v.pins, vec![gene1, Pin { role: 1, gene: 0 }]);
    press(&mut v, KeyCode::Char('X'));
    assert!(v.pins.is_empty());
}

#[test]
fn a_drag_needs_a_modifier_to_pin_move_or_unpin() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    press(&mut v, KeyCode::Down);
    let s = screen(&mut v, 180, 44);
    let site = cell_of(&s, "GENE2 50");
    let genes = (90, 38);
    let (down, up) = (
        MouseEventKind::Down(MouseButton::Left),
        MouseEventKind::Up(MouseButton::Left),
    );
    // A plain press and release over the strip only selects.
    mouse(&mut v, down, site, KeyModifiers::NONE);
    mouse(&mut v, up, genes, KeyModifiers::NONE);
    assert!(v.pins.is_empty());
    assert_eq!(v.site, 2, "but the site is selected");
    // Alt-drag a site row onto the strip: its gene is pinned.
    mouse(&mut v, down, site, KeyModifiers::ALT);
    mouse(&mut v, up, genes, KeyModifiers::ALT);
    assert_eq!(v.pins, vec![Pin { role: 0, gene: 1 }]);
    assert_eq!(v.pane, Pane::Genes);
    let s = screen(&mut v, 180, 44);
    assert!(s.contains("1 pinned"), "{s}");
    // A plain click on the pinned gene keeps it.
    let pinned = cell_of(&s, "GENE2 M1");
    mouse(&mut v, down, (pinned.0 + 2, pinned.1), KeyModifiers::NONE);
    mouse(&mut v, up, (5, 5), KeyModifiers::NONE);
    assert_eq!(v.pins.len(), 1, "a click never unpins");
    // Alt-drag it off the strip: unpinned.
    mouse(&mut v, down, (pinned.0 + 2, pinned.1), KeyModifiers::ALT);
    mouse(&mut v, up, (5, 5), KeyModifiers::ALT);
    assert!(v.pins.is_empty());
}

#[test]
fn v_cycles_a_grid_and_y_the_percent_scale_and_heatmaps_go_diagonal() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    press(&mut v, KeyCode::Char('y'));
    assert_eq!(v.ramp, Ramp::Sqrt, "y acts on a % heatmap only");
    press(&mut v, KeyCode::Char('v'));
    assert_eq!(
        v.looks,
        [Look::Percent, Look::Bars],
        "the focused grid only"
    );
    press(&mut v, KeyCode::Char('y'));
    assert_eq!(v.ramp, Ramp::Linear);
    press(&mut v, KeyCode::Char('v'));
    assert_eq!(v.looks[0], Look::Z);
    let s = screen(&mut v, 160, 40);
    assert!(s.contains("diagonal") && s.contains("z heatmap"), "{s}");
    // Numbers on the highlighted row only: its z, signed.
    let numbered = s
        .lines()
        .filter(|l| l.contains("+1.0") || l.contains("-1.0"))
        .count();
    assert_eq!(numbered, 1, "{s}");
    press(&mut v, KeyCode::Char('v'));
    assert_eq!(v.looks[0], Look::Bars);
}

#[test]
fn o_sorts_the_sites_by_each_column_then_reads_then_cells() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    press(&mut v, KeyCode::Down);
    press(&mut v, KeyCode::Enter);
    assert_eq!(v.listing(), &[0, 1, 2], "the module's order");
    press(&mut v, KeyCode::Char('o'));
    assert_eq!(v.sort, Sort::Column(0));
    // The first column shown, high % first.
    let c = m.columns[0];
    let f = |i: usize| m.members[0][i].reads[c].fraction().unwrap_or(-1.0);
    let list = v.listing();
    assert!(list.windows(2).all(|w| f(w[0]) >= f(w[1])), "{list:?}");
    press(&mut v, KeyCode::Char('O'));
    let list = v.listing();
    assert!(
        list.windows(2).all(|w| f(w[0]) <= f(w[1])),
        "reversed: {list:?}"
    );
    press(&mut v, KeyCode::Char('o'));
    press(&mut v, KeyCode::Char('o'));
    assert_eq!(v.sort, Sort::Reads);
    assert_eq!(v.listing(), &[0, 1, 2], "40, 20 and 6 reads");
    press(&mut v, KeyCode::Char('o'));
    assert_eq!(v.sort, Sort::Cells);
    press(&mut v, KeyCode::Char('o'));
    assert_eq!(v.sort, Sort::Module);
}

#[test]
fn slash_jumps_to_a_gene_as_it_is_typed_and_esc_goes_back() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    press(&mut v, KeyCode::Char('/'));
    typed(&mut v, "gene2");
    assert_eq!((v.module, v.site, v.pane), (0, 2, Pane::Sites));
    let s = screen(&mut v, 160, 40);
    assert!(s.contains("find gene: gene2"), "{s}");
    // j and k are letters while searching.
    press(&mut v, KeyCode::Backspace);
    typed(&mut v, "3");
    assert_eq!((v.module, v.site), (1, 0));
    press(&mut v, KeyCode::Esc);
    assert_eq!((v.module, v.site, v.pane), (1, 0, Pane::Modules), "back");
    press(&mut v, KeyCode::Char('/'));
    typed(&mut v, "gene1");
    press(&mut v, KeyCode::Enter);
    assert_eq!((v.module, v.site), (0, 0), "kept");
    press(&mut v, KeyCode::Char('j'));
    assert_eq!(v.site, 1, "keys are keys again");
}

#[test]
fn the_screen_shows_the_grids_the_counts_and_genes_side_by_side() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    let s = screen(&mut v, 180, 40);
    // Module M2 first, then M1; each block's methylated %.
    let row = |s: &str, name: &str| -> String {
        s.lines()
            .find(|l| l.contains(name) && l.contains("genes"))
            .unwrap()
            .to_string()
    };
    let (m2, m1) = (row(&s, "M2 "), row(&s, "M1 "));
    assert!(m2.contains(" 75") && m2.contains(" 25"), "{m2}");
    assert!(m1.contains(" 64") && m1.contains(" 10"), "{m1}");
    assert!(s.find(&m2).unwrap() < s.find(&m1).unwrap());
    // The strip starts empty.
    assert!(s.contains("Space on a site pins its gene here"), "{s}");
    // M1's sites with their counts over every cell; two genes pinned side
    // by side.
    press(&mut v, KeyCode::Down);
    press(&mut v, KeyCode::Enter);
    press(&mut v, KeyCode::Char(' '));
    press(&mut v, KeyCode::End);
    press(&mut v, KeyCode::Char(' '));
    press(&mut v, KeyCode::Home);
    let s = screen(&mut v, 180, 40);
    assert!(s.contains("reads") && s.contains("cells"), "{s}");
    let gene1 = s.lines().find(|l| l.contains("GENE1 100")).unwrap();
    assert!(gene1.contains("40") && gene1.contains("9"), "{gene1}");
    let names = s
        .lines()
        .find(|l| l.contains("GENE1") && l.contains("GENE2") && !l.contains("100"))
        .unwrap();
    assert!(
        names.find("GENE1").unwrap() < names.find("GENE2").unwrap(),
        "{names}"
    );
    assert!(s.contains("all"), "the all-cells track");
    assert!(s.contains("• GENE1 100"), "pinned genes are marked: {s}");
    assert!(
        s.contains("modules") && s.contains("faba"),
        "levels and the mark"
    );
}

#[test]
fn a_long_list_scrolls_to_keep_the_site_in_focus() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    press(&mut v, KeyCode::Down);
    press(&mut v, KeyCode::Enter);
    // So short only one site row fits in the list.
    screen(&mut v, 160, 18);
    press(&mut v, KeyCode::Down);
    press(&mut v, KeyCode::Down);
    screen(&mut v, 160, 18);
    assert_eq!(v.site, 2);
    assert!(v.from <= 2 && 2 < v.from + 3);
}

#[test]
fn cell_modules_narrow_before_any_is_hidden() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    // Room for both cell modules at full width: bars with their %.
    let s = screen(&mut v, 180, 40);
    let m2 = s
        .lines()
        .find(|l| l.contains("M2 ") && l.contains("genes"))
        .unwrap();
    assert!(m2.contains(" 75") && m2.contains(" 25"), "{m2}");
    // Too narrow for that: compact cells, both shown, no %, none hidden.
    let s = screen(&mut v, 100, 40);
    let h = heads(&s);
    assert!(h.contains("C1") && h.contains("C2"), "{h}");
    let m2 = s
        .lines()
        .find(|l| l.contains("M2 ") && l.contains("genes"))
        .unwrap();
    assert!(!m2.contains(" 75"), "{m2}");
    assert!(!s.contains('‹') && !h.contains('›'), "{s}");
}

#[test]
fn hidden_cell_modules_are_counted_and_scroll_in_every_pane_together() {
    let m = loaded();
    let name = |k: usize| format!("C{}", m.columns[k] + 1);
    let mut v = ModulesView::new("test", &m);
    press(&mut v, KeyCode::Enter);
    press(&mut v, KeyCode::Char(' '));
    // One cell module fits: the other is counted, not dropped silently.
    let s = screen(&mut v, 60, 40);
    let h = heads(&s);
    assert!(
        h.contains("1›") && h.contains(&name(0)) && !h.contains(&name(1)),
        "{h}"
    );
    assert!(s.contains("[/]"), "the footer names the keys: {s}");
    press(&mut v, KeyCode::Char(']'));
    let s = screen(&mut v, 60, 40);
    let h = heads(&s);
    assert!(
        h.contains("‹1") && h.contains(&name(1)) && !h.contains(&name(0)),
        "{h}"
    );
    // The sites grid and the strip's tracks follow the same columns.
    let marked = s.lines().filter(|l| l.contains("‹1")).count();
    assert!(marked >= 3, "modules, sites and strip each say so: {s}");
    let track = |s: &str, c: &str| s.lines().any(|l| l.trim_start_matches('│').starts_with(c));
    assert!(track(&s, &name(1)) && !track(&s, &name(0)), "{s}");
    // No further than the last; back with [ or the wheel sideways.
    press(&mut v, KeyCode::Char(']'));
    assert_eq!(v.col_from, 1);
    press(&mut v, KeyCode::Char('['));
    assert_eq!(v.col_from, 0);
    mouse(
        &mut v,
        MouseEventKind::ScrollRight,
        (2, 5),
        KeyModifiers::NONE,
    );
    assert_eq!(v.col_from, 1);
    mouse(
        &mut v,
        MouseEventKind::ScrollUp,
        (2, 5),
        KeyModifiers::SHIFT,
    );
    assert_eq!(v.col_from, 0);
    // With room for all again, nothing is hidden.
    press(&mut v, KeyCode::Char(']'));
    let s = screen(&mut v, 180, 40);
    assert_eq!(v.col_from, 0, "{s}");
}

#[test]
fn sorting_by_a_hidden_column_scrolls_it_into_view() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    press(&mut v, KeyCode::Enter);
    screen(&mut v, 60, 40);
    press(&mut v, KeyCode::Char('o'));
    assert_eq!((v.sort, v.col_from), (Sort::Column(0), 0));
    press(&mut v, KeyCode::Char('o'));
    assert_eq!((v.sort, v.col_from), (Sort::Column(1), 1));
}

#[test]
fn a_short_terminal_draws_without_panicking() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    for h in [4, 8, 15] {
        screen(&mut v, 60, h);
    }
}

#[test]
fn pinned_genes_wrap_into_rows_in_a_narrow_strip() {
    let m = loaded();
    let mut v = ModulesView::new("test", &m);
    press(&mut v, KeyCode::Enter);
    press(&mut v, KeyCode::Char(' '));
    press(&mut v, KeyCode::Esc);
    press(&mut v, KeyCode::Down);
    press(&mut v, KeyCode::Enter);
    press(&mut v, KeyCode::Char(' '));
    press(&mut v, KeyCode::End);
    press(&mut v, KeyCode::Char(' '));
    assert_eq!(v.pins.len(), 3);
    // Room for two genes across: the third goes on a row below.
    let s = screen(&mut v, 36, 80);
    let row = |name: &str| s.lines().position(|l| l.contains(name)).unwrap();
    assert_eq!(row("GENE3 M2"), row("GENE1 M1"), "{s}");
    assert!(row("GENE2 M1") > row("GENE1 M1"), "{s}");
}
