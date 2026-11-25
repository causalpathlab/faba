use super::super::figure_view;
use super::*;
use crate::site_analysis::pileup::tui::Track;
use crate::tui::ui::Scale;
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

/// `s` writes each gene's figure, and one of them all one under another.
#[test]
fn the_files_are_each_gene_and_the_genes_one_under_another() {
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    let files = s.files().unwrap();
    let names: Vec<&str> = files.iter().map(|f| f.0.as_str()).collect();
    assert_eq!(names, ["_GENE1", "_GENE2", "_column"]);
    // Each gene's on its own page, at the figure's own size.
    for (k, gene) in ["GENE1", "GENE2"].into_iter().enumerate() {
        let own = s.views[k].figure_on(&Sheet::OWN).unwrap();
        assert_eq!(
            crate::figure::svg_size(&files[k].1),
            crate::figure::svg_size(&own)
        );
        assert!(files[k].1.contains(&format!(">{gene}<")), "{gene}");
    }
    let height = |svg: &str| crate::figure::svg_size(svg).1;
    assert_eq!(
        height(&files[2].1),
        height(&files[0].1) + height(&files[1].1)
    );
    // They render, as a save does.
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().join("both");
    s.controls.save_many(&files, out.to_str().unwrap());
    for name in names {
        assert!(
            dir.path().join(format!("both{name}.pdf")).exists(),
            "{name}"
        );
    }
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
    let width = |svg: &str| crate::figure::svg_size(svg).0;
    let one = width(&s.views[0].figure_on(&Sheet::OWN).unwrap());
    let all = |s: &mut StackView| s.files().unwrap().pop().unwrap();
    assert_eq!(width(&all(&mut s).1), one, "one under another");
    press(&mut s, KeyCode::Char('|'));
    assert_eq!(s.arrange, Arrange::Row);
    let (name, svg) = all(&mut s);
    assert_eq!(
        (name.as_str(), width(&svg)),
        ("_row", 2.0 * one),
        "side by side"
    );
    // Side by side on screen, too: GENE2 to the right of GENE1, same rows.
    let text = screen(&mut s, 160, 30);
    let row_of = |i: usize| i / 160;
    let (g1, g2) = (text.find("GENE1").unwrap(), text.find("GENE2").unwrap());
    assert_eq!(row_of(g1), row_of(g2), "{text}");
    press(&mut s, KeyCode::Char('|'));
    assert_eq!(s.arrange, Arrange::Grid);
    assert_eq!(all(&mut s).0, "_grid");
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
fn zoom_keys_zoom_the_gene_with_the_keys_and_z_or_w_every_gene() {
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    screen(&mut s, 100, 40);
    let whole: Vec<_> = s.views.iter().map(|v| v.window).collect();
    for key in ['=', '+'] {
        press(&mut s, KeyCode::Char(key));
        assert_ne!(
            s.views[0].window, whole[0],
            "{key} zooms the gene with the keys"
        );
        assert_eq!(s.views[1].window, whole[1], "{key}: and only it");
        press(&mut s, KeyCode::Char('0'));
        assert_eq!(s.views[0].window, whole[0], "0: it whole");
    }
    press(&mut s, KeyCode::Char('Z'));
    assert!(
        s.views.iter().zip(&whole).all(|(v, w)| v.window != *w),
        "Z zooms every gene"
    );
    press(&mut s, KeyCode::Char('-'));
    assert_eq!(s.views[0].window, whole[0], "- the gene with the keys out");
    assert_ne!(s.views[1].window, whole[1], "and only it");
    // The layout-bound keys are gone.
    for gone in ['_', ')'] {
        press(&mut s, KeyCode::Char(gone));
        assert_ne!(s.views[1].window, whole[1], "{gone} does nothing");
    }
    press(&mut s, KeyCode::Char('W'));
    let now: Vec<_> = s.views.iter().map(|v| v.window).collect();
    assert_eq!(now, whole, "W shows every gene whole");
    assert!(s.exit.is_none());
}

/// Two genes, each with two alike read tracks with totals: the difference
/// and the mirror both there.
fn paired<'a>(m: &'a [(i64, f64)], half: &'a [(i64, f64)]) -> StackView<'a> {
    let one = |name: &str| {
        let mut a = Track::single("wt", "sum", m, false).with_total(half);
        a.mirror_with_next = true;
        let b = Track::single("mut", "sum", half, false).with_total(half);
        PileupView::new(name, "chr1", vec![a, b], EXTENT)
    };
    StackView::new(vec![one("GENE1"), one("GENE2")])
}

#[test]
fn the_footer_names_paging_the_difference_and_the_mirror_where_they_work() {
    let (m, half) = (positions(1.0), positions(2.0));
    let mut s = paired(&m, &half);
    let said = s.footer_line(false).to_string();
    for key in [
        "PgUp/PgDn half a view",
        "+/- zoom",
        "Z/W zoom all/all whole",
        "d difference/fold",
        "v split",
        "? keys",
    ] {
        assert!(said.contains(key), "{key}: {said}");
    }
    assert!(!said.contains("+/_") && !said.contains(") all"), "{said}");
    let (a, b) = (positions(1.0), positions(3.0));
    let plain = stack(&a, &b).footer_line(false).to_string();
    assert!(
        !plain.contains("d difference") && !plain.contains("v split"),
        "{plain}"
    );
    let keys = crate::tui::help::shown(&mut s, 120, 50).expect("? opens the keys");
    for key in ["Tab/Shift+Tab", "←/→ h/l", "PgUp/PgDn", "Z", "W", "|", "q"] {
        assert!(keys.contains(key), "{key}:\n{keys}");
    }
}

#[test]
fn a_drag_or_the_wheel_pans_the_gene_under_the_pointer() {
    use ratatui::crossterm::event::{MouseButton, MouseEvent, MouseEventKind};
    let (a, b) = (positions(1.0), positions(3.0));
    let mut s = stack(&a, &b);
    screen(&mut s, 100, 40);
    // Zoomed in, every gene, so there is room to pan.
    press(&mut s, KeyCode::Char('Z'));
    press(&mut s, KeyCode::Char('Z'));
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

/// `s` as the check draws it, the gene with the keys zoomed in on a site
/// mid-gene, so that every move, zoom and jump moves.
fn mid_gene(mut s: StackView) -> StackView {
    screen(&mut s, 160, 50);
    press(&mut s, KeyCode::Char('+'));
    press(&mut s, KeyCode::Char('n'));
    press(&mut s, KeyCode::Char('n'));
    s
}

/// `/`, `m` and `M` are listed too: they say the key needs one gene.
#[test]
fn the_keys_listed_are_the_keys_bound() {
    use crate::tui::help::assert_lists_its_keys;
    let (a, b) = (positions(1.0), positions(3.0));
    assert_lists_its_keys("two genes", || mid_gene(stack(&a, &b)));
    let half = positions(0.5);
    assert_lists_its_keys("compared and mirrored", || mid_gene(paired(&a, &half)));
    assert_lists_its_keys("rows hidden", || {
        let mut s = mid_gene(group_stack(&a, 0));
        press(&mut s, KeyCode::Char(']'));
        s
    });
}

/// Two genes' WT and MUT tracks over their totals, and their gene models.
fn marked<'a>(m: &'a [(i64, f64)], half: &'a [(i64, f64)]) -> Vec<(View<'a>, Vec<Track<'a>>)> {
    ["ID1_GENE1", "ID2_GENE2"]
        .into_iter()
        .map(|title| {
            let mut a = Track::single("WT", "sum", m, false).with_total(half);
            a.mirror_with_next = true;
            let b = Track::single("MUT", "sum", half, false).with_total(half);
            let view = View {
                title,
                chr: "chr1",
                extent: EXTENT,
                on: "methylated",
                off: "unmethylated",
                keys: None,
                genes: None,
                cursor: None,
                switchable: false,
                from_cluster: false,
                track_from: 0,
            };
            (view, vec![a, b])
        })
        .collect()
}

/// The genes marked in the list and saved from their stack are the files
/// `faba pileup --batch-process -q GENE1,GENE2` writes, drawn the same way
/// (its bars, y scale, difference row and key as the flags say): the same
/// names and the same figures, the same PDFs.
#[test]
fn a_stack_saves_what_the_command_line_writes_for_the_same_genes() {
    use crate::site_analysis::pileup::PileupArgs;
    #[derive(clap::Parser)]
    struct Cli {
        #[command(flatten)]
        args: PileupArgs,
    }
    // Sites in the middle of the genes, so a window at them is not whole.
    let middle = |scale| -> Vec<(i64, f64)> {
        let p = positions(scale).into_iter();
        p.filter(|p| (1_030_000..1_060_000).contains(&p.0))
            .collect()
    };
    let (m, half) = (middle(1.0), middle(2.0));
    let views = marked(&m, &half)
        .into_iter()
        .map(|(view, tracks)| {
            let mut b = super::browser(view, tracks, None);
            b.show = crate::site_analysis::show::Show::Both;
            b
        })
        .collect();
    let mut s = StackView::new(views);
    screen(&mut s, 120, 40);
    // Drawn as sqrt, the difference row hidden, the key beside, a row,
    // every gene at its sites (the command line's window).
    for key in ['y', 'y', 'D', 'K', '|', 'W'] {
        press(&mut s, KeyCode::Char(key));
    }
    screen(&mut s, 120, 40);
    let columns = s.views[0].columns;
    assert!(s.views.iter().all(|v| v.columns == columns));
    let tui = s.files().unwrap();
    let bins = columns.to_string();
    let line = [
        "pileup",
        "x.zarr",
        "-q",
        "GENE1,GENE2",
        "--yscale",
        "sqrt",
        "--key-at",
        "right",
        "--row",
        "--bins",
        &bins,
    ];
    let args = <Cli as clap::Parser>::try_parse_from(line).unwrap().args;
    let settings = args.figure_settings();
    let mut cli: Vec<PileupView> = marked(&m, &half)
        .into_iter()
        .map(|(view, tracks)| {
            let at = (args.figure_window(), columns);
            figure_view(view, tracks, &settings, at, args.key)
        })
        .collect();
    let written = export(&mut cli, &Sheet::OWN, args.combine(2).unwrap()).unwrap();
    assert_eq!(tui.len(), 3);
    assert!(cli.iter().all(|v| v.window != v.extent), "at the sites");
    assert_eq!(written, tui);
    assert!(!tui[0].1.contains("WT vs MUT"), "no difference row");
    // The PDFs too, as the browser and the command line each write them.
    let dir = tempfile::tempdir().unwrap();
    let at = |name: &str| dir.path().join(name).to_str().unwrap().to_string();
    s.controls.save_many(&tui, &at("tui"));
    for (suffix, svg) in &written {
        crate::figure::save(svg, &at(&format!("cli{suffix}"))).unwrap();
        // svg2pdf numbers a page's fonts in an order that differs from one
        // run to the next (two runs of the command line differ so too): the
        // same size, the same fonts.
        let read = |p: String| {
            let pdf = std::fs::read(format!("{p}.pdf")).unwrap();
            let text = String::from_utf8_lossy(&pdf).to_string();
            let mut fonts: Vec<String> = (text.split("/BaseFont /").skip(1))
                .filter_map(|f| {
                    f.split_whitespace()
                        .next()?
                        .split('+')
                        .nth(1)
                        .map(String::from)
                })
                .collect();
            fonts.sort();
            (pdf.len(), fonts)
        };
        let (a, b) = (
            read(at(&format!("tui{suffix}"))),
            read(at(&format!("cli{suffix}"))),
        );
        assert!(!a.1.is_empty(), "{suffix}: no fonts");
        assert_eq!(a, b, "{suffix}: the PDFs differ");
    }
}

/// `W` puts every gene at its own sites from whole, as `0` does one, and
/// every gene whole again: each by its own sites, not the first's.
#[test]
fn w_steps_every_gene_between_whole_and_its_own_sites() {
    let a: Vec<(i64, f64)> = (0..100).map(|i| (1_020_000 + i * 100, 1.0)).collect();
    let b: Vec<(i64, f64)> = (0..100).map(|i| (1_060_000 + i * 200, 1.0)).collect();
    let mut s = stack(&a, &b);
    screen(&mut s, 120, 40);
    press(&mut s, KeyCode::Char('W'));
    let sites: Vec<_> = s.views.iter().map(|v| v.sites_window()).collect();
    let now: Vec<_> = s.views.iter().map(|v| v.window).collect();
    assert_eq!(now, sites);
    assert_ne!(sites[0], sites[1]);
    assert!(
        sites[0].0 < 1_020_000 && 1_029_900 < sites[0].1,
        "{sites:?}"
    );
    press(&mut s, KeyCode::Char('W'));
    assert!(s.views.iter().all(|v| v.window == EXTENT));
}
