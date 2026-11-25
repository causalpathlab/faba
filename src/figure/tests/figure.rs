use super::*;
use ratatui::crossterm::event::KeyModifiers;

fn key(code: KeyCode) -> KeyEvent {
    KeyEvent::new(code, KeyModifiers::NONE)
}

fn panel() -> String {
    let values = [0.0, 3.0, 10.0, 4.0, 1.0];
    let front = [0.0, 1.0, 8.0, 4.0, 0.0];
    let mut c = Canvas::new(300.0, 160.0);
    Bars {
        rotated: false,
        values: &values,
        front: Some(&front),
        accent: &|i| i < 2,
        colour: &|i| if i == 3 { "#0072b2" } else { BAR },
        y_scale: Scale::Log,
        y_max: None,
        ticks: vec![(0, "0".into()), (4, "L4 & <x>".into())],
        pointer: Some(2),
        marks: vec![1, 3],
        dividers: vec![2],
        title: "PANEL1".into(),
        x_title: "value".into(),
        y_title: "n".into(),
    }
    .draw(&mut c, 0.0, 0.0, 300.0, 160.0);
    c.glyph(10.0, 10.0, 20.0, 40.0, 'A', ACCENT);
    c.finish()
}

#[test]
fn a_divider_sits_on_the_bar_edge_and_drops_tick_stubs() {
    // Five bars over a 236-wide plot from x = 52: bar 2's left edge is 146.4.
    let svg = panel();
    assert!(
        svg.contains(r#"x1="146.40" y1="20.00" x2="146.40""#),
        "{svg}"
    );
    // With dividers, tick labels name spans: no solid stub under them. The
    // stub of the tick at bar 4 would sit at its centre, x = 264.4.
    assert!(!svg.contains(r#"x1="264.40""#));
}

#[test]
fn a_bar_takes_its_own_colour_unless_accented() {
    let svg = panel();
    // Bar 3 is not accented, so it takes its colour; bars 0 and 1 stay accent.
    assert!(svg.contains(r##"fill="#0072b2""##), "{svg}");
    assert!(svg.contains(&format!(r#"fill="{ACCENT}""#)));
}

#[test]
fn svg_is_well_formed_and_escaped() {
    let svg = panel();
    assert!(svg.starts_with("<svg"));
    assert!(svg.contains("L4 &amp; &lt;x&gt;"));
    assert!(svg.contains(ACCENT));
    let tree = usvg::Tree::from_str(&svg, &usvg::Options::default()).unwrap();
    assert_eq!(tree.size().width(), 300.0);
}

#[test]
fn save_writes_a_pdf_only() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("fig1.pdf");
    let paths = save(&panel(), prefix.to_str().unwrap()).unwrap();
    assert_eq!(paths.len(), 1);
    let pdf = std::fs::read(&paths[0]).unwrap();
    assert!(pdf.starts_with(b"%PDF"), "{}", paths[0]);
    assert!(paths[0].ends_with("fig1.pdf"));
    assert!(!dir.path().join("fig1.png").exists(), "no PNG");
}

#[test]
fn prompt_confirms_edits_and_cancels() {
    let mut p = SavePrompt::new("view1");
    assert!(!p.active() && p.footer().is_none());
    p.open();
    assert!(p.active());
    for _ in 0..5 {
        p.handle(key(KeyCode::Backspace));
    }
    for ch in "out/a".chars() {
        assert_eq!(p.handle(key(KeyCode::Char(ch))), None);
    }
    assert_eq!(p.handle(key(KeyCode::Enter)).as_deref(), Some("out/a"));
    assert!(!p.active());
    p.report(Ok(vec!["out/a.pdf".into()]));
    assert!(p.footer().is_some());
    p.dismiss();
    assert!(p.footer().is_none());

    p.open();
    assert!(p.footer().is_some(), "the last name is offered again");
    p.handle(key(KeyCode::Esc));
    assert!(!p.active());
}

#[test]
fn a_page_keeps_the_views_shape_when_a_minimum_binds() {
    let (cw, ch) = CELL;
    let area = Rect::new(0, 0, 10, 10);
    assert_eq!(page(area, (0.0, 0.0), CELL), Some((10.0 * cw, 10.0 * ch)));
    // Too narrow for the least width: both sides grow alike.
    let (w, h) = page(area, (6.0 * 10.0 * cw, 0.0), CELL).unwrap();
    assert_eq!((w, h), (6.0 * 10.0 * cw, 6.0 * 10.0 * ch));
    // Too short: the same.
    let (w, h) = page(area, (0.0, 2.0 * 10.0 * ch), CELL).unwrap();
    assert_eq!((w, h), (2.0 * 10.0 * cw, 2.0 * 10.0 * ch));
    assert_eq!(page(Rect::default(), (1.0, 1.0), CELL), None);
}

#[test]
fn a_page_takes_its_cells_shape_from_the_terminals_font() {
    let area = Rect::new(0, 0, 10, 10);
    let text = Controls::new("x");
    assert_eq!(text.page(area, (0.0, 0.0)), page(area, (0.0, 0.0), CELL));
    // A font 8 by 20 pixels: cells 0.4 as wide as high, as tall as CELL.
    let (cw, ch) = font_cell(ratatui_image::FontSize {
        width: 8,
        height: 20,
    });
    assert!((cw / ch - 0.4).abs() < 1e-9 && ch == CELL.1, "{cw}x{ch}");
    assert_eq!(
        font_cell(ratatui_image::FontSize {
            width: 0,
            height: 0
        }),
        CELL
    );
    // With a picker, its font's cells.
    let picker = ratatui_image::picker::Picker::halfblocks();
    let cell = font_cell(picker.font_size());
    let images = Controls::new("x").with_picker(picker);
    assert_eq!(images.page(area, (0.0, 0.0)), page(area, (0.0, 0.0), cell));
}

#[test]
fn empty_shapes_are_not_written() {
    let mut c = Canvas::new(100.0, 200.0);
    c.rect(1.0, 1.0, 10.0, 0.0, INK);
    c.rect(1.0, 1.0, 0.0, 10.0, INK);
    c.rect(1.0, 1.0, 10.0, f64::NAN, INK);
    let values = [Some(0.0), Some(1.5), None, Some(-0.5), Some(0.0)];
    Diverging {
        rotated: false,
        values: &values,
        ticks: Vec::new(),
        pointer: None,
        title: String::new(),
        x_title: String::new(),
        y_title: String::new(),
        label: &|v| format!("{v:+.1}"),
    }
    .draw(&mut c, 0.0, 0.0, 100.0, 200.0);
    let svg = c.finish();
    assert!(!svg.contains(r#"height="0.00""#), "{svg}");
    assert!(!svg.contains("NaN"), "{svg}");
    // Two non-zero bars plus the background.
    assert_eq!(svg.matches("<rect").count(), 3, "{svg}");
}

#[test]
fn the_save_prompt_pops_up_with_the_figures_saved_here() {
    use ratatui::backend::TestBackend;
    use ratatui::Terminal;
    let draw = |c: &mut Controls| {
        let mut term = Terminal::new(TestBackend::new(100, 30)).unwrap();
        term.draw(|f| c.render_save(f, f.area())).unwrap();
        let buf = term.backend().buffer().clone();
        (0..30)
            .map(|y| (0..100).map(|x| buf[(x, y)].symbol()).collect::<String>())
            .collect::<Vec<_>>()
            .join("\n")
    };
    let mut c = Controls::new("pileup_GENE1");
    assert!(!draw(&mut c).contains("save the figure"), "closed, nothing");
    assert_eq!(c.key(key(KeyCode::Char('s'))), Key::Used);
    let text = draw(&mut c);
    assert!(
        text.contains("save the figure") && text.contains("pileup_GENE1"),
        "{text}"
    );

    // Saved, it is listed the next time the prompt opens.
    let tmp = tempfile::tempdir().unwrap();
    let prefix = tmp.path().join("fig1").to_string_lossy().into_owned();
    c.key(key(KeyCode::Esc));
    c.save(&panel(), &prefix);
    c.key(key(KeyCode::Char('s')));
    let text = draw(&mut c);
    assert!(
        text.contains("fig1.pdf") && text.contains("pileup"),
        "{text}"
    );
}

#[test]
fn the_prompt_never_replaces_a_saved_figure() {
    let dir = tempfile::tempdir().unwrap();
    let at = |n: &str| dir.path().join(n).to_string_lossy().into_owned();
    std::fs::write(at("logo.pdf"), b"").unwrap();
    std::fs::write(at("logo_2.png"), b"").unwrap();
    let mut p = SavePrompt::new(&at("logo"));
    p.open();
    assert_eq!(
        p.text(),
        Some(at("logo_3").as_str()),
        "the first name free of both"
    );

    // A taken name typed in keeps the prompt open, with the next free one.
    p.handle(key(KeyCode::Esc));
    p.open();
    for _ in 0..2 {
        p.handle(key(KeyCode::Backspace));
    }
    assert_eq!(p.text(), Some(at("logo").as_str()));
    assert_eq!(p.handle(key(KeyCode::Enter)), None, "logo.pdf is there");
    assert!(p.active() && p.taken().is_some());
    assert_eq!(p.text(), Some(at("logo_3").as_str()));
    assert_eq!(p.handle(key(KeyCode::Enter)), Some(at("logo_3")));
    assert!(p.taken().is_none());
}

#[test]
fn a_bar_too_thin_to_write_is_skipped() {
    let mut c = Canvas::new(100.0, 100.0);
    c.rect(1.0, 1.0, 10.0, 0.004, INK);
    c.rect(1.0, 1.0, 0.001, 10.0, INK);
    c.rect(1.0, 1.0, 10.0, f64::INFINITY, INK);
    c.rect(1.0, 1.0, 10.0, 0.5, INK);
    let svg = c.finish();
    assert!(
        !svg.contains(r#"height="0.00""#) && !svg.contains(r#"width="0.00""#),
        "{svg}"
    );
    assert!(!svg.contains("inf"), "{svg}");
    // The background, and the one bar thick enough to see.
    assert_eq!(svg.matches("<rect").count(), 2, "{svg}");
}

#[test]
fn figures_tile_in_rows_of_the_given_width() {
    let fig = |w: f64, h: f64| Canvas::new(w, h).finish();
    let figs = [fig(100.0, 50.0), fig(80.0, 70.0), fig(100.0, 40.0)];
    let size = |svg: &str| {
        let get = |k: &str| -> f64 {
            let at = svg.find(&format!(" {k}=\"")).unwrap() + k.len() + 3;
            svg[at..].split('"').next().unwrap().parse().unwrap()
        };
        (get("width"), get("height"))
    };
    assert_eq!(size(&tile_svgs(&figs, 1)), (100.0, 160.0), "a column");
    assert_eq!(size(&tile_svgs(&figs, 3)), (300.0, 70.0), "a row");
    assert_eq!(size(&tile_svgs(&figs, 2)), (200.0, 110.0), "two to a row");
    assert!(tile_svgs(&figs, 2).contains(r#"<svg x="100" y="0" "#));
}

#[test]
fn panels_fit_at_their_least_and_share_a_taller_page_in_proportion() {
    let panels = Panels {
        fixed: 60.0,
        least: &[280.0, 210.0],
    };
    assert_eq!(panels.height(), 550.0);
    assert_eq!(panels.heights(550.0), vec![280.0, 210.0]);
    let tall = panels.heights(1040.0 + 60.0);
    assert!((tall.iter().sum::<f64>() - 1040.0).abs() < 1e-9, "{tall:?}");
    assert!((tall[0] / tall[1] - 280.0 / 210.0).abs() < 1e-9, "{tall:?}");
    // Equal panels share alike.
    let even = Panels {
        fixed: 50.0,
        least: &[100.0; 3],
    };
    assert_eq!(even.heights(350.0), vec![100.0; 3]);
    assert_eq!(even.heights(650.0), vec![200.0; 3]);
}

#[test]
fn images_are_asked_of_the_terminal_attached_once_per_attach() {
    use ratatui_image::picker::Picker;
    // Outside a session (no terminal to ask) a detecting view draws text.
    let mut c = Controls::new("x").detect();
    c.refresh();
    assert!(c.images().is_none());
    assert_eq!(c.picked_at, Some(legume_tui::attaches()));
    // Another attach: asked again, the thumbnails made again.
    c.thumbs.insert("thumb.png".into(), None);
    c.picked_at = Some(u64::MAX);
    c.refresh();
    assert!(c.thumbs.is_empty());
    assert_eq!(c.picked_at, Some(legume_tui::attaches()));
    // The same attach: nothing asked.
    c.thumbs.insert("thumb.png".into(), None);
    c.refresh();
    assert_eq!(c.thumbs.len(), 1);
    // A picker handed in is kept.
    let mut given = Controls::new("x").with_picker(Picker::halfblocks());
    given.refresh();
    assert!(given.images().is_some());
}

#[test]
fn turning_images_off_outlives_another_terminal() {
    use ratatui_image::picker::Picker;
    let mut c = Controls::new("x");
    assert!(!c.images_on(), "no picker: text");
    c.picker = Some(Picker::halfblocks());
    assert!(c.images_on(), "a terminal with images: pictures");
    assert_eq!(c.key(key(KeyCode::Char('i'))), Key::Used);
    assert!(!c.images_on());
    c.picker = Some(Picker::halfblocks());
    assert!(!c.images_on(), "off stays off for the next terminal");
    c.set_images(true);
    assert!(c.images_on());
}

/// `--figure` writes the PDF and a PNG at two pixels per point beside it,
/// under any spelling of the prefix.
#[test]
fn save_with_png_writes_the_pdf_and_png() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("fig").to_string_lossy().into_owned();
    let svg = svg(100.0, 50.0, |c| c.rect(10.0, 10.0, 20.0, 20.0, INK));
    let paths = save_with_png(&svg, &format!("{prefix}.pdf")).unwrap();
    let want = ["pdf", "png"].map(|x| format!("{prefix}.{x}"));
    assert_eq!(paths, want);
    assert_eq!(saved_paths(&prefix), want, "what a save is checked against");
    let png = std::fs::read(&paths[1]).unwrap();
    assert!(png.starts_with(b"\x89PNG"));
    let image = image::load_from_memory(&png).unwrap();
    assert_eq!((image.width(), image.height()), (200, 100));
}

/// The PDF keeps a label as text in a font, not only as outlines, so it
/// stays editable.
#[test]
fn a_pdf_keeps_its_labels_as_text() {
    let dir = tempfile::tempdir().unwrap();
    let prefix = dir.path().join("fig").to_string_lossy().into_owned();
    let svg = svg(100.0, 50.0, |c| {
        c.text(10.0, 30.0, "LABEL1", 8.0, Anchor::Start, INK)
    });
    let pdf = std::fs::read(&save(&svg, &prefix).unwrap()[0]).unwrap();
    let has = |what: &[u8]| pdf.windows(what.len()).any(|w| w == what);
    assert!(
        has(b"/Font") && has(b"/FontDescriptor"),
        "no font in the PDF"
    );
}

/// A canvas styled for a font draws at its 8-point layout on a page
/// shrunk or grown by the font, inside one scale: its labels print scaled
/// from the 8 points they were laid out for, or every one at the size asked
/// for; its strokes scale with them, its outer size stays as asked, and a
/// sequence-logo letter keeps filling its box.
#[test]
fn a_styled_canvas_draws_its_layout_at_one_scale() {
    let draw = |style: &style::Style| {
        let mut c = Canvas::styled(100.0, 50.0, style);
        let page = c.page();
        c.text(1.0, 10.0, "T1", 8.0, Anchor::Start, INK);
        c.bold(1.0, 20.0, "T2", 11.0, Anchor::Start, INK);
        c.vtext(1.0, 30.0, "T3", 9.0, INK);
        c.rotated(1.0, 40.0, "T4", 8.0, Anchor::End, INK);
        c.line(0.0, 0.0, 10.0, 0.0, INK, 0.5);
        c.glyph(50.0, 10.0, 20.0, 30.0, 'A', INK);
        (page, c.finish())
    };
    let fit = style::Style {
        font: 6.0,
        ..Default::default()
    };
    let (page, scaled) = draw(&fit);
    assert_eq!(page, (100.0 / 0.75, 50.0 / 0.75));
    assert_eq!(svg_size(&scaled), (100.0, 50.0));
    assert_eq!(scaled.matches("<g ").count(), 1, "{scaled}");
    assert!(
        scaled.contains(r#"<g transform="scale(0.75)">"#),
        "{scaled}"
    );
    // The background fills the page as asked, outside the scale.
    let ground = r##"<rect x="0.00" y="0.00" width="100.00" height="50.00" fill="#ffffff"/>"##;
    assert!(scaled.contains(ground), "{scaled}");
    // Written at the layout's sizes and places, printed scaled.
    let bold = r#"font-size="11" font-weight="bold""#;
    assert!(scaled.contains(bold), "{scaled}");
    assert!(scaled.contains(r#"stroke-width="0.5""#), "{scaled}");
    let printed = |svg: &str| -> Vec<String> { labels(svg).into_iter().map(|l| l.0).collect() };
    assert_eq!(printed(&scaled), ["6", "8.25", "6.75", "6"], "{scaled}");
    assert_eq!(labels(&scaled)[0].2, Some(7.5), "a baseline as printed");
    let asked = style::Style {
        uniform: true,
        ..fit
    };
    let (_, uniform) = draw(&asked);
    assert!(printed(&uniform).iter().all(|l| l == "6"), "{uniform}");
    for svg in [&scaled, &uniform] {
        assert!(
            svg.contains(r#"font-size="100" font-weight="bold""#),
            "{svg}"
        );
    }
    let (page, own) = draw(&style::Style::default());
    assert_eq!(page, (100.0, 50.0));
    assert!(!own.contains("<g "), "no scale at 8 points: {own}");
    assert_eq!(own, {
        let mut c = Canvas::new(100.0, 50.0);
        c.text(1.0, 10.0, "T1", 8.0, Anchor::Start, INK);
        c.bold(1.0, 20.0, "T2", 11.0, Anchor::Start, INK);
        c.vtext(1.0, 30.0, "T3", 9.0, INK);
        c.rotated(1.0, 40.0, "T4", 8.0, Anchor::End, INK);
        c.line(0.0, 0.0, 10.0, 0.0, INK, 0.5);
        c.glyph(50.0, 10.0, 20.0, 30.0, 'A', INK);
        c.finish()
    });
}

/// Bars lay out at 8 points on the canvas' page: at 6 points the page is a
/// third larger and their margins are those at 8, all scaled by three
/// quarters; the panel alone draws no title.
#[test]
fn bars_lay_out_on_the_canvas_page() {
    let values = [1.0, 3.0, 2.0];
    let draw = |style: &style::Style| {
        let mut c = Canvas::styled(200.0, 100.0, style);
        let (w, h) = c.page();
        Bars {
            rotated: false,
            values: &values,
            front: None,
            accent: &|i| i == 2,
            colour: &|_| BAR,
            y_scale: Scale::Linear,
            y_max: None,
            ticks: vec![(0, "0".into())],
            pointer: Some(1),
            marks: Vec::new(),
            dividers: Vec::new(),
            title: "PANEL1".into(),
            x_title: "value".into(),
            y_title: "n".into(),
        }
        .draw(&mut c, 0.0, 0.0, w, h);
        c.finish()
    };
    let small = style::Style {
        font: 6.0,
        ..Default::default()
    };
    let svg = draw(&small);
    // The y axis 52 from the left, from 20 down to 38 above the bottom of
    // a page 133.33 high; scaled by 0.75.
    assert!(svg.contains(r#"<g transform="scale(0.75)">"#), "{svg}");
    assert!(
        svg.contains(r#"x1="52.00" y1="20.00" x2="52.00" y2="95.33""#),
        "{svg}"
    );
    assert!(svg.contains("PANEL1"));
    let panel = style::Style {
        panel_only: true,
        palette: style::Palette::Gray,
        ..small
    };
    let svg = draw(&panel);
    assert!(!svg.contains("PANEL1"), "{svg}");
    // Room for the top label only, not a title.
    assert!(svg.contains(r#"x1="52.00" y1="8.00" x2="52.00""#), "{svg}");
    // Gray: what a threshold drops in light gray, the pointer in ink.
    let removed = style::Palette::Gray.roles().removed;
    assert!(svg.contains(&format!(r#"fill="{removed}""#)), "{svg}");
    assert!(!svg.contains(ACCENT), "{svg}");
}

/// A figure drawn at its own size, with no font asked for, is drawn as it
/// always was.
#[test]
fn the_default_figures_are_unchanged() {
    let mut c = Canvas::new(300.0, 200.0);
    let up = [1.0, 4.0, 2.0];
    let down = [0.5, 0.0, 3.0];
    Mirror {
        up: Half {
            values: &up,
            front: Some(&down),
            colour: BAR,
            name: "up".into(),
        },
        down: Half::plain(&down, BAR),
        y_scale: Scale::Linear,
        y_max: None,
        y_labels: Default::default(),
        ticks: vec![(0, "L1-0".into()), (2, "L1-2".into())],
        pointer: Some(1),
        title: "PANEL2".into(),
        x_title: "x".into(),
        y_title: "y".into(),
        rotated: true,
    }
    .draw(&mut c, 0.0, 0.0, 300.0, 200.0);
    let shown = [fingerprint(&panel()), fingerprint(&c.finish())];
    assert_eq!(shown, [3060415989485868728, 16840006330901660691]);
}

/// A title too long for its room drops its closing parenthetical, then
/// what follows a colon, then is cut short with an ellipsis; one that fits
/// is as it is.
#[test]
fn a_long_title_is_shortened_to_fit() {
    let at = |s: &str, room: f64| super::fit_text(s, room, 8.0);
    let w = |s: &str| super::text_width(s, 8.0);
    let s = "TITLE1 TITLE2 (TITLE3)";
    assert_eq!(at(s, w(s)), s);
    assert_eq!(at(s, w(s) - 1.0), "TITLE1 TITLE2");
    let s = "TITLE1: TITLE2 TITLE3 (TITLE4)";
    assert_eq!(at(s, w("TITLE1: TITLE2")), "TITLE1");
    let cut = at("TITLE1TITLE2TITLE3", w("TITLE1TI"));
    assert!(cut.ends_with('…') && w(&cut) <= w("TITLE1TI"), "{cut}");
}
