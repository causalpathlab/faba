use super::*;
use ratatui::crossterm::event::KeyModifiers;

fn press(p: &mut InputPicker, code: KeyCode) {
    crate::tui::feed(p, KeyEvent::new(code, KeyModifiers::NONE));
}

/// `root/{plain/, out/{a_m6a_site.zarr.zip, b_m6a_site.zarr.zip, m6a_sites.parquet}}`.
fn tree() -> tempfile::TempDir {
    let tmp = tempfile::tempdir().unwrap();
    let out = tmp.path().join("out");
    std::fs::create_dir_all(tmp.path().join("plain")).unwrap();
    std::fs::create_dir_all(&out).unwrap();
    for f in [
        "a_m6a_site.zarr.zip",
        "b_m6a_site.zarr.zip",
        "m6a_sites.parquet",
    ] {
        std::fs::write(out.join(f), b"").unwrap();
    }
    tmp
}

fn site_matrix(name: &str) -> bool {
    crate::qc::layout::is_site_matrix_name(name)
}

#[test]
fn mark_files_from_a_folder_or_choose_the_folder() {
    let tmp = tree();
    let root = tmp.path().to_path_buf();
    let mut p = InputPicker::new(root.clone(), "pileup", "site matrices", site_matrix, true);
    assert!(p.browser.list.select("plain"));
    press(&mut p, KeyCode::Char(' '));
    assert!(
        p.decision.is_none() && p.error.is_some(),
        "not a faba folder"
    );

    // A faba folder, chosen whole.
    assert!(p.browser.list.select("out"));
    press(&mut p, KeyCode::Char(' '));
    let folder = Chosen {
        paths: vec![root.join("out")],
        separate: false,
    };
    assert_eq!(p.decision, Some(Some(folder)));

    // Or opened, and only its site matrices listed and marked.
    let mut p = InputPicker::new(
        root.join("out"),
        "pileup",
        "site matrices",
        site_matrix,
        true,
    );
    let names: Vec<&str> = p.browser.list.shown().map(|e| e.name.as_str()).collect();
    assert_eq!(names, ["..", "a_m6a_site.zarr.zip", "b_m6a_site.zarr.zip"]);
    for (name, key) in [
        ("b_m6a_site.zarr.zip", KeyCode::Char(' ')),
        ("a_m6a_site.zarr.zip", KeyCode::Enter),
        ("b_m6a_site.zarr.zip", KeyCode::Char(' ')),
        ("b_m6a_site.zarr.zip", KeyCode::Char(' ')),
    ] {
        assert!(p.browser.list.select(name));
        press(&mut p, key);
    }
    let out = root.join("out");
    let marked = vec![
        out.join("a_m6a_site.zarr.zip"),
        out.join("b_m6a_site.zarr.zip"),
    ];
    assert_eq!(
        p.browser.marked.as_ref(),
        Some(&marked),
        "unmarked and marked again, in the order marked"
    );
    // Tab: a track each.
    press(&mut p, KeyCode::Tab);
    crate::tui::feed(&mut p, crate::tui::apply_key());
    let chosen = Chosen {
        paths: marked,
        separate: true,
    };
    assert_eq!(
        p.decision,
        Some(Some(chosen)),
        "the marked, apart, not the folder"
    );
}

#[test]
fn draws_at_any_terminal_size() {
    let tmp = tree();
    let mut p = InputPicker::new(
        tmp.path().join("out"),
        "metagene",
        "site tables",
        |n| n.ends_with("_sites.parquet"),
        false,
    );
    p.browser
        .toggle_mark(tmp.path().join("out/m6a_sites.parquet"));
    for (w, h) in [(20, 6), (80, 24), (200, 60)] {
        let mut t = ratatui::Terminal::new(ratatui::backend::TestBackend::new(w, h)).unwrap();
        t.draw(|f| p.render(f)).unwrap();
    }
}

#[test]
fn opens_at_the_last_choice() {
    let tmp = tree();
    let root = tmp.path().to_path_buf();
    let mut p = InputPicker::new(root.clone(), "pileup", "site matrices", site_matrix, true);
    p.start_at(&root.join("out/b_m6a_site.zarr.zip"));
    assert_eq!(p.browser.cwd, root.join("out"));
    assert_eq!(
        p.browser.highlighted().map(|e| e.name.as_str()),
        Some("b_m6a_site.zarr.zip")
    );
    // Gone since: where it was.
    let mut p = InputPicker::new(root.clone(), "pileup", "site matrices", site_matrix, true);
    p.start_at(&root.join("gone/x_m6a_site.zarr.zip"));
    assert_eq!(p.browser.cwd, root);
}
