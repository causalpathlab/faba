use super::*;
use ratatui::crossterm::event::KeyModifiers;

fn press(p: &mut impl View, code: KeyCode) {
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

/// A pileup picker at `cwd`, opening at `last` when given.
fn picker(cwd: PathBuf, last: Option<&Path>) -> InputPicker {
    InputPicker::new(cwd, last, "pileup", "site matrices", site_matrix, true)
}

fn text(p: &Path) -> Box<str> {
    p.to_string_lossy().into_owned().into_boxed_str()
}

#[test]
fn mark_files_from_a_folder_or_choose_the_folder() {
    let tmp = tree();
    let root = tmp.path().to_path_buf();
    let mut p = picker(root.clone(), None);
    assert!(p.browser.select("plain"));
    press(&mut p, KeyCode::Char(' '));
    assert!(
        p.decision.is_none() && p.error.is_some(),
        "not a faba folder"
    );

    // A faba folder, chosen whole.
    assert!(p.browser.select("out"));
    press(&mut p, KeyCode::Char(' '));
    let folder = Chosen {
        paths: vec![text(&root.join("out"))],
        separate: false,
    };
    assert_eq!(p.decision, Some(Some(folder)));

    // Or opened, and only its site matrices listed and marked.
    let mut p = picker(root.join("out"), None);
    let names: Vec<&str> = p.browser.shown().iter().map(|e| e.name()).collect();
    assert_eq!(names, ["..", "a_m6a_site.zarr.zip", "b_m6a_site.zarr.zip"]);
    for (name, key) in [
        ("b_m6a_site.zarr.zip", KeyCode::Char(' ')),
        ("a_m6a_site.zarr.zip", KeyCode::Char(' ')),
        ("b_m6a_site.zarr.zip", KeyCode::Char(' ')),
        ("b_m6a_site.zarr.zip", KeyCode::Char(' ')),
    ] {
        assert!(p.browser.select(name));
        press(&mut p, key);
    }
    let out = root.join("out");
    let marked = vec![
        out.join("a_m6a_site.zarr.zip"),
        out.join("b_m6a_site.zarr.zip"),
    ];
    assert_eq!(
        p.browser.marks.iter().collect::<Vec<_>>(),
        marked,
        "unmarked and marked again, in the order marked"
    );
    // Tab: a track each.
    press(&mut p, KeyCode::Tab);
    crate::tui::feed(&mut p, crate::tui::apply_key());
    let chosen = Chosen {
        paths: marked.iter().map(|m| text(m)).collect(),
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
        None,
        "metagene",
        "site tables",
        |n| n.ends_with("_sites.parquet"),
        false,
    );
    p.browser
        .marks
        .toggle(tmp.path().join("out/m6a_sites.parquet"));
    for (w, h) in [(20, 6), (80, 24), (200, 60)] {
        let mut t = ratatui::Terminal::new(ratatui::backend::TestBackend::new(w, h)).unwrap();
        t.draw(|f| p.render(f)).unwrap();
    }
}

#[test]
fn opens_at_the_last_choice() {
    let tmp = tree();
    let root = tmp.path().to_path_buf();
    let p = picker(root.clone(), Some(&root.join("out/b_m6a_site.zarr.zip")));
    assert_eq!(p.browser.dir, root.join("out"));
    assert_eq!(
        p.browser.current().map(|e| e.name()),
        Some("b_m6a_site.zarr.zip")
    );
    // Gone since: where it was asked from.
    let p = picker(root.clone(), Some(&root.join("gone/x_m6a_site.zarr.zip")));
    assert_eq!(p.browser.dir, root);
}

#[test]
fn the_last_choice_is_kept_per_working_folder() {
    let tmp = tempfile::tempdir().unwrap();
    let memo = tmp.path().join("cache/faba/last_input");
    let (a, b) = (Path::new("/work/a"), Path::new("/work/b"));
    assert_eq!(recall(&memo, a), None, "nothing kept yet");
    remember(&memo, a, "/work/a/out");
    remember(&memo, b, "/work/b/x_m6a_site.zarr.zip");
    remember(&memo, a, "/work/a/out2");
    assert_eq!(
        recall(&memo, a),
        Some(PathBuf::from("/work/a/out2")),
        "the latest"
    );
    assert_eq!(
        recall(&memo, b),
        Some(PathBuf::from("/work/b/x_m6a_site.zarr.zip")),
        "each folder its own"
    );
    let kept = std::fs::read_to_string(&memo).unwrap();
    assert_eq!(kept.lines().count(), 2, "one line per folder");
}

#[test]
fn the_genome_picker_lists_fasta_and_offers_a_download() {
    let tmp = tempfile::tempdir().unwrap();
    for f in ["hg38.fa", "mm10.fasta.gz", "notes.txt"] {
        std::fs::write(tmp.path().join(f), b"").unwrap();
    }
    let mut p = ReferencePicker::new(
        tmp.path().to_path_buf(),
        "pwm",
        Reference::Genome,
        Some(true),
    );
    let names: Vec<&str> = p.browser.shown().iter().map(|e| e.name()).collect();
    assert!(names.contains(&"hg38.fa") && names.contains(&"mm10.fasta.gz"));
    assert!(!names.contains(&"notes.txt"), "{names:?}");

    // Tab: the sources, on the one naming chromosomes chr1; Esc backs out.
    press(&mut p, KeyCode::Tab);
    assert!(matches!(&p.fetching, Fetching::Choosing(c) if c.bam_chr == Some(true)));
    press(&mut p, KeyCode::Esc);
    assert!(matches!(p.fetching, Fetching::No) && p.decision.is_none());

    p.browser.select("hg38.fa");
    press(&mut p, KeyCode::Enter);
    assert_eq!(p.decision, Some(Some(tmp.path().join("hg38.fa"))));

    let mut p = ReferencePicker::new(tmp.path().to_path_buf(), "pwm", Reference::Genome, None);
    press(&mut p, KeyCode::Esc);
    assert_eq!(p.decision, Some(None), "Esc cancels");
}

#[test]
fn the_annotation_picker_lists_gtf_and_gff() {
    let tmp = tempfile::tempdir().unwrap();
    for f in ["genes.gtf.gz", "genes.gff3", "hg38.fa"] {
        std::fs::write(tmp.path().join(f), b"").unwrap();
    }
    let mut p = ReferencePicker::new(
        tmp.path().to_path_buf(),
        "pileup",
        Reference::Annotation,
        None,
    );
    let names: Vec<&str> = p.browser.shown().iter().map(|e| e.name()).collect();
    assert!(names.contains(&"genes.gtf.gz") && names.contains(&"genes.gff3"));
    assert!(!names.contains(&"hg38.fa"), "{names:?}");
    press(&mut p, KeyCode::Tab);
    assert!(
        matches!(p.fetching, Fetching::Choosing(_)),
        "Tab downloads here too"
    );
}

#[test]
fn enter_takes_the_marked_or_the_file_here() {
    let tmp = tree();
    let out = tmp.path().join("out");
    // Nothing marked: the file under the cursor.
    let mut p = picker(out.clone(), None);
    assert!(p.browser.select("b_m6a_site.zarr.zip"));
    press(&mut p, KeyCode::Enter);
    let here = Chosen {
        paths: vec![text(&out.join("b_m6a_site.zarr.zip"))],
        separate: false,
    };
    assert_eq!(p.decision, Some(Some(here)));

    // Marked: those, not the one under the cursor.
    let mut p = picker(out.clone(), None);
    assert!(p.browser.select("a_m6a_site.zarr.zip"));
    press(&mut p, KeyCode::Char(' '));
    assert!(p.browser.select("b_m6a_site.zarr.zip"));
    press(&mut p, KeyCode::Enter);
    let marked = Chosen {
        paths: vec![text(&out.join("a_m6a_site.zarr.zip"))],
        separate: false,
    };
    assert_eq!(p.decision, Some(Some(marked)));

    // Esc with nothing typed cancels.
    let mut p = picker(out, None);
    press(&mut p, KeyCode::Esc);
    assert_eq!(p.decision, Some(None));
}

#[test]
fn a_folder_asks_before_dropping_marks() {
    let tmp = tree();
    let root = tmp.path().to_path_buf();
    let mut p = picker(root.join("out"), None);
    assert!(p.browser.select("a_m6a_site.zarr.zip"));
    press(&mut p, KeyCode::Char(' '));
    assert_eq!(p.browser.marks.len(), 1, "Space marks");
    // Up to the root, where `out` is a faba folder.
    press(&mut p, KeyCode::Left);
    assert!(p.browser.select("out"));
    press(&mut p, KeyCode::Char(' '));
    assert!(p.decision.is_none() && p.error.as_deref().is_some_and(|e| e.contains("again")));
    press(&mut p, KeyCode::Char(' '));
    assert!(p.decision.is_some(), "the second Space takes it");
}
