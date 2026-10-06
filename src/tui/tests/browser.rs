use super::*;
use ratatui::crossterm::event::KeyModifiers;

fn key(code: KeyCode) -> KeyEvent {
    KeyEvent::new(code, KeyModifiers::NONE)
}

fn typed(b: &mut Browser, text: &str) {
    for c in text.chars() {
        assert!(matches!(b.key(key(KeyCode::Char(c))), Nav::Moved), "{c}");
    }
}

fn names(b: &Browser) -> Vec<&str> {
    b.list.shown().map(|e| e.name.as_str()).collect()
}

/// `root/{ann/{genes.gtf.gz, notes.txt}, b.gff3, a.txt, .hidden.gtf, m.zarr/}`.
fn tree() -> tempfile::TempDir {
    let tmp = tempfile::tempdir().unwrap();
    let root = tmp.path();
    std::fs::create_dir_all(root.join("ann")).unwrap();
    std::fs::create_dir_all(root.join("m.zarr")).unwrap();
    for f in [
        "ann/genes.gtf.gz",
        "ann/notes.txt",
        "b.gff3",
        "a.txt",
        ".hidden.gtf",
    ] {
        std::fs::write(root.join(f), b"").unwrap();
    }
    tmp
}

fn gff(name: &str) -> bool {
    name.ends_with(".gtf.gz") || name.ends_with(".gff3")
}

fn gff_browser(root: &Path) -> Browser {
    let tag = |p: &Path| p.extension().is_some_and(|e| e == "gff3");
    Browser::new(root.to_path_buf(), gff, tag).opened()
}

#[test]
fn lists_folders_then_the_kept_files_and_picks_a_file() {
    let tmp = tree();
    let root = normalize(tmp.path());
    let mut b = gff_browser(&root);
    // Hidden files and zarr stores are left out; other files unless kept.
    assert_eq!(names(&b), ["..", "ann", "b.gff3"]);
    assert!(b.list.shown().nth(2).unwrap().tagged && !b.list.shown().nth(1).unwrap().tagged);

    // Into `ann`, where only the annotation is listed, and pick it.
    b.list.at = 1;
    assert!(matches!(b.key(key(KeyCode::Enter)), Nav::Moved));
    assert_eq!(b.cwd, root.join("ann"));
    assert_eq!(names(&b), ["..", "genes.gtf.gz"]);
    b.list.at = 1;
    // Right only opens folders; Enter picks a file.
    assert!(matches!(b.key(key(KeyCode::Right)), Nav::Moved));
    assert!(
        matches!(b.key(key(KeyCode::Enter)), Nav::Picked(p) if p == root.join("ann/genes.gtf.gz"))
    );

    // Up lands back on the folder just left; other keys are the caller's.
    assert!(matches!(b.key(key(KeyCode::Backspace)), Nav::Moved));
    assert_eq!(b.cwd, root);
    assert_eq!(b.highlighted().unwrap().name, "ann");
    for code in [KeyCode::Esc, KeyCode::Char(' '), KeyCode::Tab] {
        assert!(matches!(b.key(key(code)), Nav::Ignored), "{code:?}");
    }
}

#[test]
fn typing_narrows_the_listing_and_esc_restores_it() {
    let tmp = tree();
    let mut b = gff_browser(tmp.path());
    // Any case, the cursor on the first name starting with the text, and
    // `..` hidden.
    typed(&mut b, "AN");
    assert_eq!(names(&b), ["ann"]);
    assert_eq!(b.highlighted().map(|e| e.dir), Some(true));
    assert!(matches!(b.key(key(KeyCode::Backspace)), Nav::Moved));
    assert_eq!(names(&b), ["ann"], "`a` still matches ann alone");
    typed(&mut b, "zz");
    assert!(b.list.is_empty());
    assert!(matches!(b.key(key(KeyCode::Esc)), Nav::Moved));
    assert_eq!(names(&b), ["..", "ann", "b.gff3"]);
    // With nothing typed, Esc is the caller's.
    assert!(matches!(b.key(key(KeyCode::Esc)), Nav::Ignored));
}

#[test]
fn typing_a_path_walks_there() {
    let tmp = tree();
    let root = normalize(tmp.path());
    let mut b = gff_browser(&std::env::temp_dir());
    // `/` first goes to the root; each name and `/` goes one folder down.
    typed(&mut b, &format!("{}/", root.display()));
    assert_eq!(b.cwd, root);
    assert!(b.list.find.text.is_empty());
    // A partial name goes to the highlighted folder; `../` back up.
    typed(&mut b, "an/");
    assert_eq!(b.cwd, root.join("ann"));
    typed(&mut b, "../");
    assert_eq!(b.cwd, root);
    // `~` is home.
    if let Some(home) = std::env::var_os("HOME") {
        typed(&mut b, "~");
        assert_eq!(b.cwd, normalize(Path::new(&home)));
    }
}

#[test]
fn a_relative_start_is_made_absolute_so_up_keeps_working() {
    let mut b = Browser::new(PathBuf::from("."), |_| false, |_| false).opened();
    assert!(b.cwd.is_absolute(), "{}", b.cwd.display());
    let start = b.cwd.clone();
    assert!(matches!(b.key(key(KeyCode::Left)), Nav::Moved));
    assert_eq!(Some(b.cwd.as_path()), start.parent());
    assert!(!b.list.is_empty());
}
