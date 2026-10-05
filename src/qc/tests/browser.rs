use super::*;
use ratatui::crossterm::event::KeyModifiers;

fn key(code: KeyCode) -> KeyEvent {
    KeyEvent::new(code, KeyModifiers::NONE)
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

#[test]
fn lists_directories_then_the_kept_files_and_picks_a_file() {
    let tmp = tree();
    let root = tmp.path().to_path_buf();
    let mut tag = |p: &Path| p.extension().is_some_and(|e| e == "gff3");
    let mut listing = Listing {
        keep: &gff,
        tag: &mut tag,
    };
    let mut b = Browser::new(root.clone());
    b.open(root.clone(), &mut listing);
    let names: Vec<&str> = b.entries.iter().map(|e| e.name.as_str()).collect();
    // Hidden files and zarr stores are left out; other files unless kept.
    assert_eq!(names, ["..", "ann", "b.gff3"]);
    assert!(b.entries[2].tagged && !b.entries[1].tagged);

    // Into `ann`, where only the annotation is listed, and pick it.
    b.at = 1;
    assert!(matches!(
        b.key(key(KeyCode::Enter), &mut listing),
        Nav::Moved
    ));
    assert_eq!(b.cwd, root.join("ann"));
    let names: Vec<&str> = b.entries.iter().map(|e| e.name.as_str()).collect();
    assert_eq!(names, ["..", "genes.gtf.gz"]);
    b.at = 1;
    // Right only opens directories; Enter picks a file.
    assert!(matches!(
        b.key(key(KeyCode::Right), &mut listing),
        Nav::Moved
    ));
    assert!(
        matches!(b.key(key(KeyCode::Enter), &mut listing), Nav::Picked(p) if p == root.join("ann/genes.gtf.gz"))
    );

    // Up lands back on the directory just left; other keys are the caller's.
    assert!(matches!(
        b.key(key(KeyCode::Backspace), &mut listing),
        Nav::Moved
    ));
    assert_eq!(b.cwd, root);
    assert_eq!(b.entries[b.at].name, "ann");
    assert!(matches!(
        b.key(key(KeyCode::Esc), &mut listing),
        Nav::Ignored
    ));
}

#[test]
fn a_relative_start_is_made_absolute_so_up_keeps_working() {
    let mut no_tag = |_: &Path| false;
    let mut listing = Listing {
        keep: &|_| false,
        tag: &mut no_tag,
    };
    let mut b = Browser::new(PathBuf::from("."));
    b.open(PathBuf::from("."), &mut listing);
    assert!(b.cwd.is_absolute(), "{}", b.cwd.display());
    let start = b.cwd.clone();
    assert!(matches!(
        b.key(key(KeyCode::Left), &mut listing),
        Nav::Moved
    ));
    assert_eq!(Some(b.cwd.as_path()), start.parent());
    assert!(!b.entries.is_empty());
}
