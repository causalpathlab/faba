use super::*;

fn picture(w: u32, h: u32) -> RgbaImage {
    RgbaImage::from_fn(w, h, |x, y| {
        image::Rgba([(x % 256) as u8, (y % 256) as u8, 90, 255])
    })
}

#[test]
fn saves_are_listed_newest_first_with_thumbnails_and_last() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().join(".faba-view");
    let file = |name: &str| {
        let p = tmp.path().join(name);
        std::fs::write(&p, b"%PDF").unwrap();
        p
    };
    let (a, b) = (file("a.pdf"), file("b.pdf"));
    let mut g = Gallery::open(&dir);
    g.add(&a, "pileup", &picture(720, 360)).unwrap();
    g.add(&b, "metagene", &picture(680, 300)).unwrap();
    let names: Vec<String> = g.entries().iter().map(Entry::name).collect();
    assert_eq!(names, ["b.pdf", "a.pdf"]);
    let thumb = image::open(g.thumb(&g.entries()[1])).unwrap();
    assert_eq!(
        (thumb.width(), thumb.height()),
        (240, 120),
        "at the picture's aspect"
    );

    // Saving a file again moves it to the top, with one thumbnail.
    g.add(&a, "pileup", &picture(720, 360)).unwrap();
    let reopened = Gallery::open(&dir);
    let names: Vec<String> = reopened.entries().iter().map(Entry::name).collect();
    assert_eq!(names, ["a.pdf", "b.pdf"], "the log lasts");
    assert_eq!(std::fs::read_dir(dir.join("thumbs")).unwrap().count(), 2);

    // A deleted file drops out.
    std::fs::remove_file(&b).unwrap();
    assert_eq!(Gallery::open(&dir).entries().len(), 1);
}

#[test]
fn ago_reads_as_a_person_would() {
    assert_eq!(ago(100, 130), "just now");
    assert_eq!(ago(0, 600), "10 min ago");
    assert_eq!(ago(0, 7200), "2 h ago");
    assert_eq!(ago(0, 3 * 86_400), "3 d ago");
}
