use super::*;

fn names(v: &[&str]) -> Vec<String> {
    v.iter().map(ToString::to_string).collect()
}

#[test]
fn listings_give_their_names() {
    let html = r#"<a href="?C=N;O=D">Name</a> <a href="/pub/">Parent</a>
        <a href="release_49/">release_49/</a> <a href="_README.TXT">x</a>
        <a href="https://elsewhere/">y</a>"#;
    assert_eq!(hrefs(html), ["release_49", "_README.TXT"]);
}

#[test]
fn gencode_releases_come_newest_first_and_ensembl_species_by_name() {
    let g = Source::Gencode("Gencode_mouse");
    let got = g.items(names(&[
        "latest_release",
        "release_M9",
        "release_M39",
        "_README.TXT",
    ]));
    assert_eq!(got, ["release_M39", "release_M9"]);
    let got = Source::Ensembl.items(names(&["mus_musculus", "danio_rerio", "README"]));
    assert_eq!(got, ["danio_rerio", "mus_musculus"]);
}

#[test]
fn the_primary_assembly_files_are_chosen() {
    let g = names(&[
        "GRCx.p1.genome.fa.gz",
        "GRCx.primary_assembly.genome.fa.gz",
        "gencode.v1.primary_assembly.basic.annotation.gtf.gz",
        "gencode.v1.primary_assembly.annotation.gff3.gz",
        "gencode.v1.primary_assembly.annotation.gtf.gz",
    ]);
    assert_eq!(
        gencode_files(&g),
        Some((
            "gencode.v1.primary_assembly.annotation.gtf.gz".into(),
            "GRCx.primary_assembly.genome.fa.gz".into()
        ))
    );
    let gtf = names(&[
        "Sp.A1.9.abinitio.gtf.gz",
        "Sp.A1.9.chr.gtf.gz",
        "Sp.A1.9.chr_patch_hapl_scaff.gtf.gz",
        "Sp.A1.9.gtf.gz",
    ]);
    assert_eq!(ensembl_gtf(&gtf).as_deref(), Some("Sp.A1.9.gtf.gz"));
    let dna = names(&[
        "Sp.A1.dna.toplevel.fa.gz",
        "Sp.A1.dna_sm.primary_assembly.fa.gz",
    ]);
    assert_eq!(
        ensembl_genome(&dna).as_deref(),
        Some("Sp.A1.dna.toplevel.fa.gz")
    );
    let dna = names(&[
        "Sp.A1.dna.toplevel.fa.gz",
        "Sp.A1.dna.primary_assembly.fa.gz",
    ]);
    assert_eq!(
        ensembl_genome(&dna).as_deref(),
        Some("Sp.A1.dna.primary_assembly.fa.gz")
    );
}

#[test]
fn files_already_there_are_kept_and_the_genome_is_unpacked_and_indexed() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().join("ref");
    std::fs::create_dir(&dir).unwrap();
    let (gtf_url, genome_url) = ("https://x/a.gtf.gz", "https://x/g.fa.gz");
    std::fs::write(dir.join("a.gtf.gz"), b"").unwrap();
    // The genome as fetched, gzipped: unpacked, indexed, and the .gz gone.
    let gz = dir.join("g.fa.gz");
    let mut enc = flate2::write::GzEncoder::new(File::create(&gz).unwrap(), Default::default());
    enc.write_all(b">chr1\nACGTACGT\n>chr2\nTTTT\n").unwrap();
    enc.finish().unwrap();
    let said = |_: &str, _| {};
    let (gtf, genome) = fetch(gtf_url, genome_url, &dir, &Stopper::default(), &said).unwrap();
    assert_eq!(
        (gtf, genome.clone()),
        (dir.join("a.gtf.gz"), dir.join("g.fa"))
    );
    assert_eq!(
        std::fs::read_to_string(&genome).unwrap(),
        ">chr1\nACGTACGT\n>chr2\nTTTT\n"
    );
    assert!(!gz.exists());
    assert!(PathBuf::from(format!("{}.fai", genome.display())).exists());
}

#[test]
fn a_plain_genome_already_there_is_kept_as_it_is() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().to_path_buf();
    std::fs::write(dir.join("a.gtf.gz"), b"").unwrap();
    std::fs::write(dir.join("g.fa"), b">chr1\nACGT\n").unwrap();
    let said = |_: &str, _| {};
    let (_, genome) = fetch(
        "https://x/a.gtf.gz",
        "https://x/g.fa",
        &dir,
        &Stopper::default(),
        &said,
    )
    .unwrap();
    assert_eq!(std::fs::read_to_string(&genome).unwrap(), ">chr1\nACGT\n");
}

#[test]
fn the_catalogue_walks_sources_then_items_and_finds_by_typing() {
    use ratatui::crossterm::event::{KeyCode, KeyModifiers};
    let k = |c| KeyEvent::new(KeyCode::Char(c), KeyModifiers::NONE);
    // BAMs naming chromosomes `1` start on the Ensembl source.
    let mut c = Catalogue::new(Some(false));
    assert_eq!(SOURCES[c.at].1, Source::Ensembl);
    assert!(!c.key(&k('x')), "no find among the sources");
    assert!(c.key(&KeyEvent::new(KeyCode::Up, KeyModifiers::NONE)));
    assert_eq!(c.at, 1);
    // Human GENCODE chosen, with its listing as it would arrive.
    c.at = 0;
    c.source = Some(SOURCES[0].1);
    let listing = std::thread::spawn(|| Ok(names(&["release_50", "release_49"])));
    c.pending = Some((listing, Arc::default()));
    while !c.poll() {
        std::thread::sleep(std::time::Duration::from_millis(5));
    }
    assert!(!c.loading());
    for ch in "49".chars() {
        assert!(c.key(&k(ch)));
    }
    assert_eq!(c.items.shown().collect::<Vec<_>>(), ["release_49"]);
    let pick = c.enter().unwrap();
    assert_eq!(pick.dir_name(), "gencode_human_release_49");
    assert!(c.back() && c.source.is_none() && c.at == 0);
    assert!(!c.back());
}

#[test]
#[ignore = "reads the GENCODE and Ensembl FTP listings"]
fn the_ftp_listings_resolve_to_files() {
    for (_, s) in SOURCES {
        let stopper = Stopper::default();
        let items = s.items(list(&s.catalogue_url(), &stopper).unwrap());
        let item = match s {
            Source::Gencode(_) => items[0].clone(),
            Source::Ensembl => "danio_rerio".into(),
        };
        let (gtf, genome) = Pick { source: s, item }.resolve(&stopper).unwrap();
        eprintln!("{gtf}\n{genome}");
        assert!(gtf.ends_with(".gtf.gz") && genome.ends_with(".fa.gz"));
    }
}
