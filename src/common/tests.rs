use super::*;

#[test]
fn a_backend_says_which_batch_it_holds_and_who_wrote_it() -> anyhow::Result<()> {
    let dir = tempfile::tempdir()?;
    let path = dir.path().join("s1_count.zarr");
    let path = path.to_str().expect("utf8");
    let triplets = TripletsRowsCols {
        triplets: vec![(0, 0, 1.0), (1, 1, 2.0)],
        rows: vec!["GENE1/count/spliced".into(), "GENE1/count/unspliced".into()],
        cols: vec!["c1".into(), "c2".into()],
    };
    let mut metadata = backend_meta("count", "s1");
    metadata.insert(meta::CONTENT.to_string(), meta::GENE_COUNT.to_string());
    drop(triplets.to_backend(path, &metadata)?);

    let data = open_sparse_matrix(path, &SparseIoBackend::Zarr)?;
    assert_eq!(data.meta(meta::SAMPLE).as_deref(), Some("s1"));
    assert_eq!(data.meta(meta::PRODUCER).as_deref(), Some("faba count"));
    assert_eq!(data.meta(meta::CONTENT).as_deref(), Some(meta::GENE_COUNT));
    Ok(())
}

#[test]
fn a_whole_write_takes_the_files_place() -> anyhow::Result<()> {
    use std::io::Write;
    let dir = tempfile::tempdir()?;
    let path = dir.path().join("a b.txt");
    std::fs::write(&path, "old")?;
    write_replacing(&path, |w| Ok(w.write_all(b"new")?))?;
    assert_eq!(std::fs::read_to_string(&path)?, "new");
    assert!(
        !part_of(&path).exists(),
        "the part file is renamed into place"
    );
    assert_eq!(part_of(&path), dir.path().join("a b.txt.part"));
    Ok(())
}

#[test]
fn a_failed_write_keeps_the_old_file_and_leaves_no_part() -> anyhow::Result<()> {
    use std::io::Write;
    let dir = tempfile::tempdir()?;
    let path = dir.path().join("a.txt");
    std::fs::write(&path, "old")?;
    let failed = write_replacing(&path, |w| {
        w.write_all(b"half")?;
        w.flush()?;
        anyhow::bail!("stopped part way")
    });
    assert!(failed.is_err());
    assert_eq!(std::fs::read_to_string(&path)?, "old");
    assert_eq!(
        std::fs::read_dir(dir.path())?.count(),
        1,
        "no part file left"
    );
    Ok(())
}
