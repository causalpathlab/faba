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
