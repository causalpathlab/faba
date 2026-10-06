use super::*;

fn touch(dir: &Path, name: &str) {
    std::fs::write(dir.join(name), b"").unwrap();
}

#[test]
fn a_directory_gives_its_site_matrices_by_batch_and_its_site_table() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    for f in [
        "b2_m6a_site.zarr.zip",
        "b1_m6a_site.zarr.zip",
        "b1_m6a.zarr.zip",
        "b1_atoi_site.zarr.zip",
        "m6a_sites.parquet",
        "pipeline_summary.json",
    ] {
        touch(dir, f);
    }
    let m6a = output_sites(dir, None).unwrap();
    assert_eq!(&*m6a.modality, "m6a", "m6A first");
    let batches: Vec<&str> = m6a.matrices.iter().map(|(b, _)| &**b).collect();
    assert_eq!(batches, ["b1", "b2"], "site matrices only, by batch");
    assert!(m6a.site_table.unwrap().ends_with("m6a_sites.parquet"));
    let atoi = output_sites(dir, Some("atoi")).unwrap();
    assert_eq!(atoi.matrices.len(), 1);
    assert!(atoi.site_table.is_none());
    assert!(output_sites(dir, Some("apa")).is_err());

    // A run record names its folder; a matrix or another file does not.
    let record = dir.join("pipeline_summary.json");
    let named = crate::run_record::output_dir_of(&record.to_string_lossy()).unwrap();
    assert_eq!(named, dir);
    let matrix = dir.join("b1_m6a.zarr.zip");
    assert!(crate::run_record::output_dir_of(&matrix.to_string_lossy()).is_none());
}
