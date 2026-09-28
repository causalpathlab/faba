use super::*;

fn read(path: &Path) -> Value {
    serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap()
}

#[test]
fn records_inputs_new_outputs_and_options() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().to_str().unwrap();
    std::fs::write(dir.path().join("old.txt"), "x").unwrap();
    let gff = dir.path().join("genes.gff");
    std::fs::write(&gff, "").unwrap();

    let record = RunRecord::start("job1", out)
        .input("gff", gff.to_str())
        .input("missing", None)
        .inputs("bam", &["a.bam", "b.bam"]);
    std::fs::write(dir.path().join("new.parquet"), "y").unwrap();
    let options = json!({"min_coverage": 3});
    let path = record.options(&options).finish(&Ok(())).unwrap();

    assert_eq!(path, dir.path().join("job1.run.json"));
    let v = read(&path);
    assert_eq!(v["job"], "job1");
    assert_eq!(v["status"], "ok");
    assert_eq!(v["options"]["min_coverage"], 3);
    assert_eq!(v["outputs"], json!(["new.parquet"]));
    assert!(v["inputs"].get("missing").is_none());
    assert_eq!(v["inputs"]["bam"].as_array().unwrap().len(), 2);
    assert!(Path::new(v["inputs"]["gff"]["path"].as_str().unwrap()).is_absolute());
}

#[cfg(unix)]
#[test]
fn a_symlinked_input_is_recorded_as_given() {
    let dir = tempfile::tempdir().unwrap();
    let real = dir.path().join("genes.v2.gff");
    std::fs::write(&real, "").unwrap();
    let link = dir.path().join("genes.gff");
    std::os::unix::fs::symlink(&real, &link).unwrap();
    let path = RunRecord::start("job4", dir.path().to_str().unwrap())
        .input("gff", link.to_str())
        .finish(&Ok(()))
        .unwrap();
    let gff = &read(&path)["inputs"]["gff"];
    assert_eq!(gff["path"], link.to_str().unwrap());
    assert_eq!(
        gff["resolved"],
        real.canonicalize().unwrap().to_str().unwrap()
    );
}

#[cfg(unix)]
#[test]
fn find_input_follows_the_file_the_run_read() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().to_str().unwrap();
    let (v1, v2) = (
        dir.path().join("genes.v1.gff"),
        dir.path().join("genes.v2.gff"),
    );
    std::fs::write(&v1, "").unwrap();
    std::fs::write(&v2, "").unwrap();
    let link = dir.path().join("genes.gff");
    std::os::unix::fs::symlink(&v1, &link).unwrap();
    RunRecord::start("job5", out)
        .input("gff", link.to_str())
        .finish(&Ok(()))
        .unwrap();
    let found = || find_recorded(out, "gff").unwrap().0;

    // The link unchanged: its own name.
    assert_eq!(found(), link.to_str().unwrap());
    // The link repointed: the file the run read.
    std::fs::remove_file(&link).unwrap();
    std::os::unix::fs::symlink(&v2, &link).unwrap();
    assert_eq!(Path::new(&found()), v1.canonicalize().unwrap());
    // The link gone: still the file the run read.
    std::fs::remove_file(&link).unwrap();
    assert_eq!(Path::new(&found()), v1.canonicalize().unwrap());
    // Both gone: nothing.
    std::fs::remove_file(&v1).unwrap();
    assert!(find_recorded(out, "gff").is_none());
}

#[test]
fn a_plain_string_input_is_still_read() {
    let dir = tempfile::tempdir().unwrap();
    let gff = dir.path().join("genes.gff");
    std::fs::write(&gff, "").unwrap();
    let record = json!({"inputs": {"gff": gff.to_str().unwrap()}});
    std::fs::write(dir.path().join("old.run.json"), record.to_string()).unwrap();
    let (found, _) = find_recorded(dir.path().to_str().unwrap(), "gff").unwrap();
    assert_eq!(found, gff.to_str().unwrap());
}

#[test]
fn a_failed_run_is_recorded_and_its_error_kept() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().to_str().unwrap();
    let outcome: anyhow::Result<()> =
        recorded(RunRecord::start("job2", out), || anyhow::bail!("no reads"));
    assert!(outcome.is_err());
    let v = read(&dir.path().join("job2.run.json"));
    assert_eq!(v["status"], "failed: no reads");
}

#[test]
fn find_input_reads_the_record_next_to_an_output() {
    let dir = tempfile::tempdir().unwrap();
    let out = dir.path().to_str().unwrap();
    let gff = dir.path().join("genes.gff");
    std::fs::write(&gff, "").unwrap();
    RunRecord::start("job3", out)
        .input("gff", gff.to_str())
        .input("genome", Some("/nowhere/genome.fa"))
        .finish(&Ok(()))
        .unwrap();
    let matrix = dir.path().join("b1_m6a.zarr");
    std::fs::create_dir(&matrix).unwrap();

    for near in [out.to_string(), matrix.to_string_lossy().into_owned()] {
        let (found, record) = find_recorded(&near, "gff").expect("gff is recorded");
        assert_eq!(Path::new(&found), gff);
        assert_eq!(record, dir.path().join("job3.run.json"));
    }
    // A recorded file that no longer exists is not offered.
    assert!(find_recorded(out, "genome").is_none());
    assert!(find_recorded(out, "known_snps").is_none());
}

#[test]
fn run_records_are_recognized_by_name() {
    assert!(is_run_record("count.run.json"));
    assert!(is_run_record("pipeline_summary.json"));
    assert!(!is_run_record("qc_summary.tsv"));
    assert!(!is_run_record("run.json.bak"));
}

#[test]
fn a_file_name_and_an_updated_option_are_written() {
    let dir = tempfile::tempdir().unwrap();
    let mut record = RunRecord::start("job6", dir.path().to_str().unwrap())
        .file_name("summary.json")
        .options(&json!({"a": 1, "b": 2}));
    record.set_option("b", &3);
    let path = record.finish(&Ok(())).unwrap();
    assert_eq!(path, dir.path().join("summary.json"));
    let v = read(&path);
    assert_eq!(v["options"], json!({"a": 1, "b": 3}));
}
