use super::*;

fn argv(bams: &[&str]) -> Vec<String> {
    let mut v = vec!["run".to_string(), "--batch-process".to_string()];
    v.extend(bams.iter().map(|s| s.to_string()));
    v.extend(
        [
            "-g",
            "/r/genes.gff",
            "-f",
            "/r/genome.fa",
            "-o",
            ".",
            "--max-threads",
            "4",
        ]
        .map(String::from),
    );
    v
}

#[test]
fn the_script_text() {
    let t = text(&argv(&["/d/sample_A.bam"]));
    assert!(t.starts_with("#!/usr/bin/env bash\n# Made by `faba run` (faba "));
    assert!(t.contains("set -euo pipefail\ncd \"$(dirname \"$0\")\"\n"));
    assert!(t.contains("if [ -f pipeline_summary.json ]; then"));
    assert!(t.contains(
        "\"${FABA:-faba}\" run --batch-process \\\n  /d/sample_A.bam \\\n  -g /r/genes.gff \\\n"
    ));
    assert!(t.trim_end().ends_with("--max-threads 4"));
}

#[test]
fn odd_paths_are_quoted() {
    assert_eq!(quote("/a b/it's.bam"), r"'/a b/it'\''s.bam'");
    assert_eq!(quote("$HOME"), "'$HOME'");
    assert_eq!(quote("/plain/x.bam"), "/plain/x.bam");
}

#[test]
fn never_overwrites_and_runs_as_bash() {
    let tmp = tempfile::tempdir().unwrap();
    let odd = tmp.path().join("a b'$x");
    std::fs::create_dir_all(&odd).unwrap();
    // A fake faba that prints its arguments, one per line.
    let fake = tmp.path().join("fake-faba");
    std::fs::write(&fake, "#!/usr/bin/env bash\nprintf '%s\\n' \"$@\"\n").unwrap();
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        std::fs::set_permissions(&fake, std::fs::Permissions::from_mode(0o755)).unwrap();
    }
    let bam = odd.join("sample_A.bam").to_string_lossy().into_owned();
    let path = write(tmp.path(), &argv(&[&bam])).unwrap();
    assert!(
        write(tmp.path(), &argv(&[&bam])).is_err(),
        "never over a script"
    );
    let out = std::process::Command::new("bash")
        .arg(&path)
        .env("FABA", &fake)
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let said: Vec<String> = String::from_utf8_lossy(&out.stdout)
        .lines()
        .map(String::from)
        .collect();
    assert_eq!(said[2], bam);

    std::fs::write(tmp.path().join(GUARD), "{}").unwrap();
    let out = std::process::Command::new("bash")
        .arg(&path)
        .env("FABA", &fake)
        .output()
        .unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("pipeline_summary.json exists"));
}

#[test]
fn several_bams_one_per_line() {
    let lines = command_lines(&vec![
        "run".to_string(),
        "--batch-process".to_string(),
        "/d/sample_A.bam".to_string(),
        "/d/sample_B.bam".to_string(),
        "/d/sample_C.bam".to_string(),
        "-g".to_string(),
        "/r/genes.gff".to_string(),
    ]);
    // First line: program and flags
    assert_eq!(lines[0], "\"${FABA:-faba}\" run --batch-process");
    // Each BAM on its own line
    assert_eq!(lines[1], "/d/sample_A.bam");
    assert_eq!(lines[2], "/d/sample_B.bam");
    assert_eq!(lines[3], "/d/sample_C.bam");
    // Flag with value
    assert_eq!(lines[4], "-g /r/genes.gff");
}

#[test]
fn flag_values_on_same_line() {
    let lines = command_lines(&vec![
        "run".to_string(),
        "--batch-process".to_string(),
        "/d/sample.bam".to_string(),
        "-x".to_string(),
        "a".to_string(),
        "b".to_string(),
        "c".to_string(),
        "-y".to_string(),
        "d".to_string(),
    ]);
    // All values for -x should be on the same line
    assert_eq!(lines[2], "-x a b c");
    // All values for -y should be on the same line
    assert_eq!(lines[3], "-y d");
}

#[test]
fn negative_number_value_stays_with_flag() {
    let lines = command_lines(&vec![
        "run".to_string(),
        "--batch-process".to_string(),
        "/d/sample.bam".to_string(),
        "--min".to_string(),
        "-1".to_string(),
        "--max".to_string(),
        "-0.5".to_string(),
    ]);
    // Negative numbers should stay on the same line as their flags
    assert_eq!(lines[2], "--min -1");
    assert_eq!(lines[3], "--max -0.5");
}

#[test]
fn flag_and_next_flag_separate_lines() {
    let lines = command_lines(&vec![
        "run".to_string(),
        "--batch-process".to_string(),
        "/d/sample.bam".to_string(),
        "-o".to_string(),
        ".".to_string(),
        "--max-threads".to_string(),
        "4".to_string(),
    ]);
    // -o with value
    assert!(lines.iter().any(|l| l == "-o ."));
    // --max-threads with value on separate line
    assert!(lines.iter().any(|l| l == "--max-threads 4"));
}

#[test]
fn empty_string_quotes() {
    assert_eq!(quote(""), "''");
}
