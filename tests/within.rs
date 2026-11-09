//! The two rounds of `faba cluster` on a synthetic output directory, through
//! the built binary: cell modules by gene counts, then m6A within each, and
//! `faba cluster-view` of what the second round records. What the chains
//! find is not checked, only that every step reads what the one before wrote.

use std::path::Path;
use std::process::Command;

use data_beans::sparse_io::{create_sparse_from_triplets_owned, SparseIoBackend};

const CELLS: usize = 60;
const GENES: usize = 30;

/// `faba` with `args`, which must succeed: what it printed.
fn faba(args: &[&str]) -> String {
    let out = Command::new(env!("CARGO_BIN_EXE_faba"))
        .args(args)
        .stdin(std::process::Stdio::null())
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "faba {args:?}\n{}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8_lossy(&out.stdout).into_owned()
}

/// Write `{dir}/{name}.zarr` from named triplets.
fn write_matrix(
    dir: &Path,
    name: &str,
    rows: &[String],
    cols: &[String],
    entries: &[(usize, usize, f32)],
) {
    let path = dir.join(format!("{name}.zarr"));
    let triplets: Vec<(u64, u64, f32)> = entries
        .iter()
        .map(|&(r, c, x)| (r as u64, c as u64, x))
        .collect();
    let mut data = create_sparse_from_triplets_owned(
        triplets,
        (rows.len(), cols.len(), entries.len()),
        Some(&path.to_string_lossy()),
        Some(&SparseIoBackend::Zarr),
    )
    .unwrap();
    let boxed =
        |names: &[String]| -> Vec<Box<str>> { names.iter().map(|s| s.as_str().into()).collect() };
    data.register_row_names_vec(&boxed(rows));
    data.register_column_names_vec(&boxed(cols));
}

/// One batch of 60 cells of two kinds, half each: a kind's half of the 30
/// genes has 8 to 10 UMIs a cell, the other half 1 to 3. Eight m6A sites on
/// four genes, 6 reads a cell each: in every other cell the first two
/// genes' sites are mostly methylated, in the rest the last two's.
fn fixture(dir: &Path) {
    let cells: Vec<String> = (0..CELLS).map(|c| format!("CELL{c:02}-1")).collect();
    let genes: Vec<String> = (0..GENES).map(|g| format!("G{g}/count/spliced")).collect();
    let mut counts = Vec::new();
    for c in 0..CELLS {
        for g in 0..GENES {
            let high = (c < CELLS / 2) == (g < GENES / 2);
            let umis = if high { 8 } else { 1 } + (c * 7 + g * 3) % 3;
            counts.push((g, c, umis as f32));
        }
    }
    write_matrix(dir, "b1_count", &genes, &cells, &counts);

    let mut sites = Vec::new();
    for s in 0..8 {
        let at = format!("G{}/m6a/chr1:{}", s / 2, 100 * (s + 1));
        sites.push(format!("{at}/methylated"));
        sites.push(format!("{at}/unmethylated"));
    }
    let mut reads = Vec::new();
    for c in 0..CELLS {
        for s in 0..8 {
            let methylated = if (c % 2 == 0) == (s < 4) { 5 } else { 1 };
            reads.push((2 * s, c, methylated as f32));
            reads.push((2 * s + 1, c, (6 - methylated) as f32));
        }
    }
    write_matrix(dir, "b1_m6a_site", &sites, &cells, &reads);
}

#[test]
fn the_second_round_fits_within_the_firsts_cell_modules() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    fixture(dir);
    let at = |name: &str| dir.join(name).to_string_lossy().into_owned();
    // Short chains over everything written: the filters are for real depth.
    let chains = [
        "--batch-process",
        "--chains",
        "2",
        "--sweeps",
        "40",
        "--burn-in",
        "20",
    ]
    .into_iter()
    .chain(["--min-cell-sites", "1", "--min-site-cells", "1"]);

    let mut count = vec!["cluster", "--modality", "count"];
    let folder = at("");
    count.push(&folder);
    count.extend(chains.clone());
    faba(&count);
    let count_record = at("count_blockmodel.run.json");
    assert!(Path::new(&count_record).is_file(), "no {count_record}");

    // The sites are those beside the count fit's matrices: none are named.
    let mut within = vec![
        "cluster",
        "--within",
        &count_record,
        "--within-min-cells",
        "5",
    ];
    within.extend(chains);
    faba(&within);
    let record = at("m6a_within.run.json");
    let json: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(&record).unwrap()).unwrap();
    let nested = &json["results"]["within"];
    assert_eq!(nested["modality"], "m6a", "{nested:#}");
    assert_eq!(nested["cells"]["read"], CELLS, "{nested:#}");
    assert!(dir.join(nested["table"].as_str().unwrap()).is_file());
    let groups = nested["groups"].as_array().unwrap();
    let cells: u64 = groups.iter().map(|g| g["cells"].as_u64().unwrap()).sum();
    assert_eq!(
        cells,
        nested["cells"]["matched"].as_u64().unwrap(),
        "{nested:#}"
    );
    let fitted: Vec<&serde_json::Value> =
        groups.iter().filter(|g| g["status"] == "fitted").collect();
    assert!(!fitted.is_empty(), "{nested:#}");
    for group in groups {
        assert!(
            group["status"] == "fitted" || group["status"] == "skipped",
            "{group:#}"
        );
    }
    for fit in &fitted {
        let sub = dir.join(fit["dir"].as_str().unwrap());
        assert!(
            sub.join(fit["record"].as_str().unwrap()).is_file(),
            "{fit:#}"
        );
        assert!(sub.join("m6a_block_sites.parquet").is_file(), "{fit:#}");
    }

    // The view of the record: every fit printed, one fit by its group, and
    // a page a fit.
    let printed = faba(&["cluster-view", &record, "--batch-process"]);
    assert!(printed.starts_with("count fit: "), "{printed}");
    assert!(printed.contains("group\tcells\tfit"), "{printed}");
    let name = fitted[0]["group"].as_str().unwrap();
    assert!(
        printed.contains(&format!("m6a within {name}\nmodule\tmain genes")),
        "{printed}"
    );
    let one = faba(&["cluster-view", &record, "--group", name, "--batch-process"]);
    assert!(
        one.starts_with("module\tmain genes") && one.contains(" %"),
        "{one}"
    );
    let page = at("fit.html");
    faba(&["cluster-view", &record, "--html", &page]);
    assert!(Path::new(&page).is_file());
    for fit in &fitted {
        let group = fit["group"].as_str().unwrap();
        assert!(dir.join(format!("fit.{group}.html")).is_file(), "{group}");
    }
}
