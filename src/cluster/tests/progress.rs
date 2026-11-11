//! A site fit's progress: its trace saved as the progress noted it and its
//! record saying the schedule, a cancel and a commit excluding each other,
//! a cancelled fit writing nothing, and the MCMC fits' flags refused before
//! anything is read.

use super::planted::planted_site_data;
use crate::cluster::count::vi::KMEANS_ITERS;
use crate::cluster::progress::{advice, Cancelled, Progress};
use crate::cluster::read::Weighted;
use crate::cluster::FitData;

/// Once the fit commits to writing, a cancel is refused; a cancel first
/// makes the commit fail.
#[test]
fn a_cancel_and_a_commit_exclude_each_other() {
    let verdict = || advice::Verdict::vi(true, false, (3, 10), 1, (0.0, "read"));
    let p = Progress::new(1);
    p.optimised(verdict()).unwrap();
    assert!(!p.cancel() && !p.cancelled());
    let p = Progress::new(1);
    assert!(p.cancel() && p.cancelled());
    assert!(p.optimised(verdict()).unwrap_err().is::<Cancelled>());
}

/// The site blockmodel's fit of planted reads (25 cells fitted, 5 placed
/// after, 6 genes) by two restarts under `progress`, into `dir`, filling
/// `record`: what it returned.
fn fit_recorded(
    progress: &Progress,
    dir: &std::path::Path,
    record: &mut crate::run_record::RunRecord,
) -> anyhow::Result<crate::cluster::Written> {
    let (data, fit_reads) = planted_site_data(30, 25, 6, 171);
    let weighted = Weighted::new(&fit_reads, None);
    let fitted: Vec<usize> = (0..25).collect();
    let kept: Vec<usize> = (0..30).collect();
    let controls: Vec<usize> = (25..30).collect();
    let args = super::parse_args(&[
        "cluster",
        "x.zarr",
        "--chains",
        "2",
        "--cell-roles",
        "3",
        "--gene-roles",
        "2",
        "--epochs",
        "8",
    ]);
    let fit = FitData {
        data: &data,
        fitted: &fitted,
        weighted: &weighted,
        all_cells: &data.cells,
        kept: &kept,
        controls: &controls,
        progress,
    };
    crate::cluster::run::sites::blocks_fit(&args, dir, record, fit)
}

/// A site fit saves each restart's ELBO beside its tables, a row per point
/// as the progress noted it, and its record says the schedule the points
/// ran in (the k-means start, then the epochs) and the verdict the views
/// show.
#[test]
fn a_site_fit_saves_its_trace() {
    use crate::common::column;
    use arrow::array::{Float64Array, Int32Array, Int64Array};
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().join("out");
    std::fs::create_dir_all(&dir).unwrap();
    let p = Progress::new(2);
    let mut record = crate::run_record::RunRecord::start("cluster", &dir.to_string_lossy());
    let written = fit_recorded(&p, &dir, &mut record).unwrap();
    record = record.file_name(&written.record);
    record.finish(&Ok(()));
    let path = dir.join(format!("{}_trace.parquet", written.modality));
    let t = crate::common::read_record_batch(&path).unwrap();
    let chain: &Int32Array = column(&t, "chain").unwrap();
    let step: &Int64Array = column(&t, "sweep").unwrap();
    let elbo: &Float64Array = column(&t, "log_joint").unwrap();
    let rows: Vec<(usize, usize, u64)> = (0..t.num_rows())
        .map(|i| {
            (
                chain.value(i) as usize,
                step.value(i) as usize,
                elbo.value(i).to_bits(),
            )
        })
        .collect();
    let want: Vec<(usize, usize, u64)> = (p.traces().iter().enumerate())
        .flat_map(|(i, t)| t.iter().map(move |&(n, v)| (i, n, v.to_bits())))
        .collect();
    assert_eq!((0..2).filter(|&i| !p.trace(i).is_empty()).count(), 2);
    assert!(want.iter().all(|r| f64::from_bits(r.2).is_finite()));
    assert_eq!(rows, want);
    let text = std::fs::read_to_string(dir.join(&written.record)).unwrap();
    let json: serde_json::Value = serde_json::from_str(&text).unwrap();
    assert_eq!(
        json["results"]["trace"],
        serde_json::json!({"sweeps": KMEANS_ITERS + 8, "burn_in": KMEANS_ITERS, "quantity": "ELBO"})
    );
    let fit = &json["results"]["fit"];
    assert_eq!(
        (&fit["method"], &fit["unit"]),
        (&"vi".into(), &"read".into())
    );
    assert_eq!(fit["restarts"].as_array().map(Vec::len), Some(2));
    assert_eq!(
        json["results"]["converged"].as_bool(),
        p.verdict().map(|v| v.converged()),
        "the verdict the views show"
    );
}

/// A cancelled fit returns [`Cancelled`] and writes nothing; one left to
/// run comes to a verdict and writes its tables, too late to cancel then.
#[test]
fn a_cancelled_fit_writes_nothing() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().join("out");
    let mut record = crate::run_record::RunRecord::start("cluster", "unused");
    let p = Progress::new(2);
    p.cancel();
    let out = fit_recorded(&p, &dir, &mut record);
    assert!(out.err().unwrap().is::<Cancelled>());
    assert!(!dir.exists(), "nothing written");

    let p = Progress::new(2);
    fit_recorded(&p, &dir, &mut record).unwrap();
    assert!(p.verdict().is_some());
    assert!(!p.cancel(), "too late to cancel");
    assert!(std::fs::read_dir(&dir).unwrap().count() > 0);
}

/// A cancel that comes while the fit reads stops it at the next step,
/// before the next input is read: here before the first, which does not
/// exist.
#[test]
fn a_cancel_while_reading_stops_before_the_next_read() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path().join("out");
    let missing = tmp.path().join("missing.zarr");
    let args = super::parse_args(&["cluster", missing.to_str().unwrap()]);
    let mut record = crate::run_record::RunRecord::start("cluster", "unused");
    let p = Progress::new(2);
    p.cancel();
    let out = crate::cluster::cluster(&args, &dir, &mut record, &p);
    assert!(out.err().unwrap().is::<Cancelled>());
    assert!(!dir.exists());
}

/// A flag of the MCMC fits is refused before any input is read or anything
/// is written: here before the input, which does not exist, is looked for.
#[test]
fn an_mcmc_flag_is_refused_before_anything_is_read() {
    let tmp = tempfile::tempdir().unwrap();
    let missing = tmp.path().join("missing.zarr");
    let out = tmp.path().join("out");
    let args = super::parse_args(&[
        "cluster",
        missing.to_str().unwrap(),
        "--sweeps",
        "100",
        "-o",
        out.to_str().unwrap(),
    ]);
    let refused = crate::cluster::run_cluster(&args).unwrap_err().to_string();
    assert!(refused.contains("--sweeps"), "{refused}");
    assert!(!out.exists() && !missing.exists());
}
