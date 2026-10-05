use super::args::*;
use super::steps::*;
use crate::common::*;
use crate::quant::check_all_bam_indices;
use crate::run_record::{recorded, RunRecord};

use log::info;
use rayon::ThreadPoolBuilder;

/// The full set of samples to QUANTIFY in every modality: the positional
/// (signal/WT) BAMs together with the `--control-bam` (MUT/YTHmut) BAMs,
/// deduplicated (a BAM may legitimately be listed in both roles). The WT-vs-MUT
/// split is used ONLY for the m6A discovery contrast (step 4); SNP, gene counts,
/// ATOI, APA and m6A per-cell matrices are all produced for EVERY one of these
/// samples, so the control background is fully quantified — not merely consumed
/// as an m6A reference. This also freezes a cell set for each control BAM in
/// step 1, so the control m6A matrices reuse it instead of the ambient superset.
pub(super) fn all_quant_bam_files(args: &PipelineArgs) -> Vec<Box<str>> {
    let (files, dropped) = unique_bam_files(
        args.bam_files
            .iter()
            .chain(args.control_bam_files.iter())
            .cloned(),
    );
    if dropped > 0 {
        log::warn!(
            "{dropped} BAM file(s) listed both positionally and in --control-bam; \
             quantifying each once to avoid double counting"
        );
    }
    files
}

pub fn run_pipeline(args: &PipelineArgs) -> anyhow::Result<()> {
    // 0. Setup
    info!("faba pipeline: unified RNA-seq analysis");
    ThreadPoolBuilder::new()
        .num_threads(args.max_threads)
        .build_global()?;
    std::fs::create_dir_all(args.out())?;
    let summary = step_record(args, "run");

    // Validate inputs
    check_all_bam_indices(&args.bam_files)?;
    check_all_bam_indices(&args.control_bam_files)?;

    let n_steps = 5;

    // Step 0: SNP genotyping (de novo discovery + optional known sites). Its
    // outputs stand alone; no later step consumes them as a mask.
    if !args.skip_snp {
        info!("Step 0/{}: SNP genotyping", n_steps);
        match recorded(step_record(args, "snp"), || run_snp_step(args)) {
            Ok(()) => info!("SNP complete"),
            Err(e) => log::warn!("SNP step failed: {}. Continuing.", e),
        }
    } else {
        info!("Step 0/{}: SKIPPED (--skip-snp)", n_steps);
    }

    // Step 1: gene counting and cell calling
    let gene_count_qc = if !args.skip_count {
        info!("Step 1/{}: gene counting and cell calling", n_steps);
        recorded(step_record(args, "count"), || run_gene_counting_step(args))?
    } else {
        info!("Step 1/{}: SKIPPED (--skip-count)", n_steps);
        None
    };

    // Step 2: per-cell read depth -- opt-in, and independent: it consumes no
    // mask and produces none, and nothing downstream reads `{batch}_depth`. It
    // sits here, directly after gene counting, only because that is where the
    // called-cell axis becomes available -- the depth matrix shares its columns
    // with every other modality rather than inventing its own. Being independent
    // it can also fail without costing anything that follows.
    if args.depth_resolution_kb.is_some() {
        info!("Step 2/{}: per-cell read depth", n_steps);
        match recorded(step_record(args, "depth"), || {
            run_read_depth_step(args, &gene_count_qc)
        }) {
            Ok(_) => info!("Read depth complete"),
            Err(e) => log::warn!("Read depth step failed: {}", e),
        }
    }

    // Step 3: ATOI Detection
    if !args.skip_atoi {
        info!("Step 3/{}: ATOI detection", n_steps);
        match recorded(step_record(args, "atoi"), || {
            run_atoi_step(args, &gene_count_qc)
        }) {
            Ok(n_sites) => info!("ATOI complete: {} putative sites", n_sites),
            Err(e) => log::warn!("ATOI step failed: {}. Continuing.", e),
        }
    } else {
        info!("Step 3/{}: SKIPPED (--skip-atoi)", n_steps);
    }

    // Step 4: m6A (DART) detection — WT-vs-MUT contrast at motif Cs (signal arm =
    // positional BAMs minus --control-bam, tested against the pooled control).
    // It runs BEFORE the heavy APA EM so the fast modalities all finish first.
    // Requires a control; skipped (not failed) when none is supplied.
    if args.skip_m6a || args.control_bam_files.is_empty() {
        if args.skip_m6a {
            info!("Step 4/{}: SKIPPED (--skip-m6a)", n_steps);
        } else {
            info!(
                "Step 4/{}: SKIPPED (m6A needs --control-bam for the WT-vs-MUT contrast)",
                n_steps
            );
        }
    } else {
        info!("Step 4/{}: m6A detection", n_steps);
        match recorded(step_record(args, "dartseq"), || {
            run_dart_step(args, &gene_count_qc)
        }) {
            Ok(_) => info!("m6A complete"),
            Err(e) => log::warn!("m6A step failed: {}", e),
        }
    }

    // Step 5: APA analysis — the heavy SCAPE EM, run LAST so it never blocks the
    // fast modalities (genes / depth / ATOI / m6A) that downstream work needs first.
    if !args.skip_apa {
        info!("Step 5/{}: APA analysis", n_steps);
        match recorded(step_record(args, "apa"), || {
            run_apa_step(args, &gene_count_qc)
        }) {
            Ok(_) => info!("APA complete"),
            Err(e) => log::warn!("APA step failed: {}", e),
        }
    } else {
        info!("Step 5/{}: SKIPPED (--skip-apa)", n_steps);
    }

    summary.finish(&Ok(()));
    info!("Pipeline complete! Results in: {}", args.out());
    Ok(())
}

/// A record for one step of the run, with the inputs that step reads, under
/// the name the standalone subcommand writes, so a step's outputs can be
/// traced to it. As `all`, the whole run: every input, written to
/// `pipeline_summary.json`.
///
/// Every option is recorded, and it is [`PipelineArgs`] itself that is
/// serialized rather than a hand-kept list: a run is defined as much by the
/// defaults it did not override as by the flags it passed, and faba's options
/// have changed between builds (`--cluster-resolution` defaulted to 0.5, then
/// to 0, and is now gone; `--n-bootstrap` existed and then did not), with a
/// version history that is not monotonic. A new option appears here the moment
/// it is added to `PipelineArgs`, with no second list to keep in sync.
fn step_record(args: &PipelineArgs, job: &str) -> RunRecord {
    let record = RunRecord::start(job, args.out()).options(args);
    let (gff, genome) = (Some(args.gff()), Some(args.genome()));
    match job {
        "run" => record
            .file_name("pipeline_summary.json")
            .inputs("bam", &args.bam_files)
            .inputs("control_bam", &args.control_bam_files)
            .input("gff", gff)
            .input("genome", genome)
            .input("known_snps", args.known_snps.as_deref()),
        "dartseq" => {
            let signal: Vec<&str> = args
                .bam_files
                .iter()
                .filter(|b| !args.control_bam_files.contains(b))
                .map(|b| &**b)
                .collect();
            record
                .inputs("bam", &signal)
                .inputs("control_bam", &args.control_bam_files)
                .input("gff", gff)
                .input("genome", genome)
        }
        _ => {
            let record = record.inputs("bam", &all_quant_bam_files(args));
            match job {
                "snp" => record
                    .input("genome", genome)
                    .input("gff", gff)
                    .input("known_snps", args.known_snps.as_deref()),
                "atoi" => record.input("gff", gff).input("genome", genome),
                "count" | "apa" => record.input("gff", gff),
                _ => record,
            }
        }
    }
}

#[cfg(test)]
mod tests;

/// `faba run`: straight through under `--batch-process`, else the setup view.
pub fn run_or_view(args: &PipelineArgs, cli: clap::Command) -> anyhow::Result<()> {
    if args.batch_process {
        args.check_batch()?;
        return run_pipeline(args);
    }
    anyhow::ensure!(
        data_beans::interactive::tui_available(),
        "`faba run` sets up the run in a full-screen view, which needs stdin and stdout \
         on a terminal; pass --batch-process with the BAMs, -g, -f and -o to run without it"
    );
    let _ = cli;
    anyhow::bail!("the setup view is not built yet")
}
