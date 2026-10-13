#[cfg(test)]
mod tests;
mod tui;

use super::output_dir;
use super::site_io::*;
use crate::common::*;
use crate::data::dna::{Dna, DnaBaseCount};
use crate::data::dna_stat_map::DnaBaseFreqMap;
use crate::data::util_htslib;

use genomic_data::bed::Bed;
use genomic_data::sam::Strand;

use std::fs::File;
use std::io::{BufWriter, Write};

#[derive(clap::ValueEnum, Clone, Debug)]
pub enum PwmSource {
    Reads,
    Reference,
}

#[derive(Args, Debug)]
pub struct ScanPwmArgs {
    #[arg(
        value_delimiter = ',',
        value_name = "SITES|DIR|BAM",
        help = "Site parquets or faba output directories, and BAM files for --source reads",
        long_help = "Site-level parquet file (from dartseq, atoi or apa output).\n\
                     Or a faba output directory (from `faba run`, a producer, or `faba qc`),\n\
                     or a run record in it (`*.run.json`, `pipeline_summary.json`): its\n\
                     `{modality}_sites.parquet` of --modality is read. Several are read\n\
                     together, as one set of sites. A `.bam` or `.cram` is a BAM file,\n\
                     read with --source reads. Without sites (here or -s), a browser\n\
                     asks: mark site tables, or choose an output directory."
    )]
    input: Vec<Box<str>>,

    #[arg(
        short = 's',
        long = "sites",
        help = "Site-level parquet file, as the positional SITES"
    )]
    sites: Option<Box<str>>,

    #[arg(
        long = "modality",
        value_name = "m6a|atoi",
        help = "Which site table to take from an output directory (default: m6a, else atoi)"
    )]
    modality: Option<Box<str>>,

    #[arg(
        short = 'f',
        long = "genome",
        help = "Reference genome FASTA for --source reference (default: the one recorded next to --sites)",
        long_help = "Reference genome FASTA, read with --source reference.\n\
                     Without it, the FASTA recorded in the `*.run.json` next to --sites is used."
    )]
    genome_file: Option<Box<str>>,

    #[arg(
        long,
        value_enum,
        default_value = "reference",
        help = "Base frequency source: reads (from BAM) or reference (from FASTA)"
    )]
    source: PwmSource,

    #[arg(
        short = 'w',
        long,
        default_value_t = 10,
        help = "Half-window size: collect +/- this many bp around each site"
    )]
    window: i64,

    #[arg(
        long,
        alias = "max-threads",
        default_value_t = 16,
        help = "Number of threads"
    )]
    threads: usize,

    #[arg(
        short,
        long,
        help = "Output file path, TSV or .gz (default with --batch-process: `{sites}.pwm.tsv` here)",
        long_help = "Output file path (TSV, or .gz for gzipped). In the view, written only\n\
                     when given; with --batch-process, without it, `{sites}.pwm.tsv` in the\n\
                     current directory, after the site table's name."
    )]
    output: Option<Box<str>>,

    #[arg(
        long = "batch-process",
        visible_alias = "batch-mode",
        default_value_t = false,
        help = "Write the PWM and stop, without showing the logo full screen",
        long_help = "Write the PWM and stop, without showing the logo full screen.\n\
                     By default `faba pwm` writes it (when -o is given) and then shows the\n\
                     sequence logo in a view (when there is a terminal)."
    )]
    batch_process: bool,

    /// The view is the default now; kept so older command lines still run.
    #[arg(
        short = 'I',
        long = "interactive",
        hide = true,
        default_value_t = false
    )]
    interactive: bool,
}

/// Swap A<->T and G<->C counts to complement a DnaBaseCount.
fn complement_base_count(count: &DnaBaseCount) -> DnaBaseCount {
    let mut out = DnaBaseCount::new();
    out.add(Some(&Dna::A), count.count_t());
    out.add(Some(&Dna::T), count.count_a());
    out.add(Some(&Dna::G), count.count_c());
    out.add(Some(&Dna::C), count.count_g());
    out
}

/// Collect PWM from reference FASTA sequence around each site.
fn collect_from_reference(
    sites: &[GenomicSite],
    genome_file: &str,
    window: i64,
) -> anyhow::Result<Vec<DnaBaseCount>> {
    let width = (2 * window + 1) as usize;
    let mut pwm: Vec<DnaBaseCount> = (0..width).map(|_| DnaBaseCount::new()).collect();

    let faidx = util_htslib::load_fasta_index(genome_file)?;

    for site in sites {
        let start = site.position - window;
        let end = site.position + window; // inclusive
        let seq = util_htslib::fetch_reference_seq(&faidx, &site.chr, start, end)?;

        if let Some(bases) = seq {
            if bases.len() != width {
                continue; // skip sites near chromosome boundaries
            }
            match site.strand {
                Strand::Forward => {
                    for (j, base) in bases.iter().enumerate() {
                        pwm[j].add(Some(base), 1);
                    }
                }
                Strand::Backward => {
                    // Reverse complement: reverse order and complement each base
                    for (j, base) in bases.iter().rev().enumerate() {
                        let comp = match base {
                            Dna::A => Dna::T,
                            Dna::T => Dna::A,
                            Dna::G => Dna::C,
                            Dna::C => Dna::G,
                        };
                        pwm[j].add(Some(&comp), 1);
                    }
                }
            }
        }
    }

    Ok(pwm)
}

/// Collect PWM from BAM read base frequencies around each site.
/// `quiet` drops the progress bar, for when a view's spinner shows instead.
fn collect_from_reads(
    sites: &[GenomicSite],
    bam_files: &[Box<str>],
    window: i64,
    quiet: bool,
) -> anyhow::Result<Vec<DnaBaseCount>> {
    let width = (2 * window + 1) as usize;

    // Ensure BAM indexes exist
    for bam_file in bam_files {
        util_htslib::check_bam_index(bam_file, None)?;
    }

    let empty = || (0..width).map(|_| DnaBaseCount::new()).collect::<Vec<_>>();
    let pwm: Vec<DnaBaseCount> = sites
        .par_iter()
        .progress_with(if quiet {
            indicatif::ProgressBar::hidden()
        } else {
            new_progress_bar(sites.len() as u64)
        })
        .map(|site| -> anyhow::Result<Vec<DnaBaseCount>> {
            let mut local_pwm = empty();

            let bed = Bed {
                chr: site.chr.clone(),
                start: (site.position - window).max(0),
                stop: site.position + window + 1,
            };

            let mut freq_map = DnaBaseFreqMap::new();
            for bam_file in bam_files {
                freq_map.update_from_region(bam_file, &bed).map_err(|e| {
                    anyhow::anyhow!("{bam_file}: {}:{}: {e}", bed.chr, site.position)
                })?;
            }

            if let Some(counts) = freq_map.marginal_frequency_map() {
                for (j, p) in ((site.position - window)..=(site.position + window)).enumerate() {
                    if let Some(count) = counts.get(&p) {
                        match site.strand {
                            Strand::Forward => {
                                local_pwm[j] += count;
                            }
                            Strand::Backward => {
                                // Reverse offset and complement
                                let rev_j = width - 1 - j;
                                let comp = complement_base_count(count);
                                local_pwm[rev_j] += &comp;
                            }
                        }
                    }
                }
            }

            Ok(local_pwm)
        })
        .try_reduce(empty, |mut acc, local| {
            for (a, b) in acc.iter_mut().zip(local.iter()) {
                *a += b;
            }
            Ok(acc)
        })?;

    Ok(pwm)
}

/// Write PWM table as TSV.
fn write_pwm(pwm: &[DnaBaseCount], window: i64, output: &str) -> anyhow::Result<()> {
    let file = File::create(output)?;
    let mut writer: Box<dyn Write> = if output.ends_with(".gz") {
        Box::new(BufWriter::new(flate2::write::GzEncoder::new(
            file,
            flate2::Compression::default(),
        )))
    } else {
        Box::new(BufWriter::new(file))
    };

    writeln!(
        writer,
        "position\tcount_A\tcount_T\tcount_G\tcount_C\ttotal\tfreq_A\tfreq_T\tfreq_G\tfreq_C"
    )?;

    for (j, count) in pwm.iter().enumerate() {
        let rel_pos = j as i64 - window;
        let a = count.count_a();
        let t = count.count_t();
        let g = count.count_g();
        let c = count.count_c();
        let total = a + t + g + c;
        let tot_f = total as f32;
        if total > 0 {
            writeln!(
                writer,
                "{}\t{}\t{}\t{}\t{}\t{}\t{:.6}\t{:.6}\t{:.6}\t{:.6}",
                rel_pos,
                a,
                t,
                g,
                c,
                total,
                a as f32 / tot_f,
                t as f32 / tot_f,
                g as f32 / tot_f,
                c as f32 / tot_f,
            )?;
        } else {
            writeln!(writer, "{}\t0\t0\t0\t0\t0\t0\t0\t0\t0", rel_pos)?;
        }
    }

    Ok(())
}

/// The genome to read every table's sites against: `explicit`, else the one
/// recorded next to each table, which must be the same for all.
fn genome_of(explicit: Option<&str>, site_files: &[Box<str>]) -> anyhow::Result<Box<str>> {
    if let Some(genome) = explicit {
        return Ok(genome.into());
    }
    let recorded = |f: &str| {
        crate::run_record::explicit_or_recorded(None, f, "genome", "genome").ok_or_else(|| {
            anyhow::anyhow!(
                "--genome is required when --source reference: \
                 no run record next to {f} names a genome"
            )
        })
    };
    let first = recorded(&site_files[0])?;
    for f in &site_files[1..] {
        let this = recorded(f)?;
        anyhow::ensure!(
            this == first,
            "the site tables were called against different genomes ({first}, {this}); \
             name one with --genome"
        );
    }
    Ok(first)
}

/// Whether `path` names an alignment file (BAM, CRAM or SAM) rather than
/// sites.
fn is_bam(path: &str) -> bool {
    let lower = path.to_ascii_lowercase();
    [".bam", ".cram", ".sam"]
        .iter()
        .any(|ext| lower.ends_with(ext))
}

/// Site tables, and BAM files.
type Inputs = (Vec<Box<str>>, Vec<Box<str>>);

impl ScanPwmArgs {
    /// The site tables and the BAM files: `-s`, then the positional inputs
    /// that are not BAMs, each of which may name an output directory (or a
    /// run record in it) holding one; with neither, those chosen in a
    /// browser. `None` when that is cancelled.
    fn inputs(&self) -> anyhow::Result<Option<Inputs>> {
        let (bams, given): (Vec<Box<str>>, Vec<Box<str>>) =
            self.input.iter().cloned().partition(|p| is_bam(p));
        let given = self.sites.iter().cloned().chain(given).collect();
        let modality = self.modality.as_deref();
        let ask = output_dir::site_tables("pwm", self.batch_process, given, modality);
        Ok(ask?.map(|sites| (sites, bams)))
    }
}

pub fn run_scan_pwm(args: &ScanPwmArgs) -> anyhow::Result<()> {
    rayon::ThreadPoolBuilder::new()
        .num_threads(args.threads)
        .build_global()
        .ok(); // ignore if already set

    let Some((site_files, bam_files)) = args.inputs()? else {
        return Ok(());
    };
    if !bam_files.is_empty() && matches!(args.source, PwmSource::Reference) {
        log::warn!(
            "{} read only with --source reads; the reference gives the bases",
            bam_files.join(", ")
        );
    }
    // A genome missing, or differing between tables, fails before the reading.
    let genome = match args.source {
        PwmSource::Reference => Some(genome_of(args.genome_file.as_deref(), &site_files)?),
        PwmSource::Reads => {
            anyhow::ensure!(
                !bam_files.is_empty(),
                "BAM file(s) required when --source reads"
            );
            None
        }
    };
    let batch = crate::tui::batch(args.batch_process, "pwm");
    let what = format!("collecting bases around {} site table(s)", site_files.len());
    let pwm = crate::tui::busy(!batch, "pwm", &what, || {
        collect(args, &site_files, genome.as_deref(), &bam_files, !batch)
    })?;

    let output = args
        .output
        .clone()
        .or_else(|| batch.then(|| output_dir::default_output(&site_files, ".pwm.tsv")));
    if let Some(output) = &output {
        write_pwm(&pwm, args.window, output)?;
        info!("wrote PWM to {output}");
    }
    if !batch {
        tui::show_pwm(&output_dir::title(&site_files), &pwm, args.window)?;
    }

    Ok(())
}

/// The base counts around the sites of `site_files`, as one set, from
/// `genome` when given, else from `bam_files`.
fn collect(
    args: &ScanPwmArgs,
    site_files: &[Box<str>],
    genome: Option<&str>,
    bam_files: &[Box<str>],
    quiet: bool,
) -> anyhow::Result<Vec<DnaBaseCount>> {
    // The tables are independent: read them together.
    let tables = site_files
        .par_iter()
        .map(|f| read_sites(f))
        .collect::<anyhow::Result<Vec<_>>>()?;
    let sites: Vec<GenomicSite> = tables.into_iter().flatten().collect();
    if site_files.len() > 1 {
        info!(
            "{} sites from {} tables, together",
            sites.len(),
            site_files.len()
        );
    }
    match genome {
        Some(genome) => collect_from_reference(&sites, genome, args.window),
        None => collect_from_reads(&sites, bam_files, args.window, quiet),
    }
}
