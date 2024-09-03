use anyhow::{self};
use bio::bio_types::strand::Strand;
use bio::io::gff::{self};

use clap::{Args, Parser, Subcommand};
use rust_htslib::bgzf;
use std::cmp::{max, min};
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Stdout, Write};
use std::path::Path;
use std::sync::{Arc, Mutex};

use rayon::prelude::*;

#[derive(Parser)]
#[command(version, about, long_about=None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    commands: Commands,
}

#[derive(Subcommand)]
enum Commands {
    /// Freq GFF file
    Freq(FreqArgs),
    /// Subset GFF file
    Subset(SubsetArgs),
}

#[derive(Args)]
struct FreqArgs {
    /// GFF file
    file: Option<Box<str>>,
}

#[derive(Args)]
struct SubsetArgs {
    /// GFF file
    file: Option<Box<str>>,

    /// genomic feature name
    #[arg(short, long)]
    feature: Option<Box<str>>,

    /// genomic region
    #[arg(short, long)]
    seqname: Option<Box<str>>,

    /// genomic region
    #[arg(short, long)]
    lb: Option<u64>,

    /// genomic region
    #[arg(short, long)]
    ub: Option<u64>,
}

fn main() -> anyhow::Result<()> {
    let cli = Cli::parse();

    match &cli.commands {
        Commands::Freq(args) => {
            let file_name = args.file.as_ref().expect("GFF file name");
            freq_main(file_name);
        }

        Commands::Subset(args) => {
            subset_main(
                args.file.as_ref().expect("GFF file name"),
                args.feature.as_ref().expect(""),
                args.seqname.clone(),
                args.lb.clone(),
                args.ub.clone(),
            )?;
        }
    }

    Ok(())
}

/// main for subset
fn subset_main(
    input_file: &str,
    select_feature: &str,
    seqname: Option<Box<str>>,
    lb: Option<u64>,
    ub: Option<u64>,
) -> std::io::Result<()> {
    type Data = Vec<Box<gff::Record>>;
    let lines = read_gff_lines(input_file).expect("failed to read lines from the GFF file");
    let records: Data = lines.into_par_iter().filter_map(parse_gff_record).collect();

    let to_bed_str = |x: &gff::Record| -> Box<str> {
        format!(
            "{}\t{}\t{}\t{}",
            x.seqname(),
            x.start(),
            x.end(),
            x.strand().map_or(".", |st| match st {
                Strand::Forward => "+",
                Strand::Reverse => "-",
                _ => ".",
            })
        )
        .into_boxed_str()
    };

    let output: Vec<Box<str>> = match &seqname {
        Some(seq) => match (&lb, &ub) {
            (Some(lb), Some(ub)) => records
                .into_iter()
                .par_bridge()
                .filter(|x| -> bool {
                    if x.seqname() != seq.to_string() {
                        return false;
                    }
                    let (a, b) = (x.start(), x.end());
                    //               a..b
                    //      lb....ub
                    // a..b
                    //
                    if b < lb || a > ub {
                        return false;
                    };
                    true
                })
                .filter(|x| x.feature_type() == select_feature)
                .map(|x| -> Box<str> { to_bed_str(&x) })
                .collect(),
            _ => records
                .into_iter()
                .par_bridge()
                .filter(|x| x.seqname() == seq.to_string())
                .filter(|x| x.feature_type() == select_feature)
                .map(|x| -> Box<str> { to_bed_str(&x) })
                .collect(),
        },
        _ => records
            .into_iter()
            .par_bridge()
            .filter(|x| x.feature_type() == select_feature)
            .map(|x| -> Box<str> { to_bed_str(&x) })
            .collect(),
    };

    let mut buf = BufWriter::new(std::io::stdout());

    for line in output {
        buf.write_all(line.as_bytes())?;
        buf.write_all(b"\n")?;
    }

    // TODO: output

    todo!("take the CLI argument struct");

    todo!("output attributes");

    // do something
    Ok(())
}

/// main for feature frequency
fn freq_main(input_file: &str) {
    let lines = read_gff_lines(input_file).expect("failed to read lines from the GFF file");
    type Data = Vec<Box<gff::Record>>;
    let records: Data = lines.into_par_iter().filter_map(parse_gff_record).collect();
    let arc_freq = Arc::new(Mutex::new(FreqMap::new()));

    records.into_iter().par_bridge().for_each(|x| {
        let mut freq = arc_freq.lock().unwrap();
        freq.push(x.as_ref());
    });

    println!("{}", arc_freq.lock().unwrap());
}

struct FreqMap {
    feature_map: HashMap<Box<str>, u64>,
    seqname_map: HashMap<Box<str>, u64>,
}

impl FreqMap {
    fn push(&mut self, record: &gff::Record) {
        // feature type
        let n = self
            .feature_map
            .entry(record.feature_type().into())
            .or_insert(0);
        *n += 1;
        // sequence name
        let n = self.seqname_map.entry(record.seqname().into()).or_insert(0);
        *n += 1;
    }

    fn new() -> Self {
        Self {
            feature_map: HashMap::new(),
            seqname_map: HashMap::new(),
        }
    }
}

impl std::fmt::Display for FreqMap {
    fn fmt(&self, fmt: &mut std::fmt::Formatter) -> std::fmt::Result {
        let (sep, nl) = ("\t", "\n");

        // fmt.write_fmt(fmt)
        fmt.write_str("# Feature Frequency\n")?;
        for (k, n) in self.feature_map.iter() {
            fmt.write_str(k.as_ref())?;
            fmt.write_str(sep)?;
            fmt.write_str(&format!("{}", n).to_string())?;
            fmt.write_str(nl)?;
        }
        fmt.write_str(nl)?;
        fmt.write_str("# Sequence Name Frequency\n")?;
        for (k, n) in self.seqname_map.iter() {
            fmt.write_str(k.as_ref())?;
            fmt.write_str(sep)?;
            fmt.write_str(&format!("{}", n).to_string())?;
            fmt.write_str(nl)?;
        }
        Ok(())
    }
}

///
/// Read every line of the input_file into memory
///
fn read_gff_lines(input_file: &str) -> anyhow::Result<Vec<Box<str>>> {
    // Any better solution without using Box?
    let buf: Box<dyn BufRead> = match Path::new(input_file).extension().and_then(|x| x.to_str()) {
        Some("gz") | Some("bgz") => {
            let _file = bgzf::Reader::from_path(input_file)?;
            Box::new(BufReader::new(_file))
        }

        _ => {
            let _file = File::open(input_file)?;
            Box::new(BufReader::new(_file))
        }
    };

    let mut lines = vec![];
    for x in buf.lines() {
        lines.push(x?.into_boxed_str());
    }
    Ok(lines)
}

/// Parse a GFF line to a record
///
/// https://en.wikipedia.org/wiki/General_feature_format
///
fn parse_gff_record(line: Box<str>) -> Option<Box<gff::Record>> {
    const SEP: char = '\t';
    const SEP_ATTR: char = ':';
    const NUM_FIELDS: usize = 9;

    // we can access immutatble &str
    // let x = line.as_ref();

    let words: Vec<_> = line.split(SEP).collect();

    if words.len() == NUM_FIELDS {
        let mut rec = gff::Record::new();
        *rec.seqname_mut() = words[0].to_string();
        *rec.source_mut() = words[1].to_string();
        *rec.feature_type_mut() = words[2].to_string();
        *rec.start_mut() = words[3].parse().unwrap_or(0);
        *rec.end_mut() = words[4].parse().unwrap_or(0);
        *rec.score_mut() = words[5].to_string();
        *rec.strand_mut() = words[6].to_string();
        *rec.phase_mut() = match words[7] {
            "." => gff::Phase::default(),
            _ => gff::Phase::from(words[7].parse().unwrap_or(0u8)),
        };

        for z in words[8].split_whitespace() {
            let kv: Vec<&str> = z.split(SEP_ATTR).collect();
            if kv.len() == 2 {
                rec.attributes_mut()
                    .insert(kv[0].to_string(), kv[1].to_string());
            }
        }
        Some(Box::new(rec))
    } else {
        None
    }
}

// #[derive(Debug, clap::ValueEnum, Clone, Copy)]
// enum FeatureType {
//     StartCodon,
//     StopCodon,
//     CDS,
//     Exon,
//     Utr3,
//     Utr5,
//     Utr,
// }

// impl

// impl FeatureType {
//     fn as_str(&self) -> &str {
//         match self {
//             FeatureType::StartCodon => "start_codon",
//             FeatureType::StopCodon => "stop_codon",
//             FeatureType::CDS => "CDS",
//             FeatureType::Exon => "exon",
//             FeatureType::Utr3 => "three_prime_UTR",
//             FeatureType::Utr5 => "five_prime_UTR",
//             FeatureType::Utr => "UTR",
//         }
//     }
// }
