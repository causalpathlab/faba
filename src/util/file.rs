use anyhow::{self};

use rust_htslib::bgzf;
use std::fs::File;
use std::io::{BufRead, BufReader}; // TODO: BufWriter, Stdout, Write
use std::path::Path;

///
/// Read every line of the input_file into memory
///
pub fn read_lines(input_file: &str) -> anyhow::Result<Vec<Box<str>>> {
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
