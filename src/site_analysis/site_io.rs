use genomic_data::sam::Strand;
use log::info;
use rustc_hash::FxHashSet;

use parquet::file::reader::{FileReader, SerializedFileReader};
use parquet::record::RowAccessor;

use std::fs::File;

pub struct GenomicSite {
    pub chr: Box<str>,
    pub position: i64,
    pub strand: Strand,
}

/// A site table's `strand` column: `+` is forward, anything else reverse.
pub fn parse_strand(s: &str) -> Strand {
    if s == "+" {
        Strand::Forward
    } else {
        Strand::Backward
    }
}

/// Each of `sites`' reads in `site_file`, as `(converted, unconverted)`:
/// from the site table's `converted` and `coverage` columns, summed over its
/// rows at the site. `None` when the table has no such columns (an APA
/// table, or one written before them).
pub fn read_site_reads(
    site_file: &str,
    sites: &[GenomicSite],
) -> anyhow::Result<Option<Vec<(f64, f64)>>> {
    let field_names = legume_numeric::matrix::parquet::peek_parquet_field_names(site_file)?;
    let idx = |name: &str| field_names.iter().position(|f| f.as_ref() == name);
    let pos = ["m6a_pos", "genomic_alpha", "primary_pos"]
        .iter()
        .find_map(|n| idx(n));
    let (Some(chr), Some(pos), Some(coverage), Some(converted)) =
        (idx("chr"), pos, idx("coverage"), idx("converted"))
    else {
        return Ok(None);
    };
    let reader = SerializedFileReader::new(File::open(site_file)?)?;
    // Only the four columns, in the file's order, so no other is decoded.
    let wanted = [chr, pos, coverage, converted];
    let schema = reader.metadata().file_metadata().schema();
    let fields: Vec<_> = schema
        .get_fields()
        .iter()
        .enumerate()
        .filter(|(i, _)| wanted.contains(i))
        .map(|(_, f)| f.clone())
        .collect();
    let at = |i: usize| wanted.iter().filter(|&&w| w < i).count();
    let (chr, pos, coverage, converted) = (at(chr), at(pos), at(coverage), at(converted));
    let projection = parquet::schema::types::Type::group_type_builder(schema.name())
        .with_fields(fields)
        .build()?;
    let count = |row: &parquet::record::Row, i: usize| -> anyhow::Result<f64> {
        Ok(match row.get_ulong(i) {
            Ok(v) => v as f64,
            Err(_) => row.get_long(i)? as f64,
        })
    };
    // By chromosome, then position: a row allocates only for a new
    // chromosome, and a site looks itself up by borrowing its name.
    type Sums = rustc_hash::FxHashMap<i64, (f64, f64)>;
    let mut reads: rustc_hash::FxHashMap<Box<str>, Sums> = Default::default();
    for record in reader.get_row_iter(Some(projection))? {
        let row = record?;
        let name = row.get_string(chr)?.as_str();
        let (c, n) = (count(&row, converted)?, count(&row, coverage)?);
        if !reads.contains_key(name) {
            reads.insert(name.into(), Sums::default());
        }
        let e = reads
            .get_mut(name)
            .and_then(|by_pos| Some(by_pos.entry(row.get_long(pos).ok()?).or_default()));
        if let Some(e) = e {
            e.0 += c;
            e.1 += n;
        }
    }
    Ok(Some(
        sites
            .iter()
            .map(|s| {
                let (c, n) = reads
                    .get(s.chr.as_ref())
                    .and_then(|by_pos| by_pos.get(&s.position))
                    .copied()
                    .unwrap_or_default();
                (c, (n - c).max(0.0))
            })
            .collect(),
    ))
}

/// Read sites from a parquet file, auto-detecting dart vs apa vs atoi format.
pub fn read_sites(site_file: &str) -> anyhow::Result<Vec<GenomicSite>> {
    let field_names = legume_numeric::matrix::parquet::peek_parquet_field_names(site_file)?;
    let has_m6a_pos = field_names.iter().any(|f| f.as_ref() == "m6a_pos");
    let has_genomic_alpha = field_names.iter().any(|f| f.as_ref() == "genomic_alpha");
    let has_primary_pos = field_names.iter().any(|f| f.as_ref() == "primary_pos");

    let file = File::open(site_file)?;
    let reader = SerializedFileReader::new(file)?;
    let row_iter = reader.get_row_iter(None)?;

    // Find column indices
    let chr_idx = field_names
        .iter()
        .position(|f| f.as_ref() == "chr")
        .ok_or_else(|| anyhow::anyhow!("missing 'chr' column in {}", site_file))?;

    let mut sites = Vec::new();

    if has_m6a_pos {
        // Dart format: chr, m6a_pos, strand
        let pos_idx = field_names
            .iter()
            .position(|f| f.as_ref() == "m6a_pos")
            .unwrap();
        let strand_idx = field_names
            .iter()
            .position(|f| f.as_ref() == "strand")
            .ok_or_else(|| anyhow::anyhow!("missing 'strand' column in {}", site_file))?;

        for record in row_iter {
            let row = record?;
            let chr: Box<str> = row.get_string(chr_idx)?.clone().into_boxed_str();
            let position = row.get_long(pos_idx)?;
            let strand = parse_strand(row.get_string(strand_idx)?);
            sites.push(GenomicSite {
                chr,
                position,
                strand,
            });
        }
    } else if has_genomic_alpha {
        // APA format: chr, genomic_alpha (no strand column)
        let pos_idx = field_names
            .iter()
            .position(|f| f.as_ref() == "genomic_alpha")
            .unwrap();

        for record in row_iter {
            let row = record?;
            let chr: Box<str> = row.get_string(chr_idx)?.clone().into_boxed_str();
            let position = row.get_long(pos_idx)?;
            sites.push(GenomicSite {
                chr,
                position,
                strand: Strand::Forward,
            });
        }
    } else if has_primary_pos {
        // ATOI format: chr, primary_pos, strand
        let pos_idx = field_names
            .iter()
            .position(|f| f.as_ref() == "primary_pos")
            .unwrap();
        let strand_idx = field_names
            .iter()
            .position(|f| f.as_ref() == "strand")
            .ok_or_else(|| anyhow::anyhow!("missing 'strand' column in {}", site_file))?;

        for record in row_iter {
            let row = record?;
            let chr: Box<str> = row.get_string(chr_idx)?.clone().into_boxed_str();
            let position = row.get_long(pos_idx)?;
            let strand = parse_strand(row.get_string(strand_idx)?);
            sites.push(GenomicSite {
                chr,
                position,
                strand,
            });
        }
    } else {
        return Err(anyhow::anyhow!(
            "unrecognized parquet format: expected 'm6a_pos' (dart), 'genomic_alpha' (apa), or 'primary_pos' (atoi) column"
        ));
    }

    // Deduplicate by (chr, position)
    let mut seen = FxHashSet::default();
    sites.retain(|s| seen.insert((s.chr.clone(), s.position)));

    info!("loaded {} unique sites from {}", sites.len(), site_file);
    Ok(sites)
}

#[cfg(test)]
#[path = "tests/site_io.rs"]
mod tests;
