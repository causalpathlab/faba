use bio::io::gff;

/// Parse a GFF line to a record
///
/// https://en.wikipedia.org/wiki/General_feature_format
///
pub fn parse(line: Box<str>) -> Option<Box<gff::Record>> {
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
