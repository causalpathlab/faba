use super::*;
use genomic_data::transcript::build_transcript_models;

//////////////
// Fixtures  //
//////////////

const CHR: &str = "chr1";

fn rec(
    tx: &str,
    gene: &str,
    feature_type: FeatureType,
    start: i64,
    stop: i64,
    strand: Strand,
) -> GffRecord {
    GffRecord {
        seqname: CHR.into(),
        feature_type,
        start,
        stop,
        strand,
        gene_id: GeneId::Ensembl(gene.into()),
        gene_name: GeneSymbol::Symbol(gene.into()),
        gene_type: GeneType::CodingGene,
        transcript_id: TranscriptId::Ensembl(tx.into()),
    }
}

/// A non-coding gene's EXON, which is what the ncRNA track is built from —
/// the `gene` row's span would carry the introns.
fn non_coding_rec(gene: &str, start: i64, stop: i64, strand: Strand) -> GffRecord {
    GffRecord {
        gene_type: GeneType::LincRNA,
        transcript_id: TranscriptId::Missing,
        ..rec("none", gene, FeatureType::Exon, start, stop, strand)
    }
}

/// A two-exon coding transcript: exons 1000..1199 and 1500..1699, CDS
/// 1100..1599 with a separate stop codon, leaving a 100 nt 5'UTR and a 97 nt
/// 3'UTR on the forward strand.
fn two_exon_tx(tx: &str, gene: &str, strand: Strand) -> Vec<GffRecord> {
    let stop_codon = match strand {
        Strand::Forward => (1600, 1602),
        Strand::Backward => (1097, 1099),
    };
    vec![
        rec(tx, gene, FeatureType::Exon, 1000, 1199, strand),
        rec(tx, gene, FeatureType::Exon, 1500, 1699, strand),
        rec(tx, gene, FeatureType::CDS, 1100, 1199, strand),
        rec(tx, gene, FeatureType::CDS, 1500, 1599, strand),
        rec(
            tx,
            gene,
            FeatureType::StopCodon,
            stop_codon.0,
            stop_codon.1,
            strand,
        ),
    ]
}

/// Sites carry 0-based positions; every coordinate in these tests is the
/// 1-based GFF one. The strand matters: placement is same-strand only, as
/// MetaPlotR's `intersectBed -s` is.
fn site(pos_1based: i64, strand: Strand) -> GenomicSite {
    GenomicSite {
        chr: CHR.into(),
        position: pos_1based - 1,
        strand,
    }
}

/// Build the index the way `run_metagene` does, then place sites on it.
fn place_on(records: &[GffRecord], positions: &[i64], strand: Strand) -> Vec<SiteAssignment> {
    let models = build_transcript_models(records);
    let nc = non_coding_bodies(records);
    let index = RegionIndex::build(&models, &nc);
    let sites: Vec<GenomicSite> = positions.iter().map(|&p| site(p, strand)).collect();
    assign_sites(&sites, &index, IsoformPolicy::All).0
}

fn place(records: &[GffRecord], positions: &[i64]) -> Vec<SiteAssignment> {
    place_on(records, positions, Strand::Forward)
}

/// `[5'UTR, CDS, 3'UTR, ncRNA]` totals for one site on `strand`.
fn classify_on(records: &[GffRecord], pos_1based: i64, strand: Strand) -> [usize; 4] {
    let mut out = [0usize; 4];
    for a in place_on(records, &[pos_1based], strand) {
        out[a.region] += 1;
    }
    out
}

fn classify(records: &[GffRecord], pos_1based: i64) -> [usize; 4] {
    classify_on(records, pos_1based, Strand::Forward)
}

/////////////////////////
// The interval index   //
/////////////////////////

fn feature(strand: Strand, intervals: &[(i64, i64)]) -> TranscriptModel {
    let total: i64 = intervals.iter().map(|&(s, e)| e - s + 1).sum();
    TranscriptModel {
        gene_id: GeneId::Ensembl("G".into()),
        gene_name: GeneSymbol::Symbol("G".into()),
        transcript_id: TranscriptId::Ensembl("T".into()),
        seqname: CHR.into(),
        strand,
        utr5: vec![],
        cds: intervals.to_vec(),
        utr3: vec![],
        utr5_size: 0,
        cds_size: total,
        utr3_size: 0,
    }
}

fn hits_on(idx: &RegionIndex, pos: i64, strand: Strand) -> Vec<IndexedInterval> {
    let mut v = Vec::new();
    idx.find_all(CHR, pos, strand, &mut v);
    v
}

fn hits(idx: &RegionIndex, pos: i64) -> Vec<IndexedInterval> {
    hits_on(idx, pos, Strand::Forward)
}

#[test]
fn find_reaches_past_a_shorter_later_interval() {
    // Sorted by start: [(100, 5000), (200, 300)]. Walking back from the
    // rightmost start <= 1000 hits (200, 300) first; `start < position` is no
    // reason to stop, because (100, 5000) still contains 1000.
    let idx = RegionIndex::build(
        &[
            feature(Strand::Forward, &[(100, 5000)]),
            feature(Strand::Forward, &[(200, 300)]),
        ],
        &[],
    );

    let h = hits(&idx, 1000);
    assert_eq!(h.len(), 1);
    assert_eq!((h[0].start, h[0].stop), (100, 5000));

    assert_eq!(hits(&idx, 5001).len(), 0);
    assert_eq!(hits(&idx, 99).len(), 0);
    let mut v = Vec::new();
    idx.find_all("chrX", 1000, Strand::Forward, &mut v);
    assert!(v.is_empty());
}

#[test]
fn find_scans_back_over_several_non_containing_intervals() {
    let idx = RegionIndex::build(
        &[
            feature(Strand::Forward, &[(100, 5000)]),
            feature(Strand::Forward, &[(200, 300)]),
            feature(Strand::Forward, &[(400, 500)]),
            feature(Strand::Forward, &[(600, 700)]),
        ],
        &[],
    );
    let h = hits(&idx, 4000);
    assert_eq!(h.len(), 1);
    assert_eq!((h[0].start, h[0].stop), (100, 5000));

    // …and the max-stop bound still lets it give up: nothing covers 5500.
    assert_eq!(hits(&idx, 5500).len(), 0);
}

#[test]
fn overlapping_transcripts_each_take_the_site() {
    // MetaPlotR emits one row per site and transcript. A site inside two
    // transcripts must be placed on both, not just the first found.
    let idx = RegionIndex::build(
        &[
            feature(Strand::Forward, &[(100, 5000)]),
            feature(Strand::Forward, &[(900, 1100)]),
        ],
        &[],
    );
    assert_eq!(hits(&idx, 1000).len(), 2);
}

///////////////////////////////
// Spliced relative position  //
///////////////////////////////

#[test]
fn relative_position_runs_along_the_spliced_feature() {
    // Two exons, 100 bp each: the intron 200..299 must consume no coordinate.
    let idx = RegionIndex::build(&[feature(Strand::Forward, &[(100, 199), (300, 399)])], &[]);
    let rel = |pos| hits(&idx, pos)[0].relative_pos(pos);

    assert_eq!(rel(100), 0);
    assert_eq!(rel(199), 99);
    // First base of exon 2 is base 100 of the transcript, not base 200.
    assert_eq!(rel(300), 100);
    assert_eq!(rel(399), 199);
}

#[test]
fn reverse_strand_relative_position_runs_five_to_three() {
    let idx = RegionIndex::build(&[feature(Strand::Backward, &[(100, 199), (300, 399)])], &[]);
    let rel = |pos| hits_on(&idx, pos, Strand::Backward)[0].relative_pos(pos);

    assert_eq!(rel(399), 0);
    assert_eq!(rel(300), 99);
    assert_eq!(rel(199), 100);
    assert_eq!(rel(100), 199);

    let b = |p: i64| {
        hits_on(&idx, p, Strand::Backward)[0]
            .place(0, p, 1.0)
            .bin(10)
    };
    assert_eq!(b(399), 0);
    assert_eq!(b(100), 9);
}

#[test]
fn bins_are_a_fraction_of_the_spliced_length() {
    let idx = RegionIndex::build(&[feature(Strand::Forward, &[(100, 199), (300, 399)])], &[]);
    let bin = |pos| hits(&idx, pos)[0].place(0, pos, 1.0).bin(10);

    assert_eq!(bin(100), 0);
    // Base 100 of 200 is the middle of the transcript, though it is 3/4 of the
    // way across the genomic span.
    assert_eq!(bin(300), 5);
    assert_eq!(bin(399), 9);
}

//////////////////////////
// Region assignment     //
//////////////////////////

#[test]
fn regions_are_disjoint_so_no_priority_order_is_needed() {
    // Under the gene-union model a site was tested 5'UTR -> CDS -> 3'UTR and an
    // oversized CDS could claim a 3'UTR site. Within one transcript the three
    // regions cannot overlap, so each site has exactly one home.
    for strand in [Strand::Forward, Strand::Backward] {
        let records = two_exon_tx("ENST1", "G", strand);
        for pos in [1050, 1150, 1550, 1650] {
            let total: usize = classify_on(&records, pos, strand).iter().sum();
            assert!(total <= 1, "position {pos} landed twice on {strand:?}");
        }
    }
}

#[test]
fn intronic_site_between_cds_exons_is_not_called_cds() {
    for strand in [Strand::Forward, Strand::Backward] {
        let records = two_exon_tx("ENST1", "G", strand);
        // 1300 sits in the intron between exons 1000..1199 and 1500..1699.
        assert_eq!(
            classify_on(&records, 1300, strand),
            [0, 0, 0, 0],
            "{strand:?}"
        );
        // The coding exons themselves are still CDS.
        assert_eq!(
            classify_on(&records, 1150, strand),
            [0, 1, 0, 0],
            "{strand:?}"
        );
        assert_eq!(
            classify_on(&records, 1550, strand),
            [0, 1, 0, 0],
            "{strand:?}"
        );
    }
}

#[test]
fn stop_codon_bins_as_the_last_cds_bin_not_the_first_utr3_bin() {
    // The whole reason the coding extent is widened over the stop codon.
    let records = two_exon_tx("ENST1", "G", Strand::Forward);
    assert_eq!(classify(&records, 1601), [0, 1, 0, 0]);

    let a = place(&records, &[1601]);
    assert_eq!(a.len(), 1);
    assert_eq!(a[0].region, CDS);
    assert_eq!(
        a[0].bin(10),
        9,
        "the stop codon is at the 3' end of the CDS"
    );
}

#[test]
fn rel_location_spans_zero_to_three() {
    let records = two_exon_tx("ENST1", "G", Strand::Forward);
    let at = |p: i64| place(&records, &[p])[0].rel_location();

    // First 5'UTR base is 0.0, first CDS base 1.0, first 3'UTR base 2.0.
    assert!((at(1000) - 0.0).abs() < 1e-9);
    assert!((at(1100) - 1.0).abs() < 1e-9);
    assert!((at(1603) - 2.0).abs() < 1e-9);
    // …and the last 3'UTR base stays below 3.
    assert!(at(1699) < 3.0 && at(1699) > 2.9);
}

#[test]
fn reverse_strand_rel_location_also_runs_five_to_three() {
    let records = two_exon_tx("ENST1", "G", Strand::Backward);
    let at = |p: i64| place_on(&records, &[p], Strand::Backward)[0].rel_location();
    // On the reverse strand the 5'UTR is the HIGH end, and every region reads
    // 5'->3' as the genomic coordinate DECREASES. So the transcript starts at
    // 1699 and the 3'UTR's first base is its highest coordinate, 1096 — the
    // low end, 1000, is where the transcript ENDS.
    assert!((at(1699) - 0.0).abs() < 1e-9);
    assert!((at(1096) - 2.0).abs() < 1e-9);
    assert!(at(1000) > 2.9 && at(1000) < 3.0);
}

//////////////////////////
// Isoform placement     //
//////////////////////////

#[test]
fn single_isoform_genes_are_kept() {
    // MetaPlotR's published `dist[duplicated(gene_name), ]` drops these.
    let records = two_exon_tx("ENST1", "G", Strand::Forward);
    assert_eq!(classify(&records, 1150), [0, 1, 0, 0]);
}

#[test]
fn a_site_is_placed_on_every_isoform_that_contains_it() {
    // ENST2 (1501 nt) is the longer and covers 1750; so does ENST1's third
    // exon. Electing the longest would place it once; both carry it.
    let mut records = two_exon_tx("ENST1", "G", Strand::Forward);
    records.push(rec(
        "ENST1",
        "G",
        FeatureType::Exon,
        1700,
        1800,
        Strand::Forward,
    ));
    records.extend(vec![
        rec("ENST2", "G", FeatureType::Exon, 1000, 2500, Strand::Forward),
        rec("ENST2", "G", FeatureType::CDS, 1100, 1599, Strand::Forward),
        rec(
            "ENST2",
            "G",
            FeatureType::StopCodon,
            1600,
            1602,
            Strand::Forward,
        ),
    ]);

    let a = place(&records, &[1750]);
    assert_eq!(a.len(), 2);
    assert_ne!(a[0].model, a[1].model);
    assert!(a.iter().all(|a| a.weight == 1.0), "unsplit under `all`");
}

//////////////////////////
// Scale factors        //
//////////////////////////

#[test]
fn scale_factors_count_each_transcript_once() {
    // Two transcripts with very different 3'UTRs. The one with the LARGE
    // 3'UTR carries nine sites, the other one. Over the site table, as
    // MetaPlotR's script takes it, the many-site transcript would set the
    // median alone (9000); each transcript counts once instead.
    let models = vec![
        sized_model(100, 1000, 100),
        sized_model(100, 1000, 9000),
        sized_model(100, 1000, 50), // carries no site: not in the medians
    ];
    let mk = |model: u32| SiteAssignment {
        site: 0,
        model: Some(model),
        region: UTR3,
        rel: 0,
        total_len: 1,
        weight: 1.0,
    };
    let mut assignments: Vec<SiteAssignment> = (0..9).map(|_| mk(1)).collect();
    assignments.push(mk(0));

    let sf = scale_factors(&assignments, &models).expect("scale factors");
    assert_eq!(
        sf.median()[UTR3],
        4550.0,
        "a site-weighted median would be 9000"
    );
    assert!((sf.utr3_sf - 4.55).abs() < 1e-9);
}

#[test]
fn scale_factors_need_a_coding_assignment() {
    let models = vec![sized_model(100, 0, 100)];
    let assignments = vec![SiteAssignment {
        site: 0,
        model: Some(0),
        region: UTR3,
        rel: 0,
        total_len: 1,
        weight: 1.0,
    }];
    assert!(scale_factors(&assignments, &models).is_none());
}

//////////////////////////
// Bin allocation       //
//////////////////////////

#[test]
fn bins_follow_the_medians_and_always_sum_to_n() {
    // A CDS twice either UTR: half the bins, a quarter each.
    assert_eq!(allocate_bins(200, &[100, 200, 100]), [50, 100, 50]);
    // Non-dividing totals still sum exactly.
    for n in [7usize, 57, 199, 200, 1000] {
        for m in [[155, 1026, 1720], [1, 1, 1], [3, 100, 7]] {
            let got = allocate_bins(n, &m);
            assert_eq!(got.iter().sum::<usize>(), n, "n={n} m={m:?} got={got:?}");
        }
    }
}

#[test]
fn a_region_with_no_length_gets_no_bins() {
    assert_eq!(allocate_bins(100, &[0, 100, 100]), [0, 50, 50]);
}

//////////////////////////
// Output                //
//////////////////////////

/// A model carrying nothing but its three region sizes, for the median tests.
fn sized_model(utr5: i64, cds: i64, utr3: i64) -> TranscriptModel {
    TranscriptModel {
        utr5_size: utr5,
        cds_size: cds,
        utr3_size: utr3,
        ..feature(Strand::Forward, &[(1, 1)])
    }
}

fn histogram(counts: [Vec<f64>; 4], sf5: f64, sf3: f64) -> GeneFeatureHistogram {
    GeneFeatureHistogram {
        counts,
        scale: ScaleFactors {
            twice_median: [310, 2052, 3440],
            utr5_sf: sf5,
            utr3_sf: sf3,
        },
    }
}

#[test]
fn bin_edges_reproduce_metaplotr_rescaling() {
    let h = histogram(
        [vec![0.0; 2], vec![0.0; 2], vec![0.0; 2], vec![]],
        0.15,
        1.68,
    );
    // The CDS keeps width 1 and spans 1..2; each UTR is scaled to its median
    // size relative to the CDS and butts against it.
    assert!((h.bin_edges(UTR5, 0).0 - 0.85).abs() < 1e-9);
    assert!((h.bin_edges(UTR5, 1).1 - 1.0).abs() < 1e-9);
    assert!((h.bin_edges(CDS, 0).0 - 1.0).abs() < 1e-9);
    assert!((h.bin_edges(CDS, 1).1 - 2.0).abs() < 1e-9);
    assert!((h.bin_edges(UTR3, 0).0 - 2.0).abs() < 1e-9);
    assert!((h.bin_edges(UTR3, 1).1 - 3.68).abs() < 1e-9);
}

/// Verbatim from `rel_and_abs_dist_calc.pl:38-40`, which prints its header in
/// three statements. If this drifts, `visualize_metagenes.R` stops reading our
/// output, which is the only reason the file exists.
#[test]
fn dist_measures_header_matches_metaplotr() {
    let metaplotr = "chr\tcoord\tgene_name\trefseqID\trel_location\t\
                     utr5_st\tutr5_end\tcds_st\tcds_end\tutr3_st\tutr3_end\t\
                     utr5_size\tcds_size\tutr3_size";
    assert!(
        DIST_MEASURES_HEADER.starts_with(metaplotr),
        "MetaPlotR's 14 columns must come first, in its order\n  ours: {DIST_MEASURES_HEADER}"
    );
    // Our extras go after theirs, so a positional reader is unaffected.
    assert_eq!(
        &DIST_MEASURES_HEADER[metaplotr.len()..],
        "\tstrand\trescaled_location\tweight"
    );
}

#[test]
fn utr3_st_is_the_signed_distance_from_the_stop_codon() {
    // The column MetaPlotR's feature-distance plot is drawn on. Positions are
    // 1-based along the mature transcript, so the first 3'UTR base sits exactly
    // on the boundary and reads 0; a CDS base reads negative.
    let records = two_exon_tx("ENST1", "G", Strand::Forward);
    let models = build_transcript_models(&records);
    let m = &models[0];

    let utr3_st = |a: &SiteAssignment| {
        let preceding = match a.region {
            UTR5 => 0,
            CDS => m.utr5_size,
            _ => m.utr5_size + m.cds_size,
        };
        (preceding + a.rel + 1) - (m.utr5_size + m.cds_size + 1)
    };

    // 1603 is the first 3'UTR base (the stop codon 1600..1602 folds into CDS).
    let first_utr3 = &place(&records, &[1603])[0];
    assert_eq!(first_utr3.region, UTR3);
    assert_eq!(utr3_st(first_utr3), 0);

    // The last CDS base is one before it.
    let last_cds = &place(&records, &[1602])[0];
    assert_eq!(last_cds.region, CDS);
    assert_eq!(utr3_st(last_cds), -1);

    // …and a base deeper into the 3'UTR is positive.
    let later = &place(&records, &[1613])[0];
    assert_eq!(utr3_st(later), 10);
}

#[test]
fn tsv_keeps_its_first_three_columns_and_integrates_to_one() {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("mg.tsv");
    let path = path.to_str().unwrap();

    let h = histogram(
        [vec![1.0, 3.0], vec![10.0, 20.5], vec![5.0, 0.5], vec![]],
        0.15,
        1.68,
    );
    h.to_tsv(path).expect("write");

    let text = std::fs::read_to_string(path).expect("read");
    let mut lines = text.lines();
    assert_eq!(
        lines.next().unwrap(),
        "#feature\tgenomic_bin\tcount\tbin_start\tbin_end\tfrac\tdensity",
        "the first three columns are an output contract"
    );

    // density integrates to 1 over the coding regions; whole counts are
    // written as integers, split ones as decimals.
    let mut integral = 0.0;
    let mut written = Vec::new();
    for line in lines {
        let f: Vec<&str> = line.split('\t').collect();
        written.push(f[2].to_string());
        let (lo, hi): (f64, f64) = (f[3].parse().unwrap(), f[4].parse().unwrap());
        integral += f[6].parse::<f64>().unwrap() * (hi - lo);
    }
    assert!(
        (integral - 1.0).abs() < 1e-4,
        "density integrated to {integral}"
    );
    assert_eq!(written, ["1", "3", "10", "20.5000", "5", "0.5000"]);
}

//////////////////////////
// Non-coding track     //
//////////////////////////

#[test]
fn the_grid_sizes_all_four_tracks_together() {
    // All four widths come from one constructor. They were split once before,
    // with the ncRNA width set by its own statement, and this is the unit that
    // can observe them agreeing.
    let scale = ScaleFactors {
        twice_median: [200, 400, 200],
        utr5_sf: 0.5,
        utr3_sf: 0.5,
    };

    let off = BinGrid::new(100, Some(&scale), false).0;
    assert_eq!(off, [25, 50, 25, 0], "no ncRNA track unless asked for");
    assert_eq!(off[..3].iter().sum::<usize>(), 100);

    let on = BinGrid::new(100, Some(&scale), true).0;
    assert_eq!(
        on,
        [25, 50, 25, 100],
        "ncRNA gets the whole budget, own axis"
    );
    // Asking for the track and getting it is the ONLY difference: the three
    // coding widths must not move when the flag flips.
    assert_eq!(off[..3], on[..3]);
}

#[test]
fn a_non_coding_gene_is_not_a_coding_transcript() {
    let records = vec![non_coding_rec("NC", 5000, 6000, Strand::Forward)];
    assert!(build_transcript_models(&records).is_empty());
}

#[test]
fn non_coding_genes_keep_whole_gene_boundaries() {
    let records = vec![non_coding_rec("NC", 5000, 6000, Strand::Forward)];
    assert_eq!(classify(&records, 5500), [0, 0, 0, 1]);
    assert_eq!(classify(&records, 4999), [0, 0, 0, 0]);
}

//////////////////////////////
// Regressions from review   //
//////////////////////////////

#[test]
fn a_site_is_not_placed_on_an_antisense_transcript() {
    // MetaPlotR intersects with `-s`. Without that filter a + strand site also
    // lands on every - strand transcript overlapping it, and `relative_pos`
    // mirrors, so the phantom copy sits at 1-p instead of p.
    let mut records = two_exon_tx("FWD", "GF", Strand::Forward);
    records.extend(two_exon_tx("REV", "GR", Strand::Backward));

    // 1150 is a CDS base of both fixtures; each strand may claim it once.
    assert_eq!(classify_on(&records, 1150, Strand::Forward), [0, 1, 0, 0]);
    assert_eq!(classify_on(&records, 1150, Strand::Backward), [0, 1, 0, 0]);
}

#[test]
fn a_coding_site_is_not_also_counted_on_the_ncrna_track() {
    // Non-coding bodies span coding genes constantly. The ncRNA track is a
    // fallback, not a parallel one, or the coding profile is reprinted on a
    // track that has no stop codon.
    let mut records = two_exon_tx("ENST1", "G", Strand::Forward);
    records.push(non_coding_rec("NC", 900, 2300, Strand::Forward));

    // A coding base goes to CDS only, even though the ncRNA body covers it.
    assert_eq!(classify(&records, 1150), [0, 1, 0, 0]);
    // A base no coding region claims still reaches the ncRNA track.
    assert_eq!(classify(&records, 1300), [0, 0, 0, 1]);
}

#[test]
fn the_ncrna_track_is_spliced_not_the_gene_span() {
    // Two exons with an intron between them. A site in the intron has no
    // transcript position, exactly as inside a coding gene.
    let records = vec![
        non_coding_rec("NC", 1000, 1099, Strand::Forward),
        non_coding_rec("NC", 1500, 1599, Strand::Forward),
    ];
    assert_eq!(classify(&records, 1050), [0, 0, 0, 1]);
    assert_eq!(classify(&records, 1550), [0, 0, 0, 1]);
    assert_eq!(classify(&records, 1300), [0, 0, 0, 0], "intronic");

    // …and the coordinate skips the intron: exon 2's first base is base 100
    // of a 200 nt body, so it bins at the midpoint.
    let a = place(&records, &[1500]);
    assert_eq!(a[0].bin(10), 5);
}

#[test]
fn a_region_with_sites_never_loses_all_its_bins() {
    // At the measured medians, `--bins 10` used to give the 5'UTR zero bins,
    // and `accumulate` then dropped every 5'UTR site while `to_tsv` took its
    // denominator from the binned counts — so the file still integrated to 1
    // with a whole track missing.
    for n in [1usize, 3, 7, 10, 57, 200] {
        let got = allocate_bins(n, &[310, 2052, 3440]);
        assert_eq!(got.iter().sum::<usize>(), n, "n={n} got={got:?}");
        if n >= 3 {
            assert!(
                got.iter().all(|&b| b > 0),
                "n={n} left a represented region with no bins: {got:?}"
            );
        }
    }
}

#[test]
fn twice_median_returns_both_middle_values_on_even_n() {
    // The even-n branch is `*mid + max(lower)`. The tempting `2 * *mid` agrees
    // whenever the two middle order statistics are equal, which every other
    // test in this file happens to satisfy — so this one uses a sample where
    // they differ and the two formulas disagree.
    let mut v = vec![1i64, 2, 3, 4]; // middles 2 and 3 -> 5
    assert_eq!(twice_median(&mut v), 5);

    let mut v = vec![10i64, 1, 7, 3]; // sorted 1,3,7,10 -> 3+7 = 10
    assert_eq!(twice_median(&mut v), 10);

    // Odd n is just twice the middle, and the empty case is 0.
    let mut v = vec![5i64, 1, 9];
    assert_eq!(twice_median(&mut v), 10);
    assert_eq!(twice_median(&mut []), 0);
}

#[test]
fn dist_measures_writes_metaplotr_columns_and_distances() {
    // The header constant is checked elsewhere; this exercises the WRITER, so
    // a lost tab or a reordered field is caught rather than assumed.
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("dist.tsv");
    let path = path.to_str().unwrap();

    let records = two_exon_tx("ENST1", "G", Strand::Forward);
    let models = build_transcript_models(&records);
    let sites = vec![site(1603, Strand::Forward)]; // first 3'UTR base
    let index = RegionIndex::build(&models, &[]);
    let (assignments, _) = assign_sites(&sites, &index, IsoformPolicy::All);
    let scale = ScaleFactors {
        twice_median: [200, 400, 200],
        utr5_sf: 0.5,
        utr3_sf: 0.5,
    };
    write_dist_measures(path, &sites, &assignments, &models, &scale).expect("write");

    let text = std::fs::read_to_string(path).expect("read");
    let mut lines = text.lines();
    assert_eq!(lines.next().unwrap(), DIST_MEASURES_HEADER);

    let row: Vec<&str> = lines.next().expect("one row").split('\t').collect();
    assert_eq!(
        row.len(),
        17,
        "14 MetaPlotR columns plus our three: {row:?}"
    );
    assert_eq!(row[0], "chr1");
    assert_eq!(row[1], "1603", "coord is 1-based, as MetaPlotR's is");
    assert_eq!(row[3], "ENST1");
    // utr3_st is column 9: the first 3'UTR base sits exactly on the boundary.
    assert_eq!(row[9], "0");
    // …and cds_end (column 8) is one past the last CDS base.
    assert_eq!(row[8], "1");
    // The three sizes follow: 5'UTR 100, CDS 203 (stop codon folded in), 3'UTR 97.
    assert_eq!((row[11], row[12], row[13]), ("100", "203", "97"));
    assert_eq!(row[16], "1.000000", "one isoform: the whole site");
    assert!(lines.next().is_none(), "exactly one assignment expected");
}

///////////////////////////
// Fixed layout, subsets //
///////////////////////////

#[test]
fn a_subset_is_counted_on_the_full_set_axis() {
    let records = two_exon_tx("T1", "GENE1", Strand::Forward);
    let models = MetaModels::from_records(&records);
    // 5'UTR, CDS (two), 3'UTR, and one off any transcript.
    let positions = [1050, 1150, 1550, 1650, 5000];
    let sites: Vec<GenomicSite> = positions
        .iter()
        .map(|&p| site(p, Strand::Forward))
        .collect();
    let layout = models
        .layout(&sites, 30)
        .expect("sites on a coding transcript");
    assert_eq!(layout.unassigned, 1);
    assert_eq!(layout.region_bins().iter().sum::<usize>(), 30);

    let all = layout.counts(|_| 1.0);
    assert_eq!(all.iter().sum::<f64>(), 4.0);
    // Keeping only the CDS sites leaves the axis as it was.
    let cds_only = layout.counts(|i| f64::from(u8::from(i == 1 || i == 2)));
    assert_eq!(cds_only.len(), all.len());
    assert_eq!(cds_only.iter().sum::<f64>(), 2.0);
    let [utr5, cds, _] = layout.region_bins();
    assert_eq!(cds_only[utr5..utr5 + cds].iter().sum::<f64>(), 2.0);
    assert!(cds_only.iter().zip(&all).all(|(k, a)| k <= a));
    assert!(layout.counts(|_| 0.0).iter().all(|&n| n == 0.0));
    // A weight sums: every site at 3 is three times the site count.
    let tripled = layout.counts(|_| 3.0);
    assert!(tripled.iter().zip(&all).all(|(t, a)| *t == 3.0 * a));
}

/// Two isoforms of one gene: ENST1 is [`two_exon_tx`]; ENST2 shares its
/// first exon and CDS start but ends in a short last exon of its own, so it
/// would lose a longest-isoform election.
fn two_isoforms() -> Vec<GffRecord> {
    let mut records = two_exon_tx("ENST1", "G", Strand::Forward);
    records.extend([
        rec("ENST2", "G", FeatureType::Exon, 1000, 1199, Strand::Forward),
        rec("ENST2", "G", FeatureType::Exon, 1300, 1349, Strand::Forward),
        rec("ENST2", "G", FeatureType::CDS, 1100, 1199, Strand::Forward),
        rec("ENST2", "G", FeatureType::CDS, 1300, 1320, Strand::Forward),
    ]);
    records
}

#[test]
fn the_layout_splits_a_site_over_the_isoforms_it_lies_in() {
    let records = two_isoforms();
    let models = MetaModels::from_records(&records);
    // 1150: CDS of both isoforms. 1330: only ENST2's last exon (3'UTR),
    // which a longest-isoform election would leave off-transcript.
    let sites = [site(1150, Strand::Forward), site(1330, Strand::Forward)];
    let layout = models.layout(&sites, 30).expect("placed");
    assert_eq!(layout.unassigned, 0);
    let weights: Vec<(u32, f64)> = layout
        .assignments
        .iter()
        .map(|a| (a.site, a.weight))
        .collect();
    assert_eq!(weights.len(), 3, "{weights:?}");
    for (site, weight) in weights {
        let want = if site == 0 { 0.5 } else { 1.0 };
        assert_eq!(weight, want, "site {site}");
    }
    // Each site adds 1 in all.
    assert_eq!(layout.counts(|_| 1.0).iter().sum::<f64>(), 2.0);
}

#[test]
fn weighted_isoforms_split_the_site_in_the_distance_table() {
    let dir = tempfile::tempdir().expect("tempdir");
    let path = dir.path().join("dist.tsv");
    let path = path.to_str().unwrap();
    let models = build_transcript_models(&two_isoforms());
    let index = RegionIndex::build(&models, &[]);
    let sites = vec![site(1150, Strand::Forward)];
    let (assignments, _) = assign_sites(&sites, &index, IsoformPolicy::Weighted);
    assert_eq!(assignments.len(), 2);
    let scale = scale_factors(&assignments, &models).expect("coding");
    write_dist_measures(path, &sites, &assignments, &models, &scale).expect("write");
    let text = std::fs::read_to_string(path).expect("read");
    let weights: Vec<&str> = text
        .lines()
        .skip(1)
        .map(|l| l.rsplit('\t').next().unwrap())
        .collect();
    assert_eq!(weights, ["0.500000", "0.500000"]);
}

#[test]
fn no_layout_without_a_coding_placement() {
    let records = two_exon_tx("T1", "GENE1", Strand::Forward);
    let models = MetaModels::from_records(&records);
    // Wrong strand: placement is same-strand only.
    let sites = [site(1150, Strand::Backward)];
    assert!(models.layout(&sites, 30).is_none());
}

#[test]
fn region_colours_reach_the_terminal() {
    use ratatui::style::Color;
    assert_eq!(region_style(UTR5).fg, Some(Color::Rgb(0x56, 0xb4, 0xe9)));
    assert_eq!(region_style(CDS).fg, Some(Color::Rgb(0x00, 0x72, 0xb2)));
    assert_eq!(region_style(UTR3).fg, Some(Color::Rgb(0x00, 0x9e, 0x73)));
}

#[test]
fn a_transcript_without_a_utr_is_left_out_of_that_utr_median() {
    // Two CDS fragments with no 5'UTR outnumber the one transcript with
    // one; a zero in the median would leave the 5'UTR no bins.
    let models = vec![
        sized_model(0, 1000, 400),
        sized_model(0, 1000, 600),
        sized_model(150, 1000, 0),
    ];
    let assignments: Vec<SiteAssignment> = (0..3)
        .map(|m| SiteAssignment {
            site: 0,
            model: Some(m),
            region: CDS,
            rel: 0,
            total_len: 1,
            weight: 1.0,
        })
        .collect();
    let sf = scale_factors(&assignments, &models).expect("coding");
    assert_eq!(sf.median(), [150.0, 1000.0, 500.0]);
    assert!(allocate_bins(48, &sf.twice_median)[UTR5] > 0);
}

#[test]
fn split_shares_that_sum_to_a_whole_print_as_one() {
    let n: f64 = [1.0 / 3.0; 3].iter().sum::<f64>() + [0.2; 5].iter().sum::<f64>();
    assert_eq!(count_text(n, 4), "2");
    assert_eq!(count_text(0.5, 4), "0.5000");
    assert_eq!(count_text(3.0, 1), "3");
}
