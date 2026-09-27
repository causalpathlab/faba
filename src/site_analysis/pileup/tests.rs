use super::*;

#[test]
fn parse_site_and_mixture_rows() {
    // site output: gene/modality/chr:pos
    assert_eq!(
        parse_row_name_full("ENSG00000139618_BRCA2/m6A/chr13:32350000"),
        Some(("ENSG00000139618_BRCA2", "m6A", "chr13", 32350000))
    );
    // mixture output: gene/modality/component
    assert_eq!(
        parse_row_name_full("ENSG00000060558_GNA15/m6A/0"),
        Some(("ENSG00000060558_GNA15", "m6A", "", 0))
    );
    // count rows (detail neither chr:pos nor an integer) don't parse
    assert_eq!(parse_row_name_full("gene_0/count/spliced"), None);
}

#[test]
fn relaxed_gene_matching() {
    let gp = "ENSG00000060558_GNA15";
    // symbol, any case
    assert!(gene_matches("GNA15", &query_symbol("GNA15"), gp));
    assert!(gene_matches("gna15", &query_symbol("gna15"), gp));
    // Ensembl ID
    assert!(gene_matches(
        "ENSG00000060558",
        &query_symbol("ENSG00000060558"),
        gp
    ));
    // full composite
    assert!(gene_matches(gp, &query_symbol(gp), gp));
    // partial substring must NOT match (consistent with aux-data scheme)
    assert!(!gene_matches("GNA", &query_symbol("GNA"), gp));
    assert!(!gene_matches(
        "RPL",
        &query_symbol("RPL"),
        "ENSG00000063177_RPL18"
    ));
}

#[test]
fn region_parsing_and_matching() {
    let r = parse_region("chr17:1000-2000").unwrap();
    assert_eq!((r.chr.as_ref(), r.lb, r.ub), ("chr17", 1000, 2000));
    // reversed bounds are normalized
    let r = parse_region("17:2000-1000").unwrap();
    assert_eq!((r.lb, r.ub), (1000, 2000));
    // malformed specs error
    assert!(parse_region("chr17").is_err());
    assert!(parse_region("chr17:1000").is_err());
    assert!(parse_region(":1-2").is_err());

    // chr-prefix tolerant (shared genomic_data::chr_eq)
    assert!(chr_eq("chr17", "17"));
    assert!(chr_eq("17", "17"));
    assert!(!chr_eq("chr1", "chr2"));
}

#[test]
fn selector_gene_or_region_union() {
    let sel = Selector::build(&["GNA15".into()], &["chr17:100-200".into()]).unwrap();
    // gene branch (mixture row: no chr, component ordinal)
    assert!(sel.selects("ENSG00000060558_GNA15", "", 0));
    // region branch (different gene, but inside the window)
    assert!(sel.selects("ENSG1_OTHER", "chr17", 150));
    // outside both
    assert!(!sel.selects("ENSG1_OTHER", "chr17", 999));
    assert!(!sel.selects("ENSG1_OTHER", "chr9", 150));

    // at least one selector is required
    assert!(Selector::build(&[], &[]).is_err());
    // region-only is fine
    assert!(Selector::build(&[], &["chr1:1-9".into()]).is_ok());
}

#[test]
fn thousands_and_axis_mapping() {
    assert_eq!(fmt_thousands(26781984), "26,781,984");
    assert_eq!(fmt_thousands(767), "767");
    assert_eq!(fmt_thousands(0), "0");
    assert_eq!(fmt_thousands(-1234), "-1,234");

    // first site -> first column, last site -> last column
    assert_eq!(pos_to_col(100, 100, 200, 10), 0);
    assert_eq!(pos_to_col(200, 100, 200, 10), 9);
    assert_eq!(pos_to_col(150, 100, 200, 10), 5);
    // degenerate single-position extent collapses to column 0
    assert_eq!(pos_to_col(100, 100, 100, 10), 0);
}

#[test]
fn distinct_positions_dedups_sorted() {
    let p = [(10, 1.0), (10, 2.0), (20, 0.5), (30, 0.0)];
    assert_eq!(distinct_positions(&p), vec![10, 20, 30]);
}

#[test]
fn aggregate_labels() {
    let mut one: FxHashMap<Box<str>, usize> = FxHashMap::default();
    one.insert("ENSG1_GNA15".into(), 2);
    assert_eq!(summarize_genes(&one).as_ref(), "ENSG1_GNA15");

    let mut many: FxHashMap<Box<str>, usize> = FxHashMap::default();
    many.insert("ENSG1_A".into(), 1);
    many.insert("ENSG2_B".into(), 3);
    let label = summarize_genes(&many);
    assert!(label.starts_with("2 genes: "), "got {label}");

    // chr label: empty -> component, single -> that chr, mixed -> *
    let mix: Vec<Box<str>> = vec!["".into(), "".into()];
    assert_eq!(summarize_chr(&mix).as_ref(), "component");
    let single: Vec<Box<str>> = vec!["chr13".into()];
    assert_eq!(summarize_chr(&single).as_ref(), "chr13");
    let multi: Vec<Box<str>> = vec!["chr1".into(), "chr2".into()];
    assert_eq!(summarize_chr(&multi).as_ref(), "*");
}

#[test]
fn channel_rows_pile_up_the_converted_channel() {
    assert_eq!(
        parse_row_name_full("ENSG1_GENE1/m6a/chr1:100/methylated"),
        Some(("ENSG1_GENE1", "m6a", "chr1", 100))
    );
    assert_eq!(
        parse_row_name_full("ENSG1_GENE1/atoi/chr1:200/edited"),
        Some(("ENSG1_GENE1", "atoi", "chr1", 200))
    );
    assert_eq!(
        parse_row_name_full("ENSG1_GENE1/m6a/chr1:100/unmethylated"),
        None
    );
    assert_eq!(
        parse_row_name_full("ENSG1_GENE1/atoi/chr1:200/unedited"),
        None
    );
}

#[test]
fn wildcard_patterns() {
    assert!(wildcard("*_wt_*", "out/rep1_wt_m6a_site.zarr.zip"));
    assert!(!wildcard("*_wt_*", "out/rep1_mut_m6a_site.zarr.zip"));
    assert!(wildcard("rep?_*", "rep2_x"));
    assert!(wildcard("*", ""));
    assert!(!wildcard("a*b", "ac"));
}

#[test]
fn track_files_group_in_given_order() {
    let files: Vec<Box<str>> = ["a_wt_1", "a_mut_1", "b_wt_2"].map(Into::into).to_vec();
    let one = track_files(&files, &[]).unwrap();
    assert_eq!(one.len(), 1);
    assert_eq!(one[0].label.as_ref(), "matrix");
    let specs: Vec<Box<str>> = ["mut=*_mut_*", "wt=*_wt_*"].map(Into::into).to_vec();
    let g = track_files(&files, &specs).unwrap();
    assert_eq!(g[0].label.as_ref(), "mut");
    assert_eq!(g[0].files.len(), 1);
    assert_eq!(g[1].files.len(), 2);
    let bad: Vec<Box<str>> = ["wt=*_wt_*"].map(Into::into).to_vec();
    assert!(
        track_files(&files, &bad).is_err(),
        "a_mut_1 matches no track"
    );
    let empty: Vec<Box<str>> = ["wt=*", "none=zzz"].map(Into::into).to_vec();
    assert!(track_files(&files, &empty).is_err(), "a track with no file");
    let malformed: Vec<Box<str>> = ["wt"].map(Into::into).to_vec();
    assert!(track_files(&files, &malformed).is_err());
}

#[test]
fn exact_selector_matches_only_its_row_key() {
    let s = Selector::exact("ENSG1_GENE1");
    assert!(s.matches_gene("ENSG1_GENE1"));
    assert!(!s.matches_gene("ENSG2_GENE2"));
    assert!(!s.matches_gene("ENSG1_GENE10"));
}

#[test]
fn queries_parse_as_locus_or_gene() {
    match parse_query(" chr1:1,000-2,000 ") {
        Some(Query::Locus(r, false)) => {
            assert_eq!((r.chr.as_ref(), r.lb, r.ub), ("chr1", 1000, 2000))
        }
        _ => panic!("window"),
    }
    match parse_query("chr2:1,500") {
        Some(Query::Locus(r, true)) => assert_eq!((r.lb, r.ub), (1500, 1500)),
        _ => panic!("position"),
    }
    assert!(matches!(parse_query("GENE1"), Some(Query::Gene(g)) if &*g == "GENE1"));
    assert!(parse_query("chr1:abc").is_none());
    assert!(parse_query("  ").is_none());
}

#[test]
fn channel_rows_carry_their_channel() {
    let row = parse_row_channel("ENSG1_GENE1/m6a/chr1:100/unmethylated");
    assert_eq!(row, Some(("ENSG1_GENE1", "m6a", "chr1", 100, Some(false))));
    let row = parse_row_channel("ENSG1_GENE1/m6a/chr1:100/methylated");
    assert_eq!(row.map(|r| r.4), Some(Some(true)));
    let row = parse_row_channel("ENSG1_GENE1/m6A/chr1:100");
    assert_eq!(row.map(|r| r.4), Some(None));
}

#[test]
fn merged_positions_sum_shared_sites() {
    let a = [(10, 1.0), (20, 2.0)];
    let b = [(5, 4.0), (20, 3.0), (30, 1.0)];
    assert_eq!(
        merge_positions(&a, &b),
        vec![(5, 4.0), (10, 1.0), (20, 5.0), (30, 1.0)]
    );
}

#[test]
fn channels_are_named_by_modality() {
    assert_eq!(channel_names("m6a"), ("methylated", "unmethylated"));
    assert_eq!(channel_names("m6A"), ("methylated", "unmethylated"));
    assert_eq!(channel_names("atoi"), ("converted", "unconverted"));
    assert_eq!(channel_names("AtoI"), ("converted", "unconverted"));
}

#[test]
fn only_the_opened_gene_is_drawn() {
    use crate::site_analysis::miami::genemodel::GeneModel;
    let model = |symbol: &str| GeneModel {
        chr: "chr1".into(),
        lo: 0,
        hi: 10,
        forward: true,
        exons: Vec::new(),
        symbol: symbol.into(),
    };
    let genes = [model("GENE1"), model("GENE2")];
    let one = genes_to_draw(&genes, "ENSG1_GENE2");
    assert_eq!(one.len(), 1);
    assert_eq!(&*one[0].symbol, "GENE2");
    assert_eq!(
        genes_to_draw(&genes, "chr1:0-100").len(),
        2,
        "a locus shows all"
    );
}

#[test]
fn depth_rows_parse_as_bins() {
    assert_eq!(parse_depth_row("chr1:0-50000"), Some(("chr1", 0, 50_000)));
    assert_eq!(
        parse_depth_row("GL000008.2:100000-150000"),
        Some(("GL000008.2", 100_000, 150_000))
    );
    assert_eq!(parse_depth_row("ENSG1_GENE1/m6a/chr1:100/methylated"), None);
}

/// The row names the producers write today (through the shared
/// `feature_row`) are the ones pileup reads: a naming change there must
/// fail here rather than leave the pileup empty.
#[test]
fn reads_the_rows_the_producers_write() {
    use data_beans::aux::feature_rows::{
        feature_row, ATOI, EDITED, M6A, METHYLATED, UNEDITED, UNMETHYLATED,
    };
    let site = |m, ch| feature_row("ENSG1_GENE1", m, ch, Some("chr1:100"));
    let expect = |m, converted| Some(("ENSG1_GENE1", m, "chr1", 100, Some(converted)));
    assert_eq!(parse_row_channel(&site(M6A, METHYLATED)), expect(M6A, true));
    assert_eq!(
        parse_row_channel(&site(M6A, UNMETHYLATED)),
        expect(M6A, false)
    );
    assert_eq!(parse_row_channel(&site(ATOI, EDITED)), expect(ATOI, true));
    assert_eq!(
        parse_row_channel(&site(ATOI, UNEDITED)),
        expect(ATOI, false)
    );
    assert_eq!(channel_names(M6A), ("methylated", "unmethylated"));
    assert_eq!(channel_names(ATOI), ("converted", "unconverted"));
    // `faba depth` names bins `{chr}:{start}-{end}`.
    let depth = format!("{}:{}-{}", "chr1", 0, 50_000);
    assert_eq!(parse_depth_row(&depth), Some(("chr1", 0, 50_000)));
}
