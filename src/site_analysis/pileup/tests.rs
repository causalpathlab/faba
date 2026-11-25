use super::*;

#[test]
fn parse_site_and_mixture_rows() {
    // site output: gene/modality/chr:pos
    assert_eq!(
        parse_row_name_full("ENSG00000000003_GENE3/m6A/chr13:32350000"),
        Some(("ENSG00000000003_GENE3", "m6A", "chr13", 32350000))
    );
    // mixture output: gene/modality/component
    assert_eq!(
        parse_row_name_full("ENSG00000000002_GENE1/m6A/0"),
        Some(("ENSG00000000002_GENE1", "m6A", "", 0))
    );
    // count rows (detail neither chr:pos nor an integer) don't parse
    assert_eq!(parse_row_name_full("gene_0/count/spliced"), None);
}

#[test]
fn relaxed_gene_matching() {
    let gp = "ENSG00000000002_GENE1";
    // symbol, any case
    assert!(gene_matches("GENE1", &query_symbol("GENE1"), gp));
    assert!(gene_matches("gene1", &query_symbol("gene1"), gp));
    // Ensembl ID
    assert!(gene_matches(
        "ENSG00000000002",
        &query_symbol("ENSG00000000002"),
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
    let sel = Selector::build(&["GENE1".into()], &["chr17:100-200".into()]).unwrap();
    // gene branch (mixture row: no chr, component ordinal)
    assert!(sel.selects("ENSG00000000002_GENE1", "", 0));
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
    one.insert("ENSG1_GENE1".into(), 2);
    assert_eq!(summarize_genes(&one).as_ref(), "ENSG1_GENE1");

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
fn channels_are_named_by_modality() {
    assert_eq!(channel_names("m6a"), ("methylated", "unmethylated"));
    assert_eq!(channel_names("m6A"), ("methylated", "unmethylated"));
    assert_eq!(channel_names("atoi"), ("converted", "unconverted"));
    assert_eq!(channel_names("AtoI"), ("converted", "unconverted"));
}

#[test]
fn only_the_opened_genes_are_drawn() {
    use crate::site_analysis::miami::genemodel::GeneModel;
    let model = |symbol: &str, key: &str| GeneModel {
        chr: "chr1".into(),
        lo: 0,
        hi: 10,
        forward: true,
        exons: Vec::new(),
        symbol: symbol.into(),
        key: key.into(),
    };
    let genes = [
        model("GENE1", "ID1_GENE1"),
        model("GENE2", "ID2_GENE2"),
        model("GENE3", "ID3.4_GENE3"),
    ];
    let keys = |k: &[&str]| k.iter().map(|&k| Box::from(k)).collect::<Vec<Box<str>>>();
    let one = genes_to_draw(&genes, Some(&keys(&["ID2_GENE2"])));
    assert_eq!(one.len(), 1);
    assert_eq!(&*one[0].symbol, "GENE2");
    let two = genes_to_draw(&genes, Some(&keys(&["ID1_GENE1", "ID3_GENE3"])));
    assert_eq!(
        two.len(),
        2,
        "every matched gene, by symbol when ids differ"
    );
    assert_eq!(genes_to_draw(&genes, None).len(), 3, "a locus shows all");
}

#[test]
fn depth_rows_parse_as_bins() {
    assert_eq!(parse_depth_row("chr1:0-50000"), Some(("chr1", 0, 50_000)));
    assert_eq!(
        parse_depth_row("GL000008.2:100000-150000"),
        Some(("GL000008.2", 100_000, 150_000))
    );
    assert_eq!(parse_depth_row("ENSG1_GENE1/m6a/chr1:100/methylated"), None);
    // Alt contigs keep their underscores; other spellings are not loci.
    assert_eq!(
        parse_depth_row("chr1_KI270706v1_random:0-50000"),
        Some(("chr1_KI270706v1_random", 0, 50_000))
    );
    assert_eq!(parse_depth_row("chr1_0_50000"), None);
    assert_eq!(parse_depth_row("chr1:50000-0"), None);
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

#[test]
fn a_gene_symbol_may_hold_a_slash() {
    let row = parse_row_channel("ID1_GENE1/B/m6a/chr1:100/methylated");
    assert_eq!(row, Some(("ID1_GENE1/B", "m6a", "chr1", 100, Some(true))));
    assert_eq!(parse_row_channel("ID1_GENE1/m6a/chr1:100/other"), None);
}

#[test]
fn an_alt_contig_position_is_not_read_as_a_gene() {
    // A BAF row's first segment is a position on an alt contig: no gene
    // query may match a piece of it.
    let gp = "chr1_KI270706v1_random:12345";
    for q in ["random:12345", "random", "KI270706v1", "chr1"] {
        assert!(!gene_matches(q, &query_symbol(q), gp), "{q}");
    }
    assert!(gene_matches(gp, &query_symbol(gp), gp));
    // The shared gene rule leaves a coordinate whole.
    assert_eq!(&*query_symbol(gp), gp);
}

/// `faba pileup` with `args`, as its command line would give them.
fn pileup_args(args: &[&str]) -> PileupArgs {
    use clap::FromArgMatches;
    let cmd = crate::faba_command()
        .find_subcommand("pileup")
        .unwrap()
        .clone();
    let m = cmd
        .try_get_matches_from(std::iter::once("pileup").chain(args.iter().copied()))
        .unwrap();
    PileupArgs::from_arg_matches(&m).unwrap()
}

#[test]
fn an_output_directory_gives_a_track_per_batch_and_its_site_table() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    for f in [
        "wt_m6a_site.zarr.zip",
        "mut_m6a_site.zarr.zip",
        "wt_m6a.zarr.zip",
        "m6a_sites.parquet",
        "pipeline_summary.json",
    ] {
        std::fs::write(dir.join(f), b"").unwrap();
    }
    let path = |f: &str| dir.join(f).to_string_lossy().into_owned();
    for given in [
        dir.to_string_lossy().into_owned(),
        path("pipeline_summary.json"),
    ] {
        let (args, groups) = resolve(&[&given]);
        let files: Vec<&str> = args.data_files.iter().map(|f| &**f).collect();
        assert_eq!(
            files,
            [path("mut_m6a_site.zarr.zip"), path("wt_m6a_site.zarr.zip")],
            "site matrices only"
        );
        assert_eq!(labels(&groups), ["mut", "wt"], "a track per batch");
        assert_eq!(args.site_file.as_deref(), Some(&*path("m6a_sites.parquet")));
    }
    let d = dir.to_string_lossy().into_owned();
    // A site table given keeps its place.
    let (args, _) = resolve(&[&d, "-s", "x.parquet"]);
    assert_eq!(args.site_file.as_deref(), Some("x.parquet"));
    // A matrix named beside the directory has its own track.
    let (_, groups) = resolve(&[&d, "extra_m6a_site.zarr.zip"]);
    assert_eq!(
        labels(&groups),
        ["mut", "wt", "extra"],
        "named after its batch"
    );
    // --track regroups them all.
    let (_, groups) = resolve(&[&d, "--track", "all=*"]);
    assert_eq!(labels(&groups), ["all"]);
    assert_eq!(groups[0].files.len(), 2);
    // Matrices named one by one make one track, as before.
    let one = path("wt_m6a_site.zarr.zip");
    let (args, groups) = resolve(&[&one]);
    assert_eq!(args.data_files, vec![Box::from(one.as_str())]);
    assert_eq!(labels(&groups), ["matrix"]);
    // --separate: a track each, by batch.
    let two = [path("wt_m6a_site.zarr.zip"), path("mut_m6a_site.zarr.zip")];
    let (_, groups) = resolve(&[&two[0], &two[1], "--separate"]);
    assert_eq!(labels(&groups), ["wt", "mut"]);
}

/// `faba pileup` with `args`, resolved: the inputs and the tracks.
fn resolve(args: &[&str]) -> (PileupArgs, Vec<TrackFiles>) {
    let args = pileup_args(args);
    args.resolved().unwrap()
}

fn labels(groups: &[TrackFiles]) -> Vec<&str> {
    groups.iter().map(|g| &*g.label).collect()
}

/// `faba pileup`'s first two tracks share a row when both are read tracks
/// of one measure, in the plain colours; else every track has its own row.
#[test]
fn the_pileup_mirrors_its_first_two_alike_tracks() {
    let args = pileup_args(&["x_m6a_site.zarr"]);
    let track = |label: &str, total: Option<Vec<(i64, f64)>>| MatrixTrack {
        label: label.into(),
        converted: vec![(10, 1.0)],
        total,
    };
    let loaded = |tracks| Loaded {
        gene: "GENE1".into(),
        chr: "chr1".into(),
        modality: "m6a".into(),
        tracks,
        sites: None,
        extent: (0, 100),
        keys: None,
        depth: Vec::new(),
        cursor: None,
    };
    let flags = |l: &Loaded| -> Vec<(bool, bool)> {
        tracks_of(&args, l, "stacked")
            .iter()
            .map(|t| (t.mirror_with_next, t.accent))
            .collect()
    };
    let alike = loaded(vec![track("a", None), track("b", None), track("c", None)]);
    assert_eq!(
        flags(&alike),
        [(true, false), (false, false), (false, false)]
    );
    // A track with a total is another measure than one without.
    let apart = loaded(vec![track("a", Some(vec![(10, 2.0)])), track("b", None)]);
    assert_eq!(flags(&apart), [(false, false); 2]);
    assert_eq!(flags(&loaded(vec![track("a", None)])), [(false, false)]);
}

/// A site matrix at `dir/name`: `GENE1`'s methylated and unmethylated rows
/// at chr1:100 then chr1:200, over cells `c1` and `c2`, holding `entries`
/// as `(row, cell, value)`.
fn site_matrix(dir: &std::path::Path, name: &str, entries: &[(u64, u64, f32)]) -> Box<str> {
    site_matrix_of(dir, name, "ENSG1_GENE1", "chr1", entries)
}

/// [`site_matrix`]'s rows for gene key `gene` on `chr`.
fn site_matrix_of(
    dir: &std::path::Path,
    name: &str,
    gene: &str,
    chr: &str,
    entries: &[(u64, u64, f32)],
) -> Box<str> {
    use crate::common::{backend_meta, ToBackend, TripletsRowsCols};
    use data_beans::aux::feature_rows::{feature_row, M6A};
    let path = dir.join(name).to_string_lossy().into_owned();
    let row = |pos: i64, ch| feature_row(gene, M6A, ch, Some(&format!("{chr}:{pos}")));
    let triplets = TripletsRowsCols {
        triplets: entries.to_vec(),
        rows: vec![
            row(100, METHYLATED),
            row(100, UNMETHYLATED),
            row(200, METHYLATED),
            row(200, UNMETHYLATED),
        ],
        cols: vec!["c1".into(), "c2".into()],
    };
    drop(
        triplets
            .to_backend(&path, &backend_meta("dartseq", "x"))
            .unwrap(),
    );
    path.into()
}

/// Tracks given by pattern give the figure a panel each, each summing its
/// own matrices' methylated sites.
#[test]
fn pattern_tracks_give_the_figure_a_panel_each() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let wt = [(0, 0, 1.0), (0, 1, 2.0), (1, 0, 7.0), (2, 0, 4.0)];
    let wt = site_matrix(dir, "b1_wt_m6a_site.zarr", &wt);
    let mt = [(0, 0, 10.0), (2, 0, 5.0), (2, 1, 20.0), (3, 1, 9.0)];
    let mt = site_matrix(dir, "b2_mut_m6a_site.zarr", &mt);
    let given: [&[&str]; 2] = [
        &["--track", "wt=*_wt_*", "--track", "mut=*_mut_*"],
        &["--track", "wt=*_wt_*,mut=*_mut_*"],
    ];
    for tracks in given {
        let mut cmd = vec![&*wt, &*mt, "--genes", "GENE1", "--svg", "--batch-process"];
        cmd.extend_from_slice(tracks);
        let (args, groups) = resolve(&cmd);
        let selector = Selector::build(&args.genes, &args.regions).unwrap();
        let sites = figure_sites(&args, &groups, &selector, None).unwrap();
        // One cell type: one colour in every track.
        assert_eq!(
            panel_rows(&sites),
            [
                ("wt", 0, 0, &[(100, 3.0), (200, 4.0)][..]),
                ("mut", 1, 0, &[(100, 10.0), (200, 25.0)][..]),
            ],
            "{tracks:?}"
        );
        assert_eq!(&*sites.chr, "chr1");
    }
}

/// A figure panel's label, track, hue and sites.
type PanelRow<'a> = (&'a str, usize, usize, &'a [(i64, f64)]);

/// The panels' [`PanelRow`]s.
fn panel_rows(sites: &FigureSites) -> Vec<PanelRow<'_>> {
    sites
        .panels
        .iter()
        .map(|p| (&*p.label, p.track, p.hue, &*p.sites))
        .collect()
}

/// Without `--track` the figure pools every file into one track, even
/// where the browser and the ASCII pileup give each its own (`--separate`).
#[test]
fn the_figure_pools_its_files_without_track() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let wt = [(0, 0, 1.0), (0, 1, 2.0), (2, 0, 4.0)];
    let wt = site_matrix(dir, "b1_wt_m6a_site.zarr", &wt);
    let mt = [(0, 0, 10.0), (2, 0, 5.0), (2, 1, 20.0)];
    let mt = site_matrix(dir, "b2_mut_m6a_site.zarr", &mt);
    let cmd = [
        &*wt,
        &*mt,
        "--separate",
        "--genes",
        "GENE1",
        "--svg",
        "--batch-process",
    ];
    let (args, groups) = resolve(&cmd);
    assert_eq!(
        labels(&groups),
        ["b1_wt", "b2_mut"],
        "--separate keeps its tracks"
    );
    let tracks = figure_tracks(&args, &groups).unwrap();
    let selector = Selector::build(&args.genes, &args.regions).unwrap();
    let sites = figure_sites(&args, &tracks, &selector, None).unwrap();
    assert_eq!(
        panel_rows(&sites),
        [("", 0, 0, &[(100, 13.0), (200, 29.0)][..])]
    );
}

/// The figure names the genes and chromosome of every track's files, not
/// the first track's alone.
#[test]
fn the_figure_is_named_by_every_track() {
    let tmp = tempfile::tempdir().unwrap();
    let dir = tmp.path();
    let one = [(0, 0, 1.0)];
    let wt = site_matrix_of(dir, "b1_wt_m6a_site.zarr", "ENSG1_GENE1", "chr1", &one);
    let mt = site_matrix_of(dir, "b2_mut_m6a_site.zarr", "ENSG2_GENE2", "chr2", &one);
    let cmd = [
        &*wt,
        &*mt,
        "--genes",
        "GENE1,GENE2",
        "--svg",
        "--batch-process",
        "--track",
        "wt=*_wt_*,mut=*_mut_*",
    ];
    let (args, groups) = resolve(&cmd);
    let selector = Selector::build(&args.genes, &args.regions).unwrap();
    let sites = figure_sites(&args, &groups, &selector, None).unwrap();
    assert_eq!(&*sites.gene, "2 genes: ENSG1_GENE1,ENSG2_GENE2");
    assert_eq!(&*sites.chr, "*");
}

/// A comma splits `--track` only where another `LABEL=` starts.
#[test]
fn a_track_pattern_may_hold_a_comma() {
    assert_eq!(split_tracks("a=*rep1,2*"), ["a=*rep1,2*"]);
    assert_eq!(
        split_tracks("wt=*_wt_*,mut=*_mut_*"),
        ["wt=*_wt_*", "mut=*_mut_*"]
    );
    assert_eq!(split_tracks("a=x,y,b=z"), ["a=x,y", "b=z"]);
    let args = pileup_args(&["s_rep1,2_x", "--track", "a=*rep1,2*"]);
    let files = ["s_rep1,2_x".into()];
    let g = track_files(&files, &args.tracks).unwrap();
    assert_eq!(labels(&g), ["a"]);
    assert_eq!(g[0].files.len(), 1);
    let args = pileup_args(&["x", "--track", "wt=*_wt_*,mut=*_mut_*"]);
    let files = ["b1_wt_x".into(), "b2_mut_x".into()];
    let g = track_files(&files, &args.tracks).unwrap();
    assert_eq!(labels(&g), ["wt", "mut"]);
}

/// Each pattern track reads depth from the BAMs its pattern matches; a
/// pooled track from all of them.
#[test]
fn each_track_reads_depth_from_its_own_bams() {
    let files: Vec<Box<str>> = ["b1_wt_x", "b2_mut_x"].map(Into::into).to_vec();
    let bams: Vec<Box<str>> = ["b1_wt.bam", "b2_mut.bam", "b3_wt.bam"]
        .map(Into::into)
        .to_vec();
    let specs: Vec<Box<str>> = ["mut=*_mut*", "wt=*_wt*"].map(Into::into).to_vec();
    let g = track_files(&files, &specs).unwrap();
    let names = |v: Vec<Vec<Box<str>>>| -> Vec<Vec<String>> {
        v.into_iter()
            .map(|b| b.iter().map(|s| s.to_string()).collect())
            .collect()
    };
    assert_eq!(
        names(track_bams(&g, &bams)),
        [vec!["b2_mut.bam"], vec!["b1_wt.bam", "b3_wt.bam"]]
    );
    let none: Vec<Box<str>> = ["b1_wt.bam"].map(Into::into).to_vec();
    assert_eq!(
        names(track_bams(&g, &none)),
        [vec![], vec!["b1_wt.bam".to_string()]]
    );
    let pooled = track_files(&files, &[]).unwrap();
    assert_eq!(names(track_bams(&pooled, &bams)).concat().len(), 3);

    // A panel draws its own track's depth of its cell type.
    let panel = |track: usize, celltype: &str| FigurePanel {
        label: celltype.into(),
        track,
        celltype: celltype.into(),
        hue: 0,
        sites: Vec::new(),
    };
    let depth_of = |ct: &str, v: f64| -> FxHashMap<Box<str>, Vec<f64>> {
        [(Box::from(ct), vec![v; 2])].into_iter().collect()
    };
    let depth = [
        depth_of("CT1", 1.0),
        FxHashMap::default(),
        depth_of("CT1", 3.0),
    ];
    let panels = [panel(0, "CT1"), panel(1, "CT1"), panel(2, "CT1")];
    let drawn = panel_data(panels.into(), &depth, 2, false);
    let bins: Vec<&[f64]> = drawn.iter().map(|p| &*p.depth_bins).collect();
    assert_eq!(bins, [&[1.0, 1.0][..], &[0.0, 0.0], &[3.0, 3.0]]);
}

/// `faba pileup x.zarr` and `line`, or why its command line refuses them.
fn figure_args(line: &[&str]) -> anyhow::Result<PileupArgs> {
    use clap::FromArgMatches;
    let cmd = crate::faba_command()
        .find_subcommand("pileup")
        .unwrap()
        .clone();
    let all = ["pileup", "x.zarr"].into_iter().chain(line.iter().copied());
    Ok(PileupArgs::from_arg_matches(
        &cmd.try_get_matches_from(all)?,
    )?)
}

/// Genes come from -q and a gene file (lines, commas or spaces, `#`
/// comments), each once whatever its case.
#[test]
fn the_genes_are_q_and_the_gene_file_each_once() {
    let dir = tempfile::tempdir().unwrap();
    let file = dir.path().join("genes.txt");
    std::fs::write(&file, "GENE2, GENE3 # a comment\n\n#GENE9\ngene1\tGENE4\n").unwrap();
    let a = figure_args(&["-q", "GENE1", "--gene-file", file.to_str().unwrap()]).unwrap();
    let genes = a.gene_list().unwrap();
    assert_eq!(
        genes,
        ["GENE1", "GENE2", "GENE3", "GENE4"].map(Box::<str>::from)
    );
    let missing = figure_args(&["--gene-file", "/no/such.txt"]).unwrap();
    assert!(missing
        .gene_list()
        .unwrap_err()
        .to_string()
        .contains("/no/such.txt"));
}

/// `--row` lays every gene side by side, `--grid COLSxROWS` COLS to a
/// row; a grid too small for the genes, or not written so, is refused.
#[test]
fn row_and_grid_lay_the_genes_out() {
    let row = figure_args(&["--row"])
        .unwrap()
        .combine(5)
        .unwrap()
        .unwrap();
    assert_eq!((row.cols, row.name), (5, "row"));
    let grid = figure_args(&["--grid", "3x2"])
        .unwrap()
        .combine(5)
        .unwrap()
        .unwrap();
    assert_eq!((grid.cols, grid.name), (3, "grid"));
    assert!(figure_args(&[]).unwrap().combine(5).unwrap().is_none());
    let small = figure_args(&["--grid", "2x2"])
        .unwrap()
        .combine(5)
        .unwrap_err();
    assert!(small.to_string().contains("holds 4"), "{small}");
    for bad in ["3", "0x2", "ax2"] {
        let e = figure_args(&["--grid", bad])
            .unwrap()
            .combine(5)
            .unwrap_err();
        assert!(e.to_string().contains("COLSxROWS"), "{bad}: {e}");
    }
    assert!(
        figure_args(&["--row", "--grid", "2x3"]).is_err(),
        "one or the other"
    );
}

/// The figure's flags are the browser's state: the y scale (log10 the
/// browser's log, log1p naming it too), both channels, the difference row, the key's place.
#[test]
fn the_figures_flags_are_the_browsers_state() {
    use crate::tui::ui::Scale;
    let plain = figure_args(&[]).unwrap().figure_settings();
    assert_eq!(plain.y_scale, Some(Scale::Linear));
    assert_eq!(plain.difference, Some(false));
    assert_eq!(plain.key_at, Some(tui::KeyAt::Inside));
    assert_eq!(plain.show, Some(crate::site_analysis::show::Show::Both));
    let a = figure_args(&[
        "--yscale",
        "log10",
        "--difference",
        "--key-at",
        "alone",
        "--key",
    ])
    .unwrap();
    let s = a.figure_settings();
    assert_eq!(s.y_scale, Some(Scale::Log));
    assert_eq!(s.difference, Some(true));
    assert_eq!(s.key_at, Some(tui::KeyAt::Alone));
    assert!(a.key);
    let sqrt = figure_args(&["--yscale", "sqrt"])
        .unwrap()
        .figure_settings();
    assert_eq!(sqrt.y_scale, Some(Scale::Sqrt));
    let log1p = figure_args(&["--yscale", "log1p"]).unwrap();
    assert_eq!(log1p.figure_settings().y_scale, Some(Scale::Log));
    assert!(
        figure_args(&["--yscale", "log"]).is_err(),
        "log10, by its name"
    );
    // Without --bins the printed pileup keeps its 80, the figure the browser's bars.
    assert_eq!(figure_args(&[]).unwrap().num_bins, None);
}
