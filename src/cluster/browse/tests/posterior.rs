//! Every level from the posterior: a gene's full pileup sums its sites'
//! reads over every gene role (from `{m}_site_cells`), split across cell
//! modules by the links' role draws; the lists above it hold each gene in
//! its main gene role only; the cell modules' sizes sum the fitted cells'
//! memberships.

use super::{count, site, write_reads};
use crate::cluster::browse::{gene, report, Count, Modules};
use crate::cluster::read::{Cell, Site};
use crate::cluster::write;
use crate::common::{column, read_record_batch};
use arrow::array::{Float32Array, Int32Array, Int64Array, StringArray};
use std::path::Path;

/// Two cell roles and two gene roles, each block's posterior-mean reads
/// `(converted, unconverted)` at each site, `SPLIT[s][c][g]`: GENE1 (main
/// gene role 1) at 100 and 200 with reads in both gene roles, one entry
/// under a read (200's all-cells reads round apart from their cell
/// modules' rounded), and at 300 with entries under a read only; GENE2
/// (main gene role 2) at 50.
const SPLIT: [[[(f64, f64); 2]; 2]; 4] = [
    [[(4.0, 2.0), (2.0, 2.0)], [(1.0, 3.0), (0.3, 0.1)]],
    [[(0.54, 0.5), (0.0, 0.0)], [(1.04, 0.0), (2.0, 0.0)]],
    [[(0.2, 0.1), (0.0, 0.0)], [(0.0, 0.0), (0.1, 0.1)]],
    [[(0.0, 0.0), (3.0, 1.0)], [(0.0, 0.0), (1.0, 1.0)]],
];

/// The sites of [`SPLIT`].
fn split_sites() -> Vec<Site> {
    vec![
        site("GENE1", 100),
        site("GENE1", 200),
        site("GENE1", 300),
        site("GENE2", 50),
    ]
}

/// Cell memberships over the two cell roles: B0 and B1 fitted, B1 mixed;
/// B2 placed, not fitted.
const SPLIT_MEMBERSHIP: [[f32; 2]; 3] = [[0.75, 0.25], [0.4, 0.6], [1.0, 0.0]];

/// [`SPLIT`] site-major over its four blocks `c * 2 + g`, as the sampler
/// tallies it.
fn split_reads() -> (Vec<f64>, Vec<f64>) {
    let (blocks, roles) = (4, 2);
    let mut on = vec![0.0; SPLIT.len() * blocks];
    let mut off = vec![0.0; SPLIT.len() * blocks];
    for (s, by_cell) in SPLIT.iter().enumerate() {
        for (c, by_gene) in by_cell.iter().enumerate() {
            for (g, &(m, u)) in by_gene.iter().enumerate() {
                on[s * blocks + c * roles + g] = m;
                off[s * blocks + c * roles + g] = u;
            }
        }
    }
    (on, off)
}

/// The tables of [`SPLIT`] and [`SPLIT_MEMBERSHIP`]; with `old`, as older
/// runs wrote them: no `{m}_site_cells`, and no `fitted` column in the
/// cell roles.
fn split_tables(dir: &Path, old: bool) {
    let (on, off) = split_reads();
    write_reads(dir, &split_sites(), (&on, &off), (4, 2), !old);
    let (blocks, roles) = (4, 2);
    let path = |name: &str| dir.join(format!("m6a_{name}.parquet"));
    let shift = [0.0, 1.0, 0.5, -1.0];
    write::blocks(&shift, &[0.1; 4], &[1.0; 4], roles, &path("blocks")).unwrap();
    let weights = [0.7, 0.3, 0.2, 0.8];
    write::gene_roles(&["GENE1", "GENE2"], &weights, roles, &path("gene_roles")).unwrap();
    let cells: Vec<Cell> = (0..3)
        .map(|i| Cell {
            barcode: format!("B{i}").into(),
            batch: "b".into(),
            fit: i < 2,
        })
        .collect();
    let rows: Vec<write::CommunityRow> = cells
        .iter()
        .zip(&SPLIT_MEMBERSHIP)
        .map(|(cell, m)| write::CommunityRow {
            cell,
            membership: Some(m.as_slice()),
            level: None,
        })
        .collect();
    write::communities(&rows, 2, "role", &path("cell_roles")).unwrap();
    if old {
        without_column(&path("cell_roles"), "fitted");
    }
    // B0's link in one block; B1's spread under half a link in every block.
    let fitted = &cells[..2];
    let mut links = vec![0.0; fitted.len() * blocks];
    links[0] = 1.0;
    links[blocks..blocks + 3].copy_from_slice(&[0.4, 0.3, 0.3]);
    let labels = [
        ("barcode", fitted.iter().map(|c| &*c.barcode).collect()),
        ("batch", fitted.iter().map(|c| &*c.batch).collect()),
    ];
    let enrichment = vec![1.0; links.len()];
    write::block_members(
        &labels,
        &links,
        &enrichment,
        blocks,
        roles,
        &path("block_cells"),
    )
    .unwrap();
}

/// The table at `path` rewritten without its column `name`, as older runs
/// wrote it.
fn without_column(path: &Path, name: &str) {
    let mut t = crate::common::read_record_batch(path).unwrap();
    t.remove_column(t.schema().index_of(name).unwrap());
    let file = std::fs::File::create(path).unwrap();
    let mut w = parquet::arrow::ArrowWriter::try_new(file, t.schema(), None).unwrap();
    w.write(&t).unwrap();
    w.close().unwrap();
}

/// The modules of [`split_tables`].
fn split(old: bool) -> Modules {
    let tmp = tempfile::tempdir().unwrap();
    split_tables(tmp.path(), old);
    Modules::load(tmp.path(), "m6a").unwrap()
}

/// Whether `a` and `b` agree to within f32 rounding.
fn close(a: Count, b: Count) -> bool {
    (a.on - b.on).abs() < 1e-5 && (a.total - b.total).abs() < 1e-5
}

/// Site `s`'s reads in cell role `c` over every gene role, from [`SPLIT`].
fn whole(s: usize, c: usize) -> Count {
    SPLIT[s][c].iter().map(|&(m, u)| count(m, m + u)).sum()
}

/// Each cell role's column in `m`.
fn col(m: &Modules, c: usize) -> usize {
    m.columns.iter().position(|&x| x == c).unwrap()
}

/// GENE1's sites, in order, as its members in its main gene role (1).
fn gene1(m: &Modules) -> &[crate::cluster::browse::data::Member] {
    &m.members[0][m.genes(0)[0].clone()]
}

/// The positions of one channel of a gene's reads.
fn positions(channel: &[(i64, f64)]) -> Vec<i64> {
    channel.iter().map(|p| p.0).collect()
}

#[test]
fn counts_sum_on_and_total() {
    let sum: Count = [count(1.0, 2.0), count(0.5, 4.0)].into_iter().sum();
    assert_eq!(sum, count(1.5, 6.0));
    assert_eq!(std::iter::empty::<Count>().sum::<Count>(), Count::default());
}

#[test]
fn site_cells_sum_each_sites_reads_over_gene_roles() {
    let tmp = tempfile::tempdir().unwrap();
    let path = tmp.path().join("site_cells.parquet");
    // Blocks `c * 2 + g`; GENE2 listed first but written after GENE1.
    let sites = [site("GENE2", 5), site("GENE1", 1)];
    let on = [0.0, 5e-4, 0.0, 0.0, 0.3, 1e-4, 0.0, 2.0];
    let off = [0.0, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0];
    write::site_cells(&sites, &on, &off, 4, 2, &path).unwrap();
    let t = read_record_batch(&path).unwrap();
    let names: Vec<String> = t
        .schema()
        .fields()
        .iter()
        .map(|f| f.name().clone())
        .collect();
    assert_eq!(
        names,
        ["gene", "chr", "pos", "cell_role", "converted", "total"]
    );
    let gene: &StringArray = column(&t, "gene").unwrap();
    let pos: &Int64Array = column(&t, "pos").unwrap();
    let role: &Int32Array = column(&t, "cell_role").unwrap();
    let (on, total): (&Float32Array, &Float32Array) = (
        column(&t, "converted").unwrap(),
        column(&t, "total").unwrap(),
    );
    let rows: Vec<(&str, i64, i32, f32, f32)> = (0..t.num_rows())
        .map(|i| {
            let (g, p, r) = (gene.value(i), pos.value(i), role.value(i));
            (g, p, r, on.value(i), total.value(i))
        })
        .collect();
    // GENE2's 5e-4 reads are under the table's threshold: left out.
    assert_eq!(
        rows,
        [("GENE1", 1, 1, 0.3001, 0.4001), ("GENE1", 1, 2, 2.0, 2.0)]
    );
}

#[test]
fn block_sites_keep_only_entries_of_a_read_or_more() {
    let tmp = tempfile::tempdir().unwrap();
    let path = tmp.path().join("block_sites.parquet");
    let sites = [site("GENE1", 1)];
    let (on, off) = ([0.3, 1e-4, 0.0, 2.0], [0.1, 0.0, 0.0, 0.0]);
    write::block_sites(&sites, &on, &off, &[0.0; 4], 4, 2, &path).unwrap();
    let t = read_record_batch(&path).unwrap();
    let total: &Float32Array = column(&t, "total").unwrap();
    let kept: Vec<f32> = total.iter().flatten().collect();
    assert_eq!(kept, [2.0], "near-empty entries never reach the lists");
    // So GENE1 at 300 (entries under a read only) is no member, and GENE1 at
    // 100 holds no 0.4 reads in its main role's lists.
    let m = split(false);
    assert!(gene1(&m).iter().all(|x| m.sites[x.site].pos != 300));
    assert!(m.sites.iter().any(|s| s.pos == 300), "in the full pileup");
}

#[test]
fn the_full_pileup_sums_every_gene_role_per_cell_module() {
    let m = split(false);
    assert!(!m.approximate.pileup);
    let reads = gene::GeneReads::of(&m, "GENE1").unwrap();
    assert_eq!(reads.chr, "chr1");
    assert_eq!(reads.groups.len(), 1 + m.cell_roles);
    let at = [(0, 100), (1, 200), (2, 300)];
    for c in 0..2 {
        let (up, down) = &reads.groups[1 + col(&m, c)];
        assert_eq!(positions(up), [100, 200, 300]);
        for (i, (s, pos)) in at.into_iter().enumerate() {
            let w = whole(s, c);
            assert!((up[i].1 - w.on).abs() < 1e-5, "C{} at {pos}", c + 1);
            assert!((down[i].1 - w.off()).abs() < 1e-5, "C{} at {pos}", c + 1);
        }
    }
    // All cells is the sum over the cell modules.
    let (up, down) = &reads.groups[0];
    for (i, (s, _)) in at.into_iter().enumerate() {
        let all: Count = [whole(s, 0), whole(s, 1)].into_iter().sum();
        assert!((up[i].1 - all.on).abs() < 1e-5 && (down[i].1 - all.off()).abs() < 1e-5);
    }
    assert_eq!(gene::note(&m, &reads), None);
    assert!(gene::GeneReads::of(&m, "GENE9").is_err());
}

#[test]
fn a_genes_share_in_its_main_module_is_a_margin_of_its_full_pileup() {
    let m = split(false);
    let members = gene1(&m);
    assert_eq!(members.len(), 2, "100 and 200 have main-role reads");
    let reads = gene::GeneReads::of(&m, "GENE1").unwrap();
    for (i, member) in members.iter().enumerate() {
        for (c, &[main, other]) in SPLIT[i].iter().enumerate() {
            let k = col(&m, c);
            let main = count(main.0, main.0 + main.1);
            assert!(close(member.reads[k], main), "site {i} C{}", c + 1);
            let sum: Count = [member.reads[k], count(other.0, other.0 + other.1)]
                .into_iter()
                .sum();
            let (up, down) = &reads.groups[1 + k];
            let full = count(up[i].1, up[i].1 + down[i].1);
            assert!(close(sum, full), "site {i} C{}", c + 1);
            assert!(close(full, whole(i, c)));
        }
    }
}

#[test]
fn cell_module_sizes_sum_the_fitted_cells_memberships() {
    let m = split(false);
    assert!(!m.approximate.sizes);
    // B2 was placed, not fitted: not in the posterior tally. B1 counts in
    // full though no block holds half a link of it.
    for c in 0..2 {
        let sum: f64 = SPLIT_MEMBERSHIP[..2].iter().map(|w| w[c] as f64).sum();
        assert!((m.cell_sizes[c] - sum).abs() < 1e-6, "{:?}", m.cell_sizes);
    }
    assert_eq!(m.cells_label(1.5), "1.5");
    // Gene modules count the genes whose main module they are.
    assert_eq!(m.gene_sizes, vec![1, 1]);
}

#[test]
fn older_outputs_load_and_say_what_is_approximate() {
    let m = split(true);
    // Sizes from the cells listed in the block cells, which miss mixed
    // cells: flagged, and shown so.
    assert!(m.approximate.sizes);
    assert_eq!(m.cells_label(1.5), "≈1.5");
    // The full pileup sums the block sites over gene roles, which hold no
    // entry under a read: flagged, and the gene's view says to refit.
    assert!(m.approximate.pileup);
    let reads = gene::GeneReads::of(&m, "GENE1").unwrap();
    assert_eq!(positions(&reads.groups[0].0), [100, 200]);
    // GENE1 at 100, all cells converted: 4 + 2 + 1, its 0.3 under a read lost.
    assert!((reads.groups[0].0[0].1 - 7.0).abs() < 1e-5);
    let note = gene::note(&m, &reads).unwrap();
    assert!(note.contains("--refit"), "{note}");
}

/// The page's level 4 holds the same numbers as the terminal's: per site
/// each cell module's reads over every gene role (`w`, else `c` when they
/// agree, as the page reads them), the sites with none in the main role
/// too, and the posterior sizes.
#[test]
fn the_pages_full_pileup_holds_the_same_posterior_numbers() {
    let tmp = tempfile::tempdir().unwrap();
    split_tables(tmp.path(), false);
    let m = Modules::load(tmp.path(), "m6a").unwrap();
    let data = serde_json::to_value(report::data(&m, &Default::default())).unwrap();
    let reads = gene::GeneReads::of(&m, "GENE1").unwrap();
    let gene = &module(&data, "M1")["genes"][0];
    let full = page_full(gene);
    let lo = gene["lo"].as_i64().unwrap();
    assert_eq!(
        full.iter().map(|s| lo + s.0).collect::<Vec<_>>(),
        positions(&reads.groups[0].0)
    );
    let tenth = |x: f64| (10.0 * x).round() / 10.0;
    for (i, (_, w, wa)) in full.iter().enumerate() {
        for (k, w) in w.iter().enumerate() {
            let (up, down) = &reads.groups[1 + k];
            assert_eq!(*w, [tenth(up[i].1), tenth(up[i].1 + down[i].1)]);
        }
        let (up, down) = &reads.groups[0];
        assert_eq!(*wa, [tenth(up[i].1), tenth(up[i].1 + down[i].1)]);
    }
    for k in 0..m.cell_roles {
        let size = m.cells_label(m.cell_sizes[m.columns[k]]);
        assert_eq!(data["cellLabels"][k], size.as_str());
    }
}

/// The page's all-cells numbers are the exact sums rounded once, as the
/// terminal shows them, not sums of numbers already rounded; it is sent
/// only where the page could not add them up to the same.
#[test]
fn the_pages_all_cells_numbers_round_the_exact_sums() {
    let tmp = tempfile::tempdir().unwrap();
    split_tables(tmp.path(), false);
    let m = Modules::load(tmp.path(), "m6a").unwrap();
    let data = serde_json::to_value(report::data(&m, &Default::default())).unwrap();
    // 1.15 and 0.85 cells: rounded first they would sum to 2.1.
    let all: f64 = m.cell_sizes.iter().sum();
    assert_eq!(data["allLabel"], m.cells_label(all).as_str());
    assert_eq!(data["allLabel"], "2.0");
    let tenth = |x: f64| (10.0 * x).round() / 10.0;
    let sites = module(&data, "M1")["genes"][0]["sites"].as_array().unwrap();
    for (site, member) in sites.iter().zip(gene1(&m)) {
        let sum: Count = member.reads.iter().copied().sum();
        assert_eq!(page_all(site, "ca", "c"), [tenth(sum.on), tenth(sum.total)]);
    }
    // GENE1 at 200: 1.04 reads in each cell module, 2.1 in all, not 2.0.
    assert_eq!(page_all(&sites[1], "ca", "c")[1], 2.1);
}

/// What the page can work out is left out of its data: a site's reads over
/// every gene role where they are its main role's, the all-cells sums
/// where adding the cell modules' gives the same.
#[test]
fn the_page_leaves_out_what_it_can_work_out() {
    let tmp = tempfile::tempdir().unwrap();
    split_tables(tmp.path(), false);
    let m = Modules::load(tmp.path(), "m6a").unwrap();
    let data = serde_json::to_value(report::data(&m, &Default::default())).unwrap();
    // GENE2 at 50: every read in its main role, sums exact.
    let site = &module(&data, "M2")["genes"][0]["sites"][0];
    for key in ["w", "ca", "wa", "o"] {
        assert!(site.get(key).is_none(), "{key} in {site}");
    }
    // GENE1 at 100 has reads in the other gene role.
    let site = &module(&data, "M1")["genes"][0]["sites"][0];
    assert!(site.get("w").is_some(), "{site}");
}

/// A site whose reads all fall in other gene roles is no member of its
/// gene's module (levels 1 to 3), yet the full pileup keeps its reads.
#[test]
fn a_site_only_in_other_gene_roles_is_in_the_full_pileup_alone() {
    let tmp = tempfile::tempdir().unwrap();
    split_tables(tmp.path(), false);
    // Blocks `c * 2 + g`; GENE1 at 200 only in gene role 2.
    let (mut on, mut off) = (vec![0.0; 8], vec![0.0; 8]);
    (on[0], off[0]) = (1.0, 1.0);
    (on[5], off[5]) = (2.0, 1.0);
    let sites = [site("GENE1", 100), site("GENE1", 200)];
    write_reads(tmp.path(), &sites, (&on, &off), (4, 2), true);
    let m = Modules::load(tmp.path(), "m6a").unwrap();
    assert_eq!(gene1(&m).len(), 1, "only 100 listed");
    let reads = gene::GeneReads::of(&m, "GENE1").unwrap();
    assert_eq!(positions(&reads.groups[0].0), [100, 200]);
    assert!((reads.groups[0].0[1].1 - 2.0).abs() < 1e-5);
    let data = serde_json::to_value(report::data(&m, &Default::default())).unwrap();
    let full = page_full(&module(&data, "M1")["genes"][0]);
    assert_eq!(full.iter().map(|s| s.0).collect::<Vec<_>>(), [50, 150]);
}

/// A gene module none of whose genes has reads in it is not listed, so no
/// empty row reaches the grids, the diagonal order or the page; its genes'
/// full pileups still open.
#[test]
fn a_gene_role_without_its_genes_reads_is_not_listed() {
    let tmp = tempfile::tempdir().unwrap();
    split_tables(tmp.path(), false);
    // GENE2 (main gene role 2) with reads only in gene role 1.
    let (mut on, mut off) = split_reads();
    for c in 0..2 {
        let at = 3 * 4 + c * 2;
        (on[at], off[at]) = (on[at + 1], off[at + 1]);
        (on[at + 1], off[at + 1]) = (0.0, 0.0);
    }
    write_reads(tmp.path(), &split_sites(), (&on, &off), (4, 2), true);
    let m = Modules::load(tmp.path(), "m6a").unwrap();
    assert!(m.members[1].is_empty());
    assert_eq!((m.order.clone(), m.diagonal.clone()), (vec![0], vec![0]));
    assert_eq!(m.genes(0).len(), 1, "GENE2 is not listed in gene role 1");
    let data = serde_json::to_value(report::data(&m, &Default::default())).unwrap();
    assert_eq!(data["modules"].as_array().unwrap().len(), 1);
    let reads = gene::GeneReads::of(&m, "GENE2").unwrap();
    assert_eq!(positions(&reads.groups[0].0), [50]);
}

#[test]
fn a_main_gene_role_past_the_blocks_is_an_error_naming_the_table() {
    let tmp = tempfile::tempdir().unwrap();
    split_tables(tmp.path(), false);
    // A gene roles table of another run: GENE1's main role is its third.
    let path = tmp.path().join("m6a_gene_roles.parquet");
    let weights = [0.1, 0.2, 0.7, 0.8, 0.1, 0.1];
    write::gene_roles(&["GENE1", "GENE2"], &weights, 3, &path).unwrap();
    let e = Modules::load(tmp.path(), "m6a").unwrap_err();
    assert!(format!("{e:#}").contains("m6a_gene_roles"), "{e:#}");
}

/// GENE1 at 100 and 200 on chr1, and at 150 on chr0 (first by name);
/// GENE2 at 50.
fn two_chromosomes(dir: &Path) {
    split_tables(dir, false);
    let mut sites = split_sites();
    sites[2] = Site {
        gene: "GENE1".into(),
        chr: "chr0".into(),
        pos: 150,
    };
    let (mut on, off) = split_reads();
    // A read at chr0:150, in the main role.
    on[2 * 4] = 1.0;
    write_reads(dir, &sites, (&on, &off), (4, 2), true);
}

#[test]
fn a_gene_on_two_chromosomes_keeps_the_one_with_most_sites() {
    let tmp = tempfile::tempdir().unwrap();
    two_chromosomes(tmp.path());
    let m = Modules::load(tmp.path(), "m6a").unwrap();
    let reads = gene::GeneReads::of(&m, "GENE1").unwrap();
    assert_eq!(reads.chr, "chr1");
    assert_eq!(positions(&reads.groups[0].0), [100, 200]);
    let note = gene::note(&m, &reads).unwrap();
    assert!(note.contains("chr0"), "{note}");
    // The page's level 4 leaves the chr0 site out too.
    let data = serde_json::to_value(report::data(&m, &Default::default())).unwrap();
    let gene = &module(&data, "M1")["genes"][0];
    assert_eq!(gene["chr"], "chr1");
    assert_eq!(page_full(gene).len(), 2);
}

/// The page's module `name`.
fn module<'a>(data: &'a serde_json::Value, name: &str) -> &'a serde_json::Value {
    data["modules"]
        .as_array()
        .unwrap()
        .iter()
        .find(|x| x["name"] == name)
        .unwrap()
}

/// A JSON `[a, b]` as numbers.
fn pair(v: &serde_json::Value) -> [f64; 2] {
    [v[0].as_f64().unwrap(), v[1].as_f64().unwrap()]
}

/// A JSON list of `[a, b]` as numbers.
fn pairs(v: &serde_json::Value) -> Vec<[f64; 2]> {
    v.as_array().unwrap().iter().map(pair).collect()
}

/// A site's all-cells pair `key` as the page reads it: sent, else its cell
/// modules' pairs `of` added up in order, as the page's script adds them.
fn page_all(site: &serde_json::Value, key: &str, of: &str) -> [f64; 2] {
    match site.get(key) {
        Some(v) => pair(v),
        None => report::page_sum(&pairs(&site[of])),
    }
}

/// A gene's level 4 as the page reads it, by position: its sites in the
/// full pileup (`w`, else `c`) and those only there (`more`), each with its
/// cell modules' pairs and all cells'.
fn page_full(gene: &serde_json::Value) -> Vec<(i64, Vec<[f64; 2]>, [f64; 2])> {
    let mut out: Vec<(i64, Vec<[f64; 2]>, [f64; 2])> = Vec::new();
    for s in gene["sites"].as_array().unwrap() {
        if s.get("o").is_some() {
            continue;
        }
        let of = if s.get("w").is_some() { "w" } else { "c" };
        let wa = match s.get("wa") {
            Some(v) => pair(v),
            None if of == "c" => page_all(s, "ca", "c"),
            None => page_all(s, "wa", "w"),
        };
        out.push((s["pos"].as_i64().unwrap(), pairs(&s[of]), wa));
    }
    for s in gene["more"].as_array().into_iter().flatten() {
        out.push((
            s["pos"].as_i64().unwrap(),
            pairs(&s["w"]),
            page_all(s, "wa", "w"),
        ));
    }
    out.sort_by_key(|s| s.0);
    out
}
