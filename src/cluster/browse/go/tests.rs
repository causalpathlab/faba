//! Gene modules' terms on made-up gene sets: the test, the universe, the
//! matching by symbol and alias, and BH.

use super::*;

/// 40 genes `ENSG000000000{i}_G{i}`: module 0 the first 20 (weights 0.9,
/// but 0.3 for the last two), module 1 the rest.
fn genes() -> Vec<(Box<str>, usize, f32)> {
    (0..40)
        .map(|i| {
            let w = if i == 18 || i == 19 { 0.3 } else { 0.9 };
            (format!("ENSG{i:011}_G{i}").into(), usize::from(i >= 20), w)
        })
        .collect()
}

/// A GMT: `UP` holds 12 of module 0's genes and 2 of module 1's, `FLAT`
/// 5 of each, `TINY` 3 genes (too few to test), `OUT` none of the fit's,
/// `FILL` 8 of module 1's alone.
fn gmt(dir: &std::path::Path) -> PathBuf {
    let set = |name: &str, g: &[usize]| {
        let genes: Vec<String> = g.iter().map(|i| format!("G{i}")).collect();
        format!("{name}\t-\t{}", genes.join("\t"))
    };
    let up: Vec<usize> = (0..12).chain([30, 31]).collect();
    let flat: Vec<usize> = (12..17).chain(25..30).collect();
    let lines = [
        set("UP", &up),
        set("FLAT", &flat),
        set("TINY", &[1, 2, 3]),
        "OUT\t-\tZZZ1\tZZZ2\tZZZ3\tZZZ4\tZZZ5".to_string(),
        set("FILL", &(32..40).collect::<Vec<_>>()),
    ];
    let path = dir.join("sets.gmt");
    std::fs::write(&path, lines.join("\n")).unwrap();
    path
}

#[test]
fn a_modules_genes_are_tested_against_the_fits_annotated_genes() {
    let tmp = tempfile::tempdir().unwrap();
    let source = Source::Gmt(gmt(tmp.path()));
    let t = terms(&genes(), 2, &source, 0.5).unwrap();
    // The universe: genes some set holds (UP, FLAT, TINY, FILL): 0..17,
    // 25..40 and 1..3 among them.
    assert_eq!(t.universe, 17 + 15);
    let m0 = t.of(0).unwrap();
    // Genes 18, 19 are outside every set; so the cut changes nothing here.
    assert_eq!((m0.all, m0.n), (17, 17));
    assert_eq!(m0.hits.len(), 1, "{:?}", m0.hits);
    let up = &m0.hits[0];
    assert_eq!(
        (up.id.as_str(), up.k, up.n, up.big_k, up.big_n),
        ("UP", 12, 17, 14, 32)
    );
    assert!(up.q <= Q && up.odds > 1.0);
    assert_eq!(up.genes.len(), 5);
    // Module 1 (M2): FILL, all its own genes; FLAT, even, in neither.
    let m1 = t.of(1).unwrap();
    assert_eq!(
        m1.hits.iter().map(|h| h.id.as_str()).collect::<Vec<_>>(),
        ["FILL"]
    );
    let lines = t.lines();
    assert!(
        lines[0].starts_with("go: gene sets sets.gmt; 32 of the fit's 40 genes"),
        "{lines:?}"
    );
    assert!(
        lines.iter().any(|l| l.starts_with("M1\t17/17\t1\tUP\t")),
        "{lines:?}"
    );
    assert!(
        lines.iter().any(|l| l.starts_with("M2\t15/15\t1\tFILL\t")),
        "{lines:?}"
    );
    // A module with too few genes says so.
    let few = terms(&genes(), 3, &source, 0.5).unwrap();
    assert!(few
        .lines()
        .iter()
        .any(|l| l.starts_with("M3\t0/0\t-\t(fewer than")));
}

/// A gene matched by its symbol keeps its row; an alias another gene has,
/// or one onto a row taken, matches nothing.
#[test]
fn aliases_never_give_one_gene_anothers_terms() {
    let genes: Vec<(Box<str>, usize, f32)> = ["ENSG1_AAA", "ENSG2_BBB", "ENSG3_CCC"]
        .iter()
        .map(|g| ((*g).into(), 0, 1.0))
        .collect();
    let mut sets = GeneSets::default();
    let alias =
        |names: &[&str]| -> FxHashSet<Box<str>> { names.iter().map(|n| (*n).into()).collect() };
    // AAA by symbol; OLD (retired) names BBB among its aliases, which BBB
    // already is; NEW by a unique alias CCC.
    sets.gene_aliases.insert("AAA".into(), alias(&["AAA"]));
    sets.gene_aliases.insert("BBB".into(), alias(&["BBB"]));
    sets.gene_aliases
        .insert("OLD".into(), alias(&["OLD", "BBB"]));
    sets.gene_aliases
        .insert("NEW".into(), alias(&["NEW", "CCC", "SHARED"]));
    sets.gene_aliases
        .insert("TWO".into(), alias(&["TWO", "SHARED"]));
    let m = matched(&sets, &genes);
    assert_eq!(m.get("AAA"), Some(&0));
    assert_eq!(m.get("BBB"), Some(&1));
    assert_eq!(m.get("OLD"), None);
    assert_eq!(m.get("NEW"), Some(&2));
    assert_eq!(m.get("TWO"), None);
}

#[test]
fn bh_is_monotone_and_at_most_one() {
    let p = [0.01f64, 0.04, 0.03, 0.5];
    let q = bh(&p.map(f64::ln));
    assert!((q[0] - 0.04).abs() < 1e-12, "{q:?}");
    assert!(q.iter().all(|&x| x <= 1.0));
    assert!(q[0] <= q[2] && q[2] <= q[1] && q[1] <= q[3]);
}

#[test]
fn the_source_says_what_is_missing() {
    let genes = genes();
    let err = Source::of(Some("a.gaf"), Some("b.gmt"), None, &genes).unwrap_err();
    assert!(err.contains("not both"));
    assert_eq!(
        Source::of(None, Some("b.gmt"), None, &genes).unwrap(),
        Source::Gmt("b.gmt".into())
    );
    let other: Vec<(Box<str>, usize, f32)> = vec![("geneA".into(), 0, 1.0)];
    assert!(Source::of(None, None, None, &other)
        .unwrap_err()
        .contains("not human or mouse"));
}
