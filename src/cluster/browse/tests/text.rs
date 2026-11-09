//! The levels as text: a gene module's genes and sites, a gene's sites in
//! every cell module, the cell modules, the settings; asked for by name,
//! wrong names refused.

use super::super::text::{modules, printed, Shown};
use super::{counted, loaded};

fn shown(f: impl FnOnce(&mut Shown)) -> Shown {
    let mut s = Shown::default();
    f(&mut s);
    s
}

/// A gene module is named as M1 or 1, its rows each a site of its genes
/// with a % and the reads in every cell module; one past the last, or not
/// a module, is refused saying which there are.
#[test]
fn a_gene_module_prints_its_sites_in_every_cell_module() {
    let m = loaded();
    for name in ["M1", "m1", "1"] {
        let lines = printed(&m, &shown(|s| s.module = Some(name.into()))).unwrap();
        assert!(lines[0].starts_with("gene module M1:"), "{lines:?}");
        let head: Vec<&str> = lines[1].split('\t').collect();
        assert_eq!(&head[..2], ["gene", "site"]);
        assert_eq!(head.len(), 2 + 2 * m.cell_roles, "{head:?}");
        assert!(
            head[2].ends_with(" %") && head[3].ends_with(" reads"),
            "{head:?}"
        );
        assert!(lines[2..]
            .iter()
            .all(|l| l.split('\t').count() == head.len()));
        assert!(
            lines[2].split('\t').nth(1).unwrap().contains(':'),
            "a site's chr:pos"
        );
    }
    let past = format!("M{}", m.gene_roles + 1);
    for bad in [past.as_str(), "x"] {
        let refused = printed(&m, &shown(|s| s.module = Some(bad.into()))).unwrap_err();
        assert!(refused.to_string().contains(bad), "{refused}");
    }
}

/// A count fit's gene is found by name in any case: its main gene module,
/// no place, and its UMIs per cell in each cell module (G3: none in the
/// first cell role's two cells, 6 over the second's two, so 3 per cell).
#[test]
fn a_count_fits_gene_prints_its_umis_per_cell() {
    let m = counted();
    let lines = printed(&m, &shown(|s| s.gene = Some("g3".into()))).unwrap();
    assert_eq!(lines[0], "gene G3: main gene module M2");
    let head: Vec<&str> = lines[1].split('\t').collect();
    assert!(
        head[1].ends_with(" UMI/cell") && head[2].ends_with(" cells"),
        "{head:?}"
    );
    let row: Vec<&str> = lines[2].split('\t').collect();
    assert_eq!(row[0], "-");
    let at = |c: usize| {
        let k = m.columns.iter().position(|&x| x == c).unwrap();
        row[1 + 2 * k]
    };
    assert_eq!((at(0), at(1)), ("0.00", "3.0"), "{row:?}");
    let missing = printed(&m, &shown(|s| s.gene = Some("G9".into()))).unwrap_err();
    assert!(missing.to_string().contains("no gene G9"), "{missing}");
}

/// The cell modules print one each with its size; the settings say there is
/// no record when the tables have none; parts asked for together come in
/// the browser's order, a blank line between them.
#[test]
fn the_cell_modules_and_settings_print_in_order() {
    let m = loaded();
    let s = shown(|s| {
        s.cells = true;
        s.settings = true;
    });
    assert!(s.any() && !Shown::default().any());
    // None asked: the modules (level 1).
    assert_eq!(printed(&m, &Shown::default()).unwrap(), modules(&m));
    let lines = printed(&m, &s).unwrap();
    assert!(lines[0].starts_with("cell modules"), "{lines:?}");
    assert_eq!(lines[1], "module\tcells");
    assert_eq!(lines[2 + m.cell_roles], "");
    assert!(
        lines[3 + m.cell_roles].starts_with("settings: "),
        "{lines:?}"
    );
}
