//! The cell modules' spread over the batches: read off the cell table,
//! shown full pane on `b`, filtered on `B`, printed on `--batches`.

use super::super::text::{printed, Shown};
use super::*;
use crate::cluster::batches::Label;

/// [`tables`] with its cells in four batches, two replicates of `a` and
/// two of `b`: cell module 1 alike in every batch, module 2 in `rep1_a`
/// alone.
fn batched() -> Modules {
    let tmp = tempfile::tempdir().unwrap();
    tables(tmp.path());
    let batches = ["rep1_a", "rep2_a", "rep1_b", "rep2_b"];
    let cells: Vec<Cell> = (0..80)
        .map(|i| Cell {
            barcode: format!("B{i}").into(),
            batch: batches[i / 20].into(),
            fit: true,
        })
        .collect();
    let rows: Vec<write::RoleRow> = (cells.iter().enumerate())
        .map(|(i, cell)| write::RoleRow {
            cell,
            membership: Some(if i < 20 { &[0.5, 0.5] } else { &[0.5, 0.0] }),
            level: None,
        })
        .collect();
    let path = tmp.path().join("m6a_cell_roles.parquet");
    write::cell_roles(&rows, 2, &path).unwrap();
    Modules::load(tmp.path(), "m6a", None).unwrap()
}

#[test]
fn the_cell_table_gives_each_modules_spread() {
    let m = batched();
    let b = m.batches.as_ref().expect("four batches");
    assert_eq!(b.design.conditions, ["a", "b"]);
    assert_eq!(b.spread[0].label, Label::Shared);
    assert_eq!(b.spread[1].label, Label::Batch(vec!["rep1_a".into()]));
    // One batch: nothing to spread over.
    assert!(loaded().batches.is_none());
}

/// `--batch-groups` regroups the batches, and is refused naming a batch
/// the fit does not have.
#[test]
fn batch_groups_regroup_and_refuse_a_batch_not_there() {
    let mut m = batched();
    let groups = [
        ("rep1_a", "x"),
        ("rep2_a", "x"),
        ("rep1_b", "x"),
        ("rep2_b", "x"),
    ]
    .map(|(b, g)| (b.to_string(), g.to_string()))
    .into();
    m.regroup(&groups).unwrap();
    assert_eq!(m.batches.as_ref().unwrap().design.conditions, ["x"]);
    let typo = [("rep1_c".to_string(), "x".to_string())].into();
    let err = m.regroup(&typo).unwrap_err().to_string();
    assert!(
        err.contains("rep1_c") && err.contains("rep1_a, rep1_b"),
        "{err}"
    );
}

/// `b` opens the pane, a row a module with its label; `B` there or in the
/// panes shows the shared modules alone, then the batch-specific ones,
/// then all; Esc or `b` closes the pane.
#[test]
fn b_shows_the_spread_and_capital_b_filters_the_panes() {
    let m = batched();
    let mut v = ModulesView::new("test", &m);
    assert!(
        screen(&mut v, 140, 30).contains("b batches"),
        "in the footer"
    );
    for close in [KeyCode::Esc, KeyCode::Char('b')] {
        press(&mut v, KeyCode::Char('b'));
        let s = screen(&mut v, 140, 30);
        assert!(
            s.contains("batches ·") && !s.contains("gene modules ·"),
            "{s}"
        );
        assert!(s.contains("batch: rep1_a") && s.contains("shared"), "{s}");
        assert!(s.contains("a: rep1_a, rep2_a"), "{s}");
        press(&mut v, close);
        assert!(screen(&mut v, 140, 30).contains("gene modules ·"));
    }
    let names = |v: &mut ModulesView| heads(&screen(v, 140, 30));
    assert!(names(&mut v).contains("C1") && names(&mut v).contains("C2"));
    press(&mut v, KeyCode::Char('B'));
    let s = screen(&mut v, 140, 30);
    assert!(s.contains("(not batch-specific)"), "{s}");
    assert!(names(&mut v).contains("C1") && !names(&mut v).contains("C2"));
    press(&mut v, KeyCode::Char('B'));
    assert!(screen(&mut v, 140, 30).contains("(batch-specific)"));
    assert!(!names(&mut v).contains("C1") && names(&mut v).contains("C2"));
    press(&mut v, KeyCode::Char('B'));
    assert!(names(&mut v).contains("C1") && names(&mut v).contains("C2"));
    // One batch: the pane says why it does not open.
    let one = loaded();
    let mut v = ModulesView::new("test", &one);
    assert!(!screen(&mut v, 140, 30).contains("b batches"));
    press(&mut v, KeyCode::Char('b'));
    assert!(screen(&mut v, 140, 30).contains("the cells are of one batch"));
}

#[test]
fn batches_prints_a_line_a_module() {
    let m = batched();
    let shown = Shown {
        batches: true,
        ..Default::default()
    };
    let lines = printed(&m, &shown).unwrap();
    assert!(lines[0].starts_with("batches:"), "{lines:?}");
    assert_eq!(
        lines[1].split('\t').collect::<Vec<_>>(),
        ["module", "label", "evenness", "rep1_a", "rep2_a", "rep1_b", "rep2_b"]
    );
    assert!(
        lines.iter().any(|l| l.starts_with("C2\tbatch: rep1_a\t")),
        "{lines:?}"
    );
    assert!(
        lines.iter().any(|l| l == "condition a: rep1_a, rep2_a"),
        "{lines:?}"
    );
}
