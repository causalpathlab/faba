//! The batches read as replicates of conditions, and modules labelled by
//! how they spread over them.

use super::*;

/// Cells per batch and shares summed per module and batch, from each
/// cell's batch and shares.
fn tally(batch_of: &[usize], shares: &[f32], nb: usize, k: usize) -> (Vec<f64>, Vec<f64>) {
    let (mut cells, mut units) = (vec![0.0; nb], vec![0.0; k * nb]);
    for (&b, row) in batch_of.iter().zip(shares.chunks_exact(k)) {
        cells[b] += 1.0;
        for (m, &s) in row.iter().enumerate() {
            units[m * nb + b] += f64::from(s);
        }
    }
    (cells, units)
}

fn names(of: &[&str]) -> Vec<String> {
    of.iter().map(|s| s.to_string()).collect()
}

/// `rep1_mut`, `mut-r2` and `wt.3` lose their replicate token; a design
/// whose conditions have one batch each is one condition of replicates.
#[test]
fn conditions_are_read_off_the_batch_names() {
    let d = Design::of(
        &names(&["rep1_mut", "rep1_wt", "rep2_mut", "rep2_wt"]),
        &BTreeMap::new(),
    );
    assert_eq!(d.conditions, ["mut", "wt"]);
    assert_eq!(d.condition_of, [0, 1, 0, 1]);
    assert_eq!(condition_name("mut-r2"), "mut");
    assert_eq!(condition_name("wt_replicate3"), "wt");
    assert_eq!(condition_name("runA/mut"), "mut");
    // Digits alone are a time point or a dose, not a replicate.
    assert_eq!(condition_name("day_3_rep1"), "day_3");
    assert_eq!(condition_name("wt.3"), "wt.3");
    assert_eq!(condition_name("brain"), "brain");
    // Each batch its own condition: no replicates to judge a module by.
    let own = Design::of(&names(&["liver", "brain"]), &BTreeMap::new());
    assert_eq!(own.conditions, ["brain", "liver"]);
    assert!(!own.replicated());
    let groups = [("a", "x"), ("b", "x"), ("c", "y")]
        .map(|(b, g)| (b.to_string(), g.to_string()))
        .into();
    let given = Design::of(&names(&["a", "b", "c"]), &groups);
    assert_eq!(given.condition_of, [0, 0, 1]);
}

/// Four batches, two of each condition, 40 cells each, and four modules:
/// even everywhere, even over mut's replicates alone, one batch's alone,
/// and mostly wt's but in one of its replicates.
#[test]
fn modules_are_shared_condition_or_batch_specific() {
    let d = Design::of(
        &names(&["rep1_mut", "rep1_wt", "rep2_mut", "rep2_wt"]),
        &BTreeMap::new(),
    );
    let mut batch_of = Vec::new();
    let mut shares = Vec::new();
    for b in 0..4 {
        for _ in 0..40 {
            batch_of.push(b);
            let mutant = d.condition_of[b] == 0;
            let row: [f32; 4] = match b {
                0 => [0.25, 0.5, 0.25, 0.0],
                2 => [0.25, 0.5, 0.0, 0.25],
                1 => [0.25, 0.0, 0.0, 0.75],
                _ => [0.25, 0.0, 0.0, 0.75 * 0.1],
            };
            assert!(mutant == (b % 2 == 0));
            shares.extend(row);
        }
    }
    let s = spread(
        &d,
        &tally(&batch_of, &shares, 4, 4).0,
        &tally(&batch_of, &shares, 4, 4).1,
        4,
    );
    assert_eq!(s[0].label, Label::Shared);
    assert!((s[0].evenness.unwrap() - 1.0).abs() < 1e-9);
    assert_eq!(s[1].label, Label::Condition("mut".into()));
    assert_eq!(s[2].label, Label::Batch(names(&["rep1_mut"])));
    assert!(s[3].label.is_batch(), "{:?}", s[3]);
    for m in &s {
        assert!((m.per_batch.iter().sum::<f64>() - 1.0).abs() < 1e-9);
    }
    // One batch a condition: nothing to tell.
    let one = Design::of(&names(&["rep1_mut"]), &BTreeMap::new());
    let s = spread(&one, &[20.0], &[20.0], 1);
    assert_eq!(s[0].label, Label::Unknown);
    // Too few cells' worth to tell.
    let s = spread(&d, &[10.0; 4], &[1.0; 4], 1);
    assert_eq!(s[0].label, Label::Few);
}

/// A condition of more replicates weighs no more than one of fewer; a
/// module mostly where nothing is replicated is not judged on the rest.
#[test]
fn unequal_designs_weigh_conditions_alike() {
    let batches: Vec<String> = (1..=8)
        .map(|r| format!("rep{r}_mut"))
        .chain(["rep1_wt".into(), "rep2_wt".into()])
        .collect();
    let d = Design::of(&batches, &BTreeMap::new());
    let s = spread(&d, &[100.0; 10], &[50.0; 10], 1);
    assert_eq!(s[0].label, Label::Shared, "{:?}", s[0]);
    let d = Design::of(&names(&["rep1_mut", "rep2_mut", "wt"]), &BTreeMap::new());
    let s = spread(&d, &[100.0; 3], &[5.0, 0.0, 95.0], 1);
    assert_eq!(s[0].label, Label::Unknown, "{:?}", s[0]);
}
