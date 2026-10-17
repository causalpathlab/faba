//! The strip's layout: sites along a gene and pinned genes in rows.

use crate::cluster::browse::view::{columns, pack, segment_width};

#[test]
fn a_gene_spreads_its_sites_over_its_extent() {
    assert_eq!(columns(&[], (0, 10), 10), Vec::<u16>::new());
    assert_eq!(columns(&[5], (5, 5), 10), vec![5]);
    let c = columns(&[100, 150, 300], (100, 300), 12);
    assert_eq!((c[0], c[2]), (1, 10));
    assert!(c[0] < c[1] && c[1] < c[2]);
    // Outside the extent: clipped to its ends, never past the segment.
    assert_eq!(columns(&[0, 999], (100, 300), 12), vec![1, 10]);
    assert!(segment_width(1) >= 10 && segment_width(1_000_000) <= 32);
}

#[test]
fn pinned_genes_pack_into_rows_as_many_as_fit() {
    // 32 columns: three of 10 and two dividers a row.
    assert_eq!(pack(&[10; 5], 32), vec![0..3, 3..5]);
    assert_eq!(pack(&[10; 3], 100), vec![0..3]);
    assert_eq!(
        pack(&[40, 10], 32),
        vec![0..1, 1..2],
        "too wide: one anyway"
    );
    assert!(pack(&[], 32).is_empty());
}
