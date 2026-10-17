//! One gene's full pileup: where it opens again (its reads are in
//! [`super::posterior`]).

use crate::cluster::browse::gene;

#[test]
fn a_gene_opened_again_in_place_keeps_its_scroll() {
    use crate::site_analysis::pileup::tui::Exit;
    let mut place = gene::Place::default();
    let exits = [
        Exit::Locus {
            lo: 1,
            hi: 9,
            at: 5,
        },
        Exit::Search("GENE2".into()),
        Exit::Modality { at: 5 },
    ];
    for (k, exit) in exits.into_iter().enumerate() {
        assert_eq!(place.left(exit, 2 * k + 4), None);
        assert_eq!(place.track_from, 2 * k + 4, "scrolled as it was");
    }
    assert_eq!(place.extent, (1, 9));
    assert_eq!(place.left(Exit::Genes, 0), Some(false));
    assert_eq!(place.left(Exit::Quit, 0), Some(true));
}
