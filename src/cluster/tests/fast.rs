//! The sampler's `ln Gamma` against the library's.

use crate::cluster::fast;
use statrs::function::gamma::ln_gamma;

/// Arguments from `lo` to `hi`, each 1% past the last.
fn grid(lo: f64, hi: f64) -> Vec<f64> {
    std::iter::successors(Some(lo), |x| Some(x * 1.01))
        .take_while(|&x| x < hi)
        .collect()
}

/// Differences of the sampler's `ln Gamma` at arguments a few apart, which
/// the predictives are, stay within `2e-3` of the library's at every size,
/// from `1e-3` to `1e7`, across the switch to the library's included.
#[test]
fn the_samplers_ln_gamma_differences_stay_near_the_librarys() {
    for x in grid(1e-3, 1e7) {
        for n in [0.4, 1.0, 3.0, 17.5] {
            let want = ln_gamma(x + n) - ln_gamma(x);
            let got = fast::ln_gamma(x + n) - fast::ln_gamma(x);
            assert!(
                (want - got).abs() <= 2e-3,
                "ln Gamma({x} + {n}) - ln Gamma({x}): {got} against {want}"
            );
        }
    }
}
