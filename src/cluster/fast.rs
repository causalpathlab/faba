//! A fast, approximate `ln Gamma`: below [`CUTOFF`] Paul Mineiro's
//! `fastlgamma` (fastapprox) in f32, returned as f64 so that sums of its
//! terms are taken in f64; from [`CUTOFF`] on the library's. About `1e-3`
//! error per term: the sampler's target at the sites that take it (see the
//! unified model) is the Beta-binomial through this `ln Gamma`, and the
//! moves are exact for that target.
//!
//! `ln Gamma(x) ~ -2.0810615 - x + 1/(12(x+3)) - ln(x (x+1) (x+2)) + (x +
//! 2.5) ln(x+3)`, every log by `fastlog` (the float's bits read as its
//! `log2`, corrected).

use statrs::function::gamma;

/// Arguments from which the library's `ln Gamma` is taken.
const CUTOFF: f64 = 16.0;

/// Mineiro's `fastlog`: `ln x` for a positive normal `x`, about `1e-4`.
#[inline(always)]
fn fast_ln(x: f32) -> f32 {
    let i = x.to_bits();
    let m = f32::from_bits((i & 0x007F_FFFF) | 0x3f00_0000);
    let y = i as f32 * 1.192_092_9e-7;
    std::f32::consts::LN_2 * (y - 124.225_52 - 1.498_030_3 * m - 1.725_88 / (0.352_088_7 + m))
}

/// Mineiro's `fastlgamma`, for `x > 0`.
#[inline(always)]
fn fast_ln_gamma(x: f32) -> f32 {
    let xp3 = 3.0 + x;
    -2.081_061_5 - x + 0.083_333_3 / xp3 - fast_ln(x * (1.0 + x) * (2.0 + x))
        + (2.5 + x) * fast_ln(xp3)
}

/// `ln Gamma(x)` for `x > 0`.
#[inline]
pub(super) fn ln_gamma(x: f64) -> f64 {
    if x < CUTOFF {
        f64::from(fast_ln_gamma(x as f32))
    } else {
        gamma::ln_gamma(x)
    }
}

/// `ln B(a, b)` through [`ln_gamma`].
pub(super) fn ln_beta(a: f64, b: f64) -> f64 {
    ln_gamma(a) + ln_gamma(b) - ln_gamma(a + b)
}
