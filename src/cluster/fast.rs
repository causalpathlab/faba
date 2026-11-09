//! The sampler's fast, approximate logs, each returned as f64 so that sums
//! of its terms are taken in f64: `ln` by `legume_numeric::fast_math` (f32,
//! Paul Mineiro's fastapprox) on the mantissa, and an
//! absolute `ln Gamma`, below [`CUTOFF`] Mineiro's `fastlgamma` and from
//! [`CUTOFF`] on the library's. About `1e-4` error per `ln`. The urns'
//! `ln Gamma` differences are sums of such logs (see [`ln_urn`]).
//!
//! `ln Gamma(x) ~ -2.0810615 - x + 1/(12(x+3)) - ln(x (x+1) (x+2)) + (x +
//! 2.5) ln(x+3)`, every log by `fastlog` (the float's bits read as its
//! `log2`, corrected).

use legume_numeric::fast_math::fast_ln;
use statrs::function::gamma;

/// Arguments from which the library's `ln Gamma` is taken.
const CUTOFF: f64 = 16.0;

/// `ln x` for `x > 0`: the f64's exponent split off exactly and its
/// mantissa's log in f32, so that no `x` is out of f32's range (an
/// infinite, NaN or subnormal `x` by the library's `ln`).
#[inline(always)]
pub(super) fn ln(x: f64) -> f64 {
    let bits = x.to_bits();
    let e = ((bits >> 52) & 0x7ff) as i64;
    if e == 0 || e == 0x7ff {
        return x.ln();
    }
    let mantissa = f64::from_bits((bits & !(0x7ff << 52)) | (1023 << 52));
    f64::from(fast_ln(mantissa as f32)) + (e - 1023) as f64 * std::f64::consts::LN_2
}

/// `ln (1 + e)` for `e >= 0`: below 0.1 its series to the fifth power
/// (within `2e-6` of itself), from there [`ln`].
#[inline(always)]
pub(super) fn ln_1p(e: f64) -> f64 {
    if e < 0.1 {
        e * (1.0 - e * (0.5 - e * (1.0 / 3.0 - e * (0.25 - e * 0.2))))
    } else {
        ln(1.0 + e)
    }
}

/// `ln Gamma(x + n) - ln Gamma(x)` for `x > 0` as the sum of its
/// rising factorial's logs, `ln x + ln (x + 1) + ...`, each by [`ln`]: the
/// urns' terms are defined so, so that a count's step is one [`ln`] and a
/// Gibbs step weighs exactly what the log joint does.
fn ln_rising(x: f64, n: u32) -> f64 {
    (0..n).map(|j| ln(x + f64::from(j))).sum()
}

/// The Dirichlet-multinomial log marginal of counts `row` under
/// `Dirichlet(conc * base)`, through [`ln_rising`].
pub(super) fn ln_urn(row: &[u32], conc: f64, base: &[f64]) -> f64 {
    let total: u32 = row.iter().sum();
    row.iter()
        .zip(base)
        .map(|(&n, &b)| ln_rising(conc * b, n))
        .sum::<f64>()
        - ln_rising(conc, total)
}

/// The change in [`ln_urn`] as one entry, pseudo-count `x` and count `n`,
/// gains a count, the total's change aside (the same for every entry).
#[inline(always)]
pub(super) fn urn_step(x: f64, n: u32) -> f64 {
    ln(x + f64::from(n))
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
