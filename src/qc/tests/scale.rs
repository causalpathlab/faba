//! One scale per quantity: its defaults and `--scale`, where values sit
//! on its ruler, its ends and its bins.

use super::*;

/// The defaults read well; `--scale` sets each quantity named, refuses a
/// log of the signed log odds and names it cannot read.
#[test]
fn scales_default_and_are_set_by_name() {
    let d = Scales::default();
    assert_eq!(d.of(Quantity::Pvalue), Scale::Linear);
    assert_eq!(d.of(Quantity::Fold), Scale::Log);
    for q in [Quantity::Coverage, Quantity::Converted, Quantity::Cells] {
        assert_eq!(d.of(q), Scale::Log, "{q:?}");
    }
    for q in [Quantity::EditRatio, Quantity::LogOdds] {
        assert_eq!(d.of(q), Scale::Linear, "{q:?}");
    }
    let s = Scales::parse(&["pv=log", "fold=linear", "edit-ratio=log"]).unwrap();
    assert_eq!(s.of(Quantity::Pvalue), Scale::Log);
    assert_eq!(s.of(Quantity::Fold), Scale::Linear);
    assert_eq!(s.of(Quantity::EditRatio), Scale::Log);
    assert_eq!(s.of(Quantity::Coverage), Scale::Log, "left as it was");
    let e = Scales::parse(&["log_odds=log"]).unwrap_err().to_string();
    assert!(e.contains("log_odds") && e.contains("linear"), "{e}");
    for bad in ["bogus=log", "pv=sqrt", "pv"] {
        assert!(Scales::parse(&[bad]).is_err(), "{bad}");
    }
}

/// A value's place on its line is what the ruler's axis reads it at:
/// the axis' own fraction, on a log and a linear ruler alike.
#[test]
fn a_ruler_places_values_as_its_axis_does() {
    let shown: Vec<f32> = (0..500).map(|i| (i % 97) as f32 * 1.7).collect();
    for scale in [Scale::Log, Scale::Linear] {
        let r = Ruler::new(Quantity::Fold, scale, &shown, &[], 40);
        let (a, b) = r.span();
        for v in [0.0, 1.0, 3.3, 50.0, 160.0] {
            let along = (r.along(v) - a) / (b - a);
            assert!((along - r.axis.frac(v)).abs() < 1e-12, "{scale:?} {v}");
        }
        assert_eq!(r.axis.log, scale == Scale::Log);
        let ticks = r.ticks();
        assert!(!ticks.is_empty());
        for (at, _) in ticks {
            assert!(at >= a - 1e-9 && at <= b + 1e-9, "{at} in {a}..{b}");
        }
    }
}

/// An infinity is clamped to its end and counted in the end bin: none is
/// lost. A p-value of 0 is shown as an infinite `-log10 p`.
#[test]
fn infinities_land_at_the_ends() {
    let shown = [0.5, 1.0, 2.0, f32::INFINITY, 3.0];
    let r = Ruler::new(Quantity::Fold, Scale::Log, &shown, &[], 10);
    let c = r.clamped(&shown);
    assert_eq!(c[3], r.axis.hi as f32);
    let counts = r.axis.counts(&c, r.bins, |_| true);
    assert_eq!(counts.iter().sum::<usize>(), shown.len());
    assert!(counts[r.bins - 1] >= 1);
    assert_eq!(Quantity::Pvalue.shown(0.0), f64::INFINITY);
    assert!((Quantity::Pvalue.shown(0.01) - 2.0).abs() < 1e-12);
    let lo = Ruler::new(
        Quantity::LogOdds,
        Scale::Linear,
        &[f32::NEG_INFINITY, 1.0],
        &[],
        5,
    );
    assert_eq!(lo.clamped(&[f32::NEG_INFINITY])[0], lo.axis.lo as f32);
}

/// A count on a linear scale gets whole-number bins, each centred on its
/// count; edges given are inside the range.
#[test]
fn counts_on_a_linear_scale_bin_whole_numbers() {
    let shown: Vec<f32> = (0..=37).map(|k| k as f32).collect();
    let r = Ruler::new(Quantity::Converted, Scale::Linear, &shown, &[50.0], 10);
    assert_eq!(r.axis.lo, -0.5);
    let width = (r.axis.hi - r.axis.lo) / r.bins as f64;
    assert_eq!(width, width.round(), "{width}");
    assert!(r.axis.hi >= 50.0);
    for k in 0..=50 {
        let k = k as f64;
        assert_eq!(
            r.axis.bin(k - 0.49, r.bins),
            r.axis.bin(k + 0.49, r.bins),
            "{k}"
        );
    }
}
