//! Oracle: the conical shadow of an occulting body, from positions.
//!
//! - `tests/golden/eclipse.csv` (scripts/gen_eclipse_golden.py): satellites
//!   behind the Earth (LEO to 1.5 million km) and behind the Moon, with the
//!   kind of eclipse and the visible fraction of the Sun's disk computed by
//!   mpmath at 40 digits, the disk overlap integrated numerically (not from
//!   the circle-intersection formula used here).
//! - Closed-form geometry: on the Sun–Earth axis the umbra ends where the
//!   Earth and the Sun have the same apparent radius,
//!   `L = R_E · d / (R_S − R_E)` beyond the Earth (≈ 1.384 × 10⁶ km); just
//!   inside is total, just outside annular. Off the axis the edge of the
//!   umbra and of the penumbra are where the separation equals `b − a` and
//!   `a + b`.
//!
//! Author: Moroya Sakamoto

mod common_eclipse {
    pub fn rows() -> Vec<Vec<String>> {
        let path = format!("{}/tests/golden/eclipse.csv", env!("CARGO_MANIFEST_DIR"));
        let text = std::fs::read_to_string(&path).unwrap();
        let out: Vec<Vec<String>> = text
            .lines()
            .filter(|l| !l.starts_with('#'))
            .skip(1)
            .map(|l| l.split(',').map(str::to_owned).collect())
            .collect();
        assert!(!out.is_empty(), "{path}: no rows");
        out
    }
}

use alice_det_math::asin64;
use alice_space::eclipse::{
    conical_illumination, earth_illumination, Eclipse, EclipseError, EARTH_EQUATORIAL_RADIUS_KM,
    SUN_RADIUS_KM,
};

fn f(s: &str) -> f64 {
    s.parse().unwrap()
}

#[test]
fn illumination_matches_the_integrated_overlap() {
    let mut compared = 0;
    let mut kinds = [0usize; 4];
    for r in common_eclipse::rows() {
        let v = |i: usize| f(&r[i]);
        let sat = [v(1), v(2), v(3)];
        let sun = [v(4), v(5), v(6)];
        let body = [v(7), v(8), v(9)];
        let got = conical_illumination(sat, sun, body, v(10)).unwrap();
        let want = match r[11].as_str() {
            "none" => Eclipse::None,
            "partial" => Eclipse::Partial,
            "annular" => Eclipse::Annular,
            "total" => Eclipse::Total,
            k => panic!("{k}"),
        };
        let ctx = format!("{} {sat:?}", r[0]);
        assert_eq!(got.eclipse, want, "{ctx}");
        kinds[want as usize] += 1;
        // angles: asin / atan2 of f64 inputs, a few ulp
        for (g, w, name) in [
            (got.sun_angular_radius, v(13), "a"),
            (got.body_angular_radius, v(14), "b"),
            (got.separation, v(15), "c"),
        ] {
            assert!(
                (g - w).abs() <= 1e-15 + 4.0 * f64::EPSILON * w,
                "{ctx} {name}: {g:e} vs {w:e}"
            );
        }
        // the fraction: the closed form's acos lose digits as the disks
        // approach tangency (d acos(x)/dx = 1/√(1−x²)); 1e-9 holds over the table
        assert!(
            (got.fraction - v(12)).abs() <= 1e-9,
            "{ctx}: {} vs {}",
            got.fraction,
            v(12)
        );
        compared += 1;
    }
    assert_eq!(compared, 475);
    assert!(kinds.iter().all(|&k| k > 0), "{kinds:?}");
}

#[test]
fn the_umbra_ends_where_the_earth_and_the_sun_look_the_same_size() {
    let d = 149_597_870.7;
    let sun = [d, 0.0, 0.0];
    let l = EARTH_EQUATORIAL_RADIUS_KM * d / (SUN_RADIUS_KM - EARTH_EQUATORIAL_RADIUS_KM);
    assert!((l - 1.384e6).abs() < 1e3, "{l}");
    let inside = earth_illumination([-(l * (1.0 - 1e-9)), 0.0, 0.0], sun).unwrap();
    let outside = earth_illumination([-(l * (1.0 + 1e-9)), 0.0, 0.0], sun).unwrap();
    assert_eq!(inside.eclipse, Eclipse::Total);
    assert_eq!(inside.fraction, 0.0);
    assert_eq!(outside.eclipse, Eclipse::Annular);
    // annular: the visible ring is 1 − b²/a², tiny just past the tip
    assert!(
        outside.fraction > 0.0 && outside.fraction < 1e-8,
        "{}",
        outside.fraction
    );
    // on the far side of the umbra tip, at 2L, the ring is 1 − (b/a)²
    let far = earth_illumination([-2.0 * l, 0.0, 0.0], sun).unwrap();
    let a = asin64(SUN_RADIUS_KM / (d + 2.0 * l));
    let b = asin64(EARTH_EQUATORIAL_RADIUS_KM / (2.0 * l));
    assert!(
        (far.fraction - (1.0 - (b / a) * (b / a))).abs() < 1e-12,
        "{}",
        far.fraction
    );
}

#[test]
fn the_shadow_edges_are_where_the_separation_is_b_minus_a_and_a_plus_b() {
    let d = 149_597_870.7;
    let sun = [d, 0.0, 0.0];
    let x = -7000.0;
    let a = asin64(SUN_RADIUS_KM / (d - x));
    // offset y such that the separation c equals `target`: solve with the
    // same geometry (the satellite at (x, y), body at the origin)
    let edge = |target: f64| {
        let (mut lo, mut hi) = (0.0, 3.0 * EARTH_EQUATORIAL_RADIUS_KM);
        for _ in 0..200 {
            let y = 0.5 * (lo + hi);
            let i = conical_illumination([x, y, 0.0], sun, [0.0; 3], EARTH_EQUATORIAL_RADIUS_KM)
                .unwrap();
            if i.separation - (i.body_angular_radius + target) < 0.0 {
                lo = y;
            } else {
                hi = y;
            }
        }
        lo
    };
    let y_umbra = edge(-a);
    let y_penumbra = edge(a);
    let at = |y: f64| {
        conical_illumination([x, y, 0.0], sun, [0.0; 3], EARTH_EQUATORIAL_RADIUS_KM).unwrap()
    };
    assert_eq!(at(y_umbra * (1.0 - 1e-9)).eclipse, Eclipse::Total);
    assert_eq!(at(y_umbra * (1.0 + 1e-9)).eclipse, Eclipse::Partial);
    assert_eq!(at(y_penumbra * (1.0 - 1e-9)).eclipse, Eclipse::Partial);
    assert_eq!(at(y_penumbra * (1.0 + 1e-9)).eclipse, Eclipse::None);
    // continuity of the fraction across both edges
    assert!(at(y_umbra * (1.0 + 1e-9)).fraction < 1e-6);
    assert!(at(y_penumbra * (1.0 - 1e-9)).fraction > 1.0 - 1e-6);
    // the fraction grows monotonically from the umbra to the penumbra edge
    let mut prev = -1.0;
    for k in 0..=100 {
        let y = y_umbra + (y_penumbra - y_umbra) * f64::from(k) / 100.0;
        let fr = at(y).fraction;
        assert!(fr >= prev, "k={k}: {fr} < {prev}");
        prev = fr;
    }
}

#[test]
fn degenerate_input_is_an_error() {
    let sun = [1.5e8, 0.0, 0.0];
    assert_eq!(
        earth_illumination([1000.0, 0.0, 0.0], sun).unwrap_err(),
        EclipseError::InsideBody
    );
    assert_eq!(
        earth_illumination([0.0, 0.0, 0.0], sun).unwrap_err(),
        EclipseError::InsideBody
    );
    assert_eq!(
        earth_illumination([1.5e8, 0.0, 0.0], sun).unwrap_err(),
        EclipseError::InsideSun
    );
    assert_eq!(
        earth_illumination([f64::NAN, 0.0, 7000.0], sun).unwrap_err(),
        EclipseError::NonFinite
    );
    assert_eq!(
        earth_illumination([7000.0, 0.0, 0.0], [f64::INFINITY, 0.0, 0.0]).unwrap_err(),
        EclipseError::NonFinite
    );
    assert_eq!(
        conical_illumination([7000.0, 0.0, 0.0], sun, [0.0; 3], -1.0).unwrap_err(),
        EclipseError::NonFinite
    );
    // the satellite exactly between the Earth and the Sun is sunlit
    let i = earth_illumination([7000.0, 0.0, 0.0], sun).unwrap();
    assert_eq!((i.eclipse, i.fraction), (Eclipse::None, 1.0));
}

/// First-order bound on the error of the visible fraction for one row,
/// derived from the arithmetic of `conical_illumination` (ε = 2⁻⁵²):
///
/// - inputs: `u = sun − sat`, `v = body − sat` are rounded (ε/2 per
///   component, an angle of at most ε/2 for each direction); `|u×v|` is
///   within √3·ε·|u||v| and `u·v` within 2ε·|u||v| of their exact values, so
///   `c = atan2(|u×v|, u·v)` is within (|δS| + |δD|)/(|u||v|) plus the
///   directions and atan2's own rounding: δc ≤ (8 + 2c)ε, in any frame;
///   `a = asin(R/|u|)` with `|u|` within 2ε relative: δa ≤ 3ε·tan a + 2ε·a,
///   and δb likewise;
/// - propagation through F = 1 − A/(πa²) (or the visible segments, the same
///   derivatives): ∂A/∂c = −2y, ∂A/∂b = 2θ_b·b, ∂A/∂a = 2θ_a·a, and
///   ∂F/∂a has the extra 2A/(πa³) of the denominator;
/// - evaluation: x_a = (c² + a² − b²)/(2c) within 3ε·(c² + a² + b²)/(2c)
///   plus ε|x_a|; θ = acos(x/a) moves the
///   segment by 2·sin²θ·δθ = 2·sinθ·δx/a + 4ε·θ·sin²θ; the segment form
///   (u − sin u)/2 is within 4ε·θ; the sum and the division ε·A, and the
///   final subtraction from 1 (or the visible fraction) ε.
///
/// The c term dominates near the edges of a large body: an earlier bound that
/// took δc relative to c (and omitted 2A/(πa³)) was exceeded by 0.6 % on a
/// rotated row; the cancellation in x_a is ~8 times smaller than it there.
fn fraction_bound(a: f64, b: f64, c: f64, kind: Eclipse) -> f64 {
    use alice_det_math::{acos64, sin64, sqrt64};
    let e = f64::EPSILON;
    let tan = |x: f64| sin64(x) / alice_det_math::cos64(x);
    let da = 3.0 * e * tan(a) + 2.0 * e * a;
    let db = 3.0 * e * tan(b) + 2.0 * e * b;
    let dc = (8.0 + 2.0 * c) * e;
    let pi = core::f64::consts::PI;
    match kind {
        Eclipse::None | Eclipse::Total => 0.0,
        Eclipse::Annular => (2.0 * b * b / (a * a * a)) * da + (2.0 * b / (a * a)) * db + 3.0 * e,
        Eclipse::Partial => {
            let xa = (c * c + a * a - b * b) / (2.0 * c);
            let xb = (c * c - a * a + b * b) / (2.0 * c);
            let (ta, tb) = (acos64(xa / a), acos64(xb / b));
            let seg = |t: f64| (2.0 * t - sin64(2.0 * t)) / 2.0;
            let area = a * a * seg(ta) + b * b * seg(tb);
            let y = sqrt64(a * a - xa * xa);
            let prop = (2.0 * y * dc + 2.0 * tb * b * db + (2.0 * ta * a + 2.0 * area / a) * da)
                / (pi * a * a);
            let dxa = 3.0 * e * (c * c + a * a + b * b) / (2.0 * c) + e * xa.abs();
            let dxb = 3.0 * e * (c * c + a * a + b * b) / (2.0 * c) + e * xb.abs();
            let (sa, sb) = (sin64(ta), sin64(tb));
            let eval_a = a * a * (2.0 * sa * dxa / a + 4.0 * e * ta * sa * sa + 4.0 * e * ta);
            let eval_b = b * b * (2.0 * sb * dxb / b + 4.0 * e * tb * sb * sb + 4.0 * e * tb);
            let eval = (eval_a + eval_b + e * area) / (pi * a * a) + e;
            prop + eval
        }
    }
}

/// Grazing configurations (`tests/golden/eclipse_grazing.csv`): satellites
/// 1e-13 to 1e-1 (relative offset) from every edge with the Sun on +x, and
/// every Earth configuration again under two seeded random rotations of the
/// scene, against the circle-intersection formula at 50 digits. The kind must
/// match, and the visible fraction must be inside [0, 1] with no clamping and
/// within [`fraction_bound`] of the exact value.
#[test]
fn grazing_configurations_stay_within_the_derived_bound() {
    let path = format!(
        "{}/tests/golden/eclipse_grazing.csv",
        env!("CARGO_MANIFEST_DIR")
    );
    let text = std::fs::read_to_string(&path).unwrap();
    let mut compared = 0;
    let mut worst_ratio: f64 = 0.0;
    for line in text.lines().filter(|l| !l.starts_with('#')).skip(1) {
        let r: Vec<&str> = line.split(',').collect();
        let v = |i: usize| f(r[i]);
        let got = conical_illumination(
            [v(1), v(2), v(3)],
            [v(4), v(5), v(6)],
            [v(7), v(8), v(9)],
            v(10),
        )
        .unwrap();
        let want = match r[11] {
            "none" => Eclipse::None,
            "partial" => Eclipse::Partial,
            "annular" => Eclipse::Annular,
            "total" => Eclipse::Total,
            k => panic!("{k}"),
        };
        assert_eq!(got.eclipse, want, "{line}");
        let err = (got.fraction - v(12)).abs();
        let bound = fraction_bound(
            got.sun_angular_radius,
            got.body_angular_radius,
            got.separation,
            want,
        );
        assert!(
            err <= bound,
            "{line}: {} off by {err:e} > {bound:e}",
            got.fraction
        );
        assert!(
            (0.0..=1.0).contains(&got.fraction),
            "{line}: {}",
            got.fraction
        );
        if bound > 0.0 {
            worst_ratio = alice_det_math::max64(worst_ratio, err / bound);
        }
        compared += 1;
    }
    assert_eq!(compared, 5600);
    assert!(worst_ratio > 0.0, "the comparison saw no rounding at all");
}
