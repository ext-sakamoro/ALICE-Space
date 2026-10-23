//! Oracle: atmospheric drag (`drag::atmospheric_drag`).
//!
//! a = −½ ρ Cd (A/m) |v_r| v_r with v_r = v − ω × r, the velocity relative to
//! an atmosphere that rotates with the Earth about z. The expected values are
//! that expression written out here in SI units (m, m/s, kg/m³), independent
//! of the crate's km-based arithmetic; the direction is opposite to v_r; a
//! satellite that moves with the atmosphere feels none.
//!
//! Author: Moroya Sakamoto

use alice_det_math::{sin_cos64, sqrt64};
use alice_space::drag::{atmospheric_drag, DragError, EARTH_ROTATION_RATE};

const EPS: f64 = f64::EPSILON;

fn norm(v: [f64; 3]) -> f64 {
    sqrt64(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
}

/// The drag in m/s² from SI inputs.
fn expected_si(r_m: [f64; 3], v_ms: [f64; 3], rho: f64, cd: f64, am: f64, w: f64) -> [f64; 3] {
    let vr = [v_ms[0] + w * r_m[1], v_ms[1] - w * r_m[0], v_ms[2]];
    let s = norm(vr);
    let k = -0.5 * rho * cd * am * s;
    [k * vr[0], k * vr[1], k * vr[2]]
}

#[test]
fn drag_matches_the_si_expression() {
    let mut compared = 0;
    for k in 0..60 {
        let t = f64::from(k);
        let (s1, c1) = sin_cos64(0.7 * t);
        let (s2, c2) = sin_cos64(0.3 * t);
        let rk = 6_578.0 + 10.0 * t;
        let r = [rk * c1 * c2, rk * s1 * c2, rk * s2];
        let v = [-7.7 * s1, 7.7 * c1 * c2, 0.4 * s2 + 0.1];
        let (rho, cd, am) = (1e-12 * (1.0 + t), 2.2, 0.01 + 0.001 * t);
        let got = atmospheric_drag(r, v, rho, cd, am, EARTH_ROTATION_RATE).unwrap();
        let want_si = expected_si(
            [r[0] * 1e3, r[1] * 1e3, r[2] * 1e3],
            [v[0] * 1e3, v[1] * 1e3, v[2] * 1e3],
            rho,
            cd,
            am,
            EARTH_ROTATION_RATE,
        );
        let want = [want_si[0] / 1e3, want_si[1] / 1e3, want_si[2] / 1e3];
        let tol = 16.0 * EPS * norm(want);
        for (c, (g, w)) in got.iter().zip(&want).enumerate() {
            assert!((g - w).abs() <= tol, "{k} {c}: {g} vs {w}");
        }
        // opposite to the velocity relative to the air
        let vr = [
            v[0] + EARTH_ROTATION_RATE * r[1],
            v[1] - EARTH_ROTATION_RATE * r[0],
            v[2],
        ];
        let cos = (got[0] * vr[0] + got[1] * vr[1] + got[2] * vr[2]) / (norm(got) * norm(vr));
        assert!((cos + 1.0).abs() <= 8.0 * EPS, "{k}: cos {cos}");
        compared += 1;
    }
    assert_eq!(compared, 60);
}

/// A point that co-rotates with the atmosphere (v = ω × r) feels no drag, and
/// a non-rotating atmosphere gives the drag of v itself.
#[test]
fn corotation_is_subtracted() {
    let r = [6_700.0, 1_200.0, 300.0];
    let w = EARTH_ROTATION_RATE;
    let v = [-w * r[1], w * r[0], 0.0];
    assert_eq!(
        atmospheric_drag(r, v, 1e-11, 2.2, 0.02, w).unwrap(),
        [0.0; 3]
    );
    let still = atmospheric_drag(r, [0.0, 7.5, 0.0], 1e-11, 2.2, 0.02, 0.0).unwrap();
    let want = -0.5 * 1e-11 * 2.2 * 0.02 * 7_500.0 * 7_500.0 / 1e3;
    assert!(
        (still[1] - want).abs() <= 8.0 * EPS * want.abs(),
        "{still:?}"
    );
    assert_eq!((still[0], still[2]), (0.0, 0.0));
}

/// Zero density gives zero; a negative or non-finite density, coefficient or
/// area-to-mass ratio, or a non-finite state gives an error.
#[test]
fn degenerate_inputs() {
    let (r, v) = ([7_000.0, 0.0, 0.0], [0.0, 7.5, 0.0]);
    assert_eq!(
        atmospheric_drag(r, v, 0.0, 2.2, 0.02, EARTH_ROTATION_RATE).unwrap(),
        [0.0; 3]
    );
    for (rho, cd, am) in [
        (-1e-12, 2.2, 0.02),
        (1e-12, -2.2, 0.02),
        (1e-12, 2.2, -0.02),
        (f64::NAN, 2.2, 0.02),
        (1e-12, f64::INFINITY, 0.02),
    ] {
        assert_eq!(
            atmospheric_drag(r, v, rho, cd, am, EARTH_ROTATION_RATE),
            Err(DragError::InvalidInput)
        );
    }
    assert_eq!(
        atmospheric_drag(
            [f64::NAN, 0.0, 0.0],
            v,
            1e-12,
            2.2,
            0.02,
            EARTH_ROTATION_RATE
        ),
        Err(DragError::InvalidInput)
    );
}
