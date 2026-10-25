//! Oracle: the time-dependent force model (`forces::ForceModel`) and its RK4
//! propagation.
//!
//! `tests/golden/forces.csv` (scripts/gen_forces_golden.py) holds the same
//! propagations computed by an independent reference: ERFA through pyerfa for
//! the Sun, the Moon and the GCRS → CIRS matrix, the forces written from their
//! definitions (third body as the direct difference, the conical shadow of
//! Montenbruck & Gill), the same RK4 steps and time arithmetic. The two differ
//! only by rounding: of ERFA's libm against alice-det-math, of the direct
//! difference against Battin's form, and of the order of the operations.
//!
//! Author: Moroya Sakamoto

use alice_det_math::{exp64, max64, sin_cos64, sqrt64};
use alice_ephemeris::Jd2;
use alice_space::forces::{DragParams, ForceModel, SrpParams};
use alice_space::gravity::ZonalGravity;
use alice_space::orbit::SpacecraftState;

const J: [f64; 5] = [1.082_626_68e-3, -2.5327e-6, -1.6196e-6, -2.273e-7, 5.407e-7];

fn density(_: Jd2, r: [f64; 3]) -> f64 {
    let h = sqrt64(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]) - 6_378.137;
    3.725e-12 * exp64(-(h - 400.0) / 58.515)
}

struct Case {
    name: &'static str,
    epoch: Jd2,
    r: [f64; 3],
    v: [f64; 3],
    dt: f64,
    steps: usize,
    srp: Option<SrpParams>,
    drag: bool,
}

fn cases() -> [Case; 2] {
    let (si, ci) = sin_cos64(51.6f64.to_radians());
    [
        Case {
            name: "leo",
            epoch: Jd2 {
                jd1: 2_461_100.5,
                jd2: 0.25,
            },
            r: [6_778.137, 0.0, 0.0],
            v: [0.0, 7.668_6 * ci, 7.668_6 * si],
            dt: 20.0,
            steps: 300,
            srp: Some(SrpParams {
                cr: 1.3,
                area_to_mass: 0.02,
            }),
            drag: true,
        },
        Case {
            name: "geo",
            epoch: Jd2 {
                jd1: 2_461_100.5,
                jd2: 0.0,
            },
            r: [42_164.0, 0.0, 0.0],
            v: [0.0, 3.074_6, 0.0],
            dt: 600.0,
            steps: 144,
            srp: Some(SrpParams {
                cr: 1.5,
                area_to_mass: 0.01,
            }),
            drag: false,
        },
    ]
}

fn model<'a>(c: &Case, dens: &'a dyn Fn(Jd2, [f64; 3]) -> f64) -> ForceModel<'a> {
    ForceModel {
        epoch: c.epoch,
        gravity: ZonalGravity {
            mu: 398_600.441_5,
            radius: 6_378.136_3,
            zonals: &J,
        },
        sun_mu: Some(132_712_440_041.279_42),
        moon_mu: Some(4_902.800_118),
        srp: c.srp,
        drag: c.drag.then_some(DragParams {
            cd: 2.2,
            area_to_mass: 0.02,
            density: dens,
        }),
    }
}

fn golden() -> Vec<(String, usize, [f64; 6])> {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/golden/forces.csv"
    ))
    .unwrap();
    text.lines()
        .filter(|l| l.starts_with("leo,") || l.starts_with("geo,"))
        .map(|l| {
            let f: Vec<&str> = l.split(',').collect();
            let v: Vec<f64> = f[2..].iter().map(|x| x.parse().unwrap()).collect();
            (
                f[0].to_owned(),
                f[1].parse().unwrap(),
                core::array::from_fn(|i| v[i]),
            )
        })
        .collect()
}

/// Each RK4 stage rounds the state at a few ulp of |r| (two per stage here),
/// so N steps accumulate at most N · 4 · 2 ulp(|r|); along-track errors then
/// grow by about 3 · 2π · t/T over an arc of t against the period T (the
/// drift of the mean motion). The velocity bound is the position bound times
/// the mean motion 2π/T.
fn tolerance(c: &Case) -> (f64, f64) {
    let r = sqrt64(c.r[0] * c.r[0] + c.r[1] * c.r[1] + c.r[2] * c.r[2]);
    let period = 2.0 * core::f64::consts::PI * sqrt64(r * r * r / 398_600.441_5);
    let arc = c.dt * c.steps as f64;
    let growth = 3.0 * 2.0 * core::f64::consts::PI * (arc / period + 1.0);
    let p = c.steps as f64 * 4.0 * 2.0 * f64::EPSILON * r * growth;
    (p, p * 2.0 * core::f64::consts::PI / period)
}

#[test]
fn propagation_matches_the_independent_reference() {
    let rows = golden();
    let mut compared = 0;
    for c in cases() {
        let m = model(&c, &density);
        let start = SpacecraftState {
            position_km: c.r,
            velocity_km_s: c.v,
            timestamp_ns: 0,
            fuel_kg: 0.0,
        };
        let path = m.propagate(&start, c.dt, c.steps).unwrap();
        let (mut worst_p, mut worst_v): (f64, f64) = (0.0, 0.0);
        for (_, k, want) in rows.iter().filter(|(n, _, _)| n == c.name) {
            let s = &path[*k];
            let dp = sqrt64(
                (0..3)
                    .map(|i| (s.position_km[i] - want[i]) * (s.position_km[i] - want[i]))
                    .sum(),
            );
            let dv = sqrt64(
                (0..3)
                    .map(|i| {
                        (s.velocity_km_s[i] - want[3 + i]) * (s.velocity_km_s[i] - want[3 + i])
                    })
                    .sum(),
            );
            worst_p = max64(worst_p, dp);
            worst_v = max64(worst_v, dv);
            compared += 1;
        }
        let (tol_p, tol_v) = tolerance(&c);
        assert!(
            worst_p <= tol_p,
            "{}: position {worst_p:e} km > {tol_p:e}",
            c.name
        );
        assert!(
            worst_v <= tol_v,
            "{}: velocity {worst_v:e} km/s > {tol_v:e}",
            c.name
        );
        // the comparison saw the rounding of two different implementations
        assert!(worst_p > 0.0, "{}: identical to the last bit", c.name);
    }
    assert_eq!(compared, 47);
}

/// A model at an epoch outside the span of an ephemeris it uses is refused, not
/// answered with an extrapolated position: the Sun (`epv00`, 1900–2100) when the
/// Sun's perturbation or the radiation pressure is on, the Moon (`moon98`,
/// 1950–2100) when the Moon's perturbation or the radiation pressure is on. A
/// model that uses neither body takes any epoch.
#[test]
fn an_epoch_outside_an_ephemeris_span_is_refused() {
    use alice_space::forces::{Body, ForceError};
    let dens = density;
    let r = [7_000.0, 0.0, 0.0];
    let v = [0.0, 7.5, 0.0];
    let at = |jd1: f64, sun: bool, moon: bool, srp: bool| {
        let c = &cases()[0];
        let mut m = model(c, &dens);
        m.epoch = Jd2 { jd1, jd2: 0.0 };
        m.sun_mu = if sun { m.sun_mu } else { None };
        m.moon_mu = if moon { m.moon_mu } else { None };
        m.srp = if srp {
            Some(SrpParams {
                cr: 1.5,
                area_to_mass: 0.01,
            })
        } else {
            None
        };
        m.drag = None;
        m.acceleration(0.0, r, v)
    };
    // 2026: inside both spans
    assert!(at(2_461_100.5, true, true, true).is_ok());
    // 2300: outside the Sun's span (and the Moon's)
    assert_eq!(
        at(2_561_118.5, true, false, false),
        Err(ForceError::OutOfEphemerisRange(Body::Sun))
    );
    assert_eq!(
        at(2_561_118.5, false, true, false),
        Err(ForceError::OutOfEphemerisRange(Body::Moon))
    );
    assert_eq!(
        at(2_561_118.5, false, false, true),
        Err(ForceError::OutOfEphemerisRange(Body::Sun))
    );
    // 1925: inside the Sun's span, outside the Moon's
    assert!(at(2_424_151.5, true, false, false).is_ok());
    assert_eq!(
        at(2_424_151.5, true, true, false),
        Err(ForceError::OutOfEphemerisRange(Body::Moon))
    );
    // 2075: inside both
    assert!(at(2_478_938.5, true, true, true).is_ok());
    // neither body: no ephemeris to be outside of
    assert!(at(2_561_118.5, false, false, false).is_ok());
    // a non-finite epoch is outside every span
    assert_eq!(
        at(f64::NAN, true, false, false),
        Err(ForceError::OutOfEphemerisRange(Body::Sun))
    );
}
