//! Oracle: the zonal gravity field (`gravity::ZonalGravity`) and its use by
//! the RK4 propagator.
//!
//! - `tests/golden/zonal_gravity.csv` (scripts/gen_zonal_golden.py): the
//!   potential U = (μ/r)(1 − Σ Jₙ (R/r)ⁿ Pₙ(z/r)) at 60 digits and its
//!   gradient by mpmath's numerical differentiation, so the expected
//!   acceleration shares no formula with the implementation. Three
//!   coefficient sets (J2; J2..J6; J2..J20), 50 positions from the poles to the
//!   equator at 1 to 10 radii.
//! - The closed-form J2 acceleration (e.g. Vallado, *Fundamentals of
//!   Astrodynamics and Applications*, 4th ed., eq. 8-30).
//! - Through the propagator: in an axisymmetric field that does not rotate,
//!   the energy v²/2 − U and the polar angular momentum x·vy − y·vx are
//!   conserved; RK4's drift falls by 2⁴ when the step halves.
//!
//! Author: Moroya Sakamoto

use alice_det_math::{max64, sin_cos64, sqrt64};
use alice_space::gravity::ZonalGravity;
use alice_space::orbit::SpacecraftState;
use alice_space::propagator::propagate_rk4_field;

const EPS: f64 = f64::EPSILON;

struct Row {
    set: String,
    mu: f64,
    radius: f64,
    r: [f64; 3],
    potential: f64,
    a: [f64; 3],
}

fn golden() -> (Vec<(String, Vec<f64>)>, Vec<Row>) {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/golden/zonal_gravity.csv"
    ))
    .unwrap();
    let sets_line = text
        .lines()
        .find_map(|l| l.strip_prefix("# sets: "))
        .expect("coefficient sets");
    let sets = sets_line
        .split("; ")
        .map(|s| {
            let (name, js) = s.split_once('=').unwrap();
            (
                name.to_owned(),
                js.split(' ').map(|j| j.parse().unwrap()).collect(),
            )
        })
        .collect();
    let rows = text
        .lines()
        .filter(|l| !l.starts_with('#') && !l.starts_with("set,"))
        .map(|l| {
            let f: Vec<&str> = l.split(',').collect();
            let v: Vec<f64> = f[1..].iter().map(|x| x.parse().unwrap()).collect();
            Row {
                set: f[0].to_owned(),
                mu: v[0],
                radius: v[1],
                r: [v[2], v[3], v[4]],
                potential: v[5],
                a: [v[6], v[7], v[8]],
            }
        })
        .collect();
    (sets, rows)
}

/// The f64 evaluation rounds the central term in a handful of operations
/// (≤ 8ε relative) and each zonal term, whose Legendre recursion and power of
/// R/r add a few operations per degree, carries at most 8ε per degree of the
/// central magnitude μ/r (|Jₙ (R/r)ⁿ Pₙ| < 1, its derivative < n²).
fn bound(n_zonals: usize, scale: f64) -> f64 {
    let n = n_zonals as f64 + 1.0;
    (16.0 + 8.0 * n * n) * EPS * scale
}

#[test]
fn acceleration_and_potential_match_mpmath() {
    let (sets, rows) = golden();
    let mut compared = 0;
    for row in &rows {
        let js = &sets.iter().find(|(n, _)| *n == row.set).unwrap().1;
        let field = ZonalGravity {
            mu: row.mu,
            radius: row.radius,
            zonals: js,
        };
        let r2 = row.r.iter().map(|c| c * c).sum::<f64>();
        let a = field.acceleration(&row.r);
        let tol_a = bound(js.len(), row.mu / r2);
        for (c, (got, want)) in a.iter().zip(&row.a).enumerate() {
            assert!(
                (got - want).abs() <= tol_a,
                "{} {:?} a{c}: {got} vs {want} (tol {tol_a:e})",
                row.set,
                row.r,
            );
        }
        let u = field.potential(&row.r);
        let tol_u = bound(js.len(), row.potential.abs());
        assert!(
            (u - row.potential).abs() <= tol_u,
            "{} {:?} U: {u} vs {}",
            row.set,
            row.r,
            row.potential
        );
        compared += 1;
    }
    assert_eq!(compared, 150);
}

/// a = −μ r/r³ + the J2 term: aₓ = −μx/r³ · (3/2) J2 (R/r)² (1 − 5z²/r²)
/// (likewise y), a_z = −μz/r³ · (3/2) J2 (R/r)² (3 − 5z²/r²).
#[test]
fn j2_matches_the_closed_form() {
    let (mu, re, j2) = (398_600.441_5, 6_378.136_3, 1.082_63e-3);
    let field = ZonalGravity {
        mu,
        radius: re,
        zonals: &[j2],
    };
    let mut compared = 0;
    for k in 0..200 {
        let t = f64::from(k);
        let (s1, c1) = sin_cos64(0.37 * t);
        let (s2, c2) = sin_cos64(0.11 * t);
        let r = [7_000.0 * c1, 7_000.0 * s1 * c2, 7_500.0 * s2];
        let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        let rn = sqrt64(r2);
        let (z2, rho2) = (r[2] * r[2] / r2, re * re / r2);
        let k2 = 1.5 * j2 * rho2;
        let want = [
            -mu * r[0] / (rn * r2) * (1.0 + k2 * (1.0 - 5.0 * z2)),
            -mu * r[1] / (rn * r2) * (1.0 + k2 * (1.0 - 5.0 * z2)),
            -mu * r[2] / (rn * r2) * (1.0 + k2 * (3.0 - 5.0 * z2)),
        ];
        let got = field.acceleration(&r);
        for c in 0..3 {
            assert!(
                (got[c] - want[c]).abs() <= 32.0 * EPS * mu / r2,
                "{r:?} {c}: {} vs {}",
                got[c],
                want[c]
            );
        }
        compared += 1;
    }
    assert_eq!(compared, 200);
}

fn energy(field: &ZonalGravity, s: &SpacecraftState) -> f64 {
    let v = s.velocity_km_s;
    (v[0] * v[0] + v[1] * v[1] + v[2] * v[2]) / 2.0 - field.potential(&s.position_km)
}

fn hz(s: &SpacecraftState) -> f64 {
    s.position_km[0] * s.velocity_km_s[1] - s.position_km[1] * s.velocity_km_s[0]
}

/// A 7000 km orbit inclined 51.6°, slightly eccentric, for 3 hours in the
/// J2..J6 field: the energy and the polar angular momentum stay constant to
/// RK4's truncation, which falls by about 2⁴ = 16 when the step halves.
#[test]
fn the_propagator_conserves_energy_and_polar_momentum_in_the_field() {
    let js = [1.082_63e-3, -2.53e-6, -1.62e-6, -2.3e-7, 5.4e-7];
    let field = ZonalGravity {
        mu: 398_600.441_5,
        radius: 6_378.136_3,
        zonals: &js,
    };
    let (si, ci) = sin_cos64(51.6f64.to_radians());
    let start = SpacecraftState {
        position_km: [7_000.0, 0.0, 0.0],
        velocity_km_s: [0.0, 7.6 * ci, 7.6 * si],
        timestamp_ns: 0,
        fuel_kg: 0.0,
    };
    let drift = |dt: f64| {
        let steps = (10_800.0 / dt) as usize;
        let path = propagate_rk4_field(&start, &field, dt, steps);
        assert_eq!(path.len(), steps + 1);
        let (e0, h0) = (energy(&field, &start), hz(&start));
        path.iter().fold((0.0f64, 0.0f64), |(de, dh), s| {
            (
                max64(de, (energy(&field, s) - e0).abs() / e0.abs()),
                max64(dh, (hz(s) - h0).abs() / h0.abs()),
            )
        })
    };
    let (e20, h20) = drift(20.0);
    let (e10, h10) = drift(10.0);
    assert!(e10 < 1e-10 && h10 < 1e-10, "{e10:e} {h10:e}");
    // fourth order: halving the step divides the drift by 2⁴ (within a
    // factor of 2 for the higher-order terms)
    let ratio = e20 / e10;
    assert!((8.0..32.0).contains(&ratio), "energy drift ratio {ratio}");
    assert!(h20 > 0.0 && e20 > 0.0);
    // a propagation that ignored the zonal terms would drift by the size of
    // J2 itself (10⁻³ relative); this field's energy is conserved far below
    assert!(e20 < 1e-8, "{e20:e}");
    let _ = h20;
}

/// With no zonal coefficient the field is the point mass of `TwoBodyAccel`
/// (within the rounding of the two orders of operations, 4 ulp); at the
/// centre the field is singular and both quantities are NaN, not a finite
/// value; a NaN position gives NaN and an infinite distance a zero potential.
#[test]
fn degenerate_fields_and_positions() {
    use alice_space::propagator::TwoBodyAccel;
    let mu = 398_600.441_5;
    let point = ZonalGravity {
        mu,
        radius: 6_378.136_3,
        zonals: &[],
    };
    let two_body = TwoBodyAccel::new(mu);
    for r in [
        [7_000.0, -300.0, 20.0],
        [0.0, 0.0, 42_164.0],
        [-1.0, 2.0, -3.0],
    ] {
        let (a, b) = (point.acceleration(&r), two_body.acceleration(&r));
        for (x, y) in a.iter().zip(&b) {
            assert!(
                (x - y).abs() <= 4.0 * EPS * y.abs(),
                "{r:?}: {a:?} vs {b:?}"
            );
        }
        let rn = sqrt64(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
        assert!((point.potential(&r) - mu / rn).abs() <= 2.0 * EPS * mu / rn);
    }
    let earth = ZonalGravity {
        mu,
        radius: 6_378.136_3,
        zonals: &[1.082_63e-3, -2.53e-6],
    };
    assert!(earth.acceleration(&[0.0; 3]).iter().all(|v| v.is_nan()));
    assert!(earth.potential(&[0.0; 3]).is_nan());
    assert!(earth
        .acceleration(&[f64::NAN, 1.0, 1.0])
        .iter()
        .all(|v| v.is_nan()));
    // at an infinite distance ρ = 0, the sum is 1 and μ/∞ = 0
    assert_eq!(earth.potential(&[f64::INFINITY, 0.0, 0.0]), 0.0);
}

/// The J2 regression of the node: over 30 orbits of a = 7000 km, e = 0.01,
/// i = 50° in the J2 field, the right ascension of the ascending node
/// (from h = r × v) drifts at the first-order secular rate
/// Ω̇ = −(3/2) n J2 (R/p)² cos i (e.g. Vallado, 4th ed., eq. 9-41). First-order
/// theory leaves out the short-periodic difference between osculating and
/// mean elements (of order J2 (R/a)² ≈ 9 × 10⁻⁴ rad, 0.5 % of the 0.16 rad
/// drift) and the O(J2) change of the mean motion: within 1 %.
#[test]
fn j2_regresses_the_node_at_the_secular_rate() {
    let (mu, re, j2) = (398_600.441_5, 6_378.136_3, 1.082_63e-3);
    let field = ZonalGravity {
        mu,
        radius: re,
        zonals: &[j2],
    };
    let (a, e, inc) = (7_000.0, 0.01, 50f64.to_radians());
    // perigee on the x axis, ascending node there too (Ω = 0, ω = 0)
    let rp = a * (1.0 - e);
    let vp = sqrt64(mu * (1.0 + e) / rp);
    let (si, ci) = sin_cos64(inc);
    let start = SpacecraftState {
        position_km: [rp, 0.0, 0.0],
        velocity_km_s: [0.0, vp * ci, vp * si],
        timestamp_ns: 0,
        fuel_kg: 0.0,
    };
    let n = sqrt64(mu / (a * a * a));
    let period = 2.0 * core::f64::consts::PI / n;
    let dt = 10.0;
    let steps = (30.0 * period / dt) as usize;
    let end = propagate_rk4_field(&start, &field, dt, steps)
        .pop()
        .unwrap();
    let (r, v) = (end.position_km, end.velocity_km_s);
    let h = [
        r[1] * v[2] - r[2] * v[1],
        r[2] * v[0] - r[0] * v[2],
        r[0] * v[1] - r[1] * v[0],
    ];
    let node = alice_det_math::atan2_64(h[0], -h[1]);
    let t = steps as f64 * dt;
    let p = a * (1.0 - e * e);
    let rate = -1.5 * n * j2 * (re / p) * (re / p) * ci;
    let want = rate * t;
    assert!(want < -0.15, "the drift is large enough to measure: {want}");
    assert!(
        ((node - want) / want).abs() <= 0.01,
        "node {node} rad against the secular {want} rad"
    );
}
