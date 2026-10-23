//! Oracles for `interp`, from closed forms.
//!
//! - Lagrange of degree n through n + 1 nodes reproduces every polynomial of
//!   degree ≤ n, and its derivative reproduces the polynomial's derivative.
//! - Hermite through k nodes with values and derivatives reproduces every
//!   polynomial of degree ≤ 2k − 1, value and derivative.
//! - On a circular orbit r(t) = R (cos ωt, sin ωt, 0) sampled every h, the
//!   Lagrange error at t is below R ωⁿ⁺¹ · |∏(t − tᵢ)| / (n + 1)! over the nodes of
//!   the window centred on t (each component's derivatives are bounded by R ωᵐ),
//!   and Hermite's below R ω²ᵏ · ∏(t − tᵢ)² / (2k)! — the bound of the centred
//!   window, so a window that is not centred (a larger product) fails it.
//!
//! The polynomials have integer coefficients and are evaluated directly (Horner)
//! for the expected values; the tolerance is a few ulps of the largest term.
//!
//! Author: Moroya Sakamoto

// the index pairs the same axis of two 3-vectors
#![allow(clippy::needless_range_loop)]

use alice_space::interp::{self, InterpError};

/// `c[0] + c[1] t + c[2] t² + …` and its derivative.
fn poly(c: &[f64], t: f64) -> (f64, f64) {
    let mut v = 0.0;
    let mut d = 0.0;
    for &ci in c.iter().rev() {
        d = d * t + v;
        v = v * t + ci;
    }
    (v, d)
}

/// Three polynomials of degree `deg` (one per axis), integer coefficients.
fn polys(deg: usize) -> [Vec<f64>; 3] {
    let gen = |seed: i64| -> Vec<f64> {
        (0..=deg)
            .map(|i| {
                let i = i64::try_from(i).unwrap();
                #[allow(clippy::cast_precision_loss)]
                let v = (((seed * 7 + i * 13) % 11) - 5) as f64;
                v
            })
            .collect()
    };
    [gen(1), gen(2), gen(3)]
}

/// Non-uniform node times on [−1.7, 2.9].
fn nodes(count: usize) -> Vec<f64> {
    #[allow(clippy::cast_precision_loss)]
    (0..count)
        .map(|i| -1.7 + 4.6 * (i as f64 / (count - 1) as f64) + 0.01 * ((i * i) % 3) as f64)
        .collect()
}

fn close(a: f64, b: f64, scale: f64) -> bool {
    (a - b).abs() <= 1e-11 * scale.max(1.0)
}

#[test]
fn lagrange_reproduces_polynomials_up_to_its_degree() {
    let mut compared = 0;
    for degree in 1..=10 {
        let t = nodes(degree + 4);
        for deg in 0..=degree {
            let c = polys(deg);
            let values: Vec<[f64; 3]> = t
                .iter()
                .map(|&ti| [poly(&c[0], ti).0, poly(&c[1], ti).0, poly(&c[2], ti).0])
                .collect();
            let last = t[t.len() - 1];
            for q in [t[0], -1.1, 0.0, 0.37, 1.5, 2.3, last] {
                let (p, dp) = interp::lagrange(&t, &values, degree, q).unwrap();
                for k in 0..3 {
                    let (v, d) = poly(&c[k], q);
                    let scale = c[k].iter().map(|x| x.abs()).sum::<f64>() * 3f64.powi(deg as i32);
                    assert!(
                        close(p[k], v, scale),
                        "deg {degree} poly {deg} t={q}: {} vs {v}",
                        p[k]
                    );
                    assert!(
                        close(dp[k], d, 10.0 * scale),
                        "deg {degree} poly {deg} t={q}: d {} vs {d}",
                        dp[k]
                    );
                    compared += 1;
                }
            }
        }
    }
    assert!(compared > 1000, "compared {compared}");
}

#[test]
fn hermite_reproduces_polynomials_up_to_twice_its_nodes_minus_one() {
    let mut compared = 0;
    for nodes_used in 1..=6 {
        let degree = 2 * nodes_used - 1;
        let t = nodes(nodes_used + 3);
        for deg in 0..=degree {
            let c = polys(deg);
            let pos: Vec<[f64; 3]> = t
                .iter()
                .map(|&ti| [poly(&c[0], ti).0, poly(&c[1], ti).0, poly(&c[2], ti).0])
                .collect();
            let vel: Vec<[f64; 3]> = t
                .iter()
                .map(|&ti| [poly(&c[0], ti).1, poly(&c[1], ti).1, poly(&c[2], ti).1])
                .collect();
            let last = t[t.len() - 1];
            for q in [t[0], -0.4, 0.37, 1.5, last] {
                let (p, dp) = interp::hermite(&t, &pos, &vel, degree, q).unwrap();
                for k in 0..3 {
                    let (v, d) = poly(&c[k], q);
                    let scale = c[k].iter().map(|x| x.abs()).sum::<f64>() * 3f64.powi(deg as i32);
                    assert!(
                        close(p[k], v, scale),
                        "Hermite {degree} poly {deg} t={q}: {} vs {v}",
                        p[k]
                    );
                    assert!(
                        close(dp[k], d, 10.0 * scale),
                        "Hermite {degree} poly {deg} t={q}: d {} vs {d}",
                        dp[k]
                    );
                    compared += 1;
                }
            }
        }
    }
    assert!(compared > 500, "compared {compared}");
}

/// The table starts at 10^6 s, far from the origin, so the interpolation has to
/// work relative to its window.
const T0: f64 = 1.0e6;

/// A LEO-like circular orbit: R = 7000 km, period 5800 s, sampled every 60 s.
fn circle(t: f64) -> ([f64; 3], [f64; 3]) {
    let r = 7000.0;
    let w = 2.0 * core::f64::consts::PI / 5800.0;
    // the phase from the table start: w * 1e6 would round the angle to 2e-13 rad
    // (1.6e-9 km at 7000 km), larger than the interpolation error being measured
    let (s, c) = (w * (t - T0)).sin_cos();
    ([r * c, r * s, 0.0], [-r * w * s, r * w * c, 0.0])
}

fn factorial(n: usize) -> f64 {
    #[allow(clippy::cast_precision_loss)]
    (1..=n).map(|k| k as f64).product()
}

#[test]
fn interpolating_a_circular_orbit_stays_within_the_error_bound() {
    let r = 7000.0;
    let w = 2.0 * core::f64::consts::PI / 5800.0;
    let h = 60.0;
    #[allow(clippy::cast_precision_loss)]
    let times: Vec<f64> = (0..60).map(|i| T0 + h * i as f64).collect();
    let pos: Vec<[f64; 3]> = times.iter().map(|&t| circle(t).0).collect();
    let vel: Vec<[f64; 3]> = times.iter().map(|&t| circle(t).1).collect();
    let mut compared = 0;
    for degree in [3usize, 5, 7, 9] {
        let count = degree + 1;
        let k = (degree + 2) / 2;
        for i in 5..50 {
            let q = times[i] + 0.37 * h;
            let exact = circle(q).0;
            // the window CCSDS 502.0 describes: the nodes centred on q (q lies between
            // node i and i + 1, so the window starts (count - 1) / 2 nodes before i)
            let lag_nodes = &times[i - (count - 1) / 2..i - (count - 1) / 2 + count];
            let her_nodes = &times[i - (k - 1) / 2..i - (k - 1) / 2 + k];
            // error bounds of the centred windows at q: |f^(m)| <= R w^m for each axis
            let lag_bound = r * w.powi(count as i32) / factorial(count)
                * lag_nodes.iter().map(|&x| (q - x).abs()).product::<f64>();
            let her_bound = r * w.powi(2 * k as i32) / factorial(2 * k)
                * her_nodes
                    .iter()
                    .map(|&x| (q - x) * (q - x))
                    .product::<f64>();
            let (pl, _) = interp::lagrange(&times, &pos, degree, q).unwrap();
            let (ph, _) = interp::hermite(&times, &pos, &vel, degree, q).unwrap();
            for c in 0..3 {
                // plus the rounding of sums of 7000 km values (a few tens of ulps)
                assert!(
                    (pl[c] - exact[c]).abs() <= lag_bound + 5e-11,
                    "Lagrange {degree} t={q}: {} > {lag_bound}",
                    (pl[c] - exact[c]).abs()
                );
                assert!(
                    (ph[c] - exact[c]).abs() <= her_bound + 5e-11,
                    "Hermite {degree} t={q}: {} > {her_bound}",
                    (ph[c] - exact[c]).abs()
                );
                compared += 2;
            }
        }
    }
    assert!(compared > 1000);
}

#[test]
fn degenerate_tables_are_errors() {
    let t = [0.0, 1.0, 2.0];
    let v = [[0.0; 3]; 3];
    assert_eq!(
        interp::lagrange(&t, &v, 3, 1.0),
        Err(InterpError::NotEnoughSamples(4, 3))
    );
    assert_eq!(
        interp::lagrange(&t, &v, 2, 2.5),
        Err(InterpError::OutOfRange(2.5))
    );
    assert!(matches!(
        interp::lagrange(&t, &v, 2, f64::NAN),
        Err(InterpError::OutOfRange(x)) if x.is_nan()
    ));
    assert_eq!(
        interp::lagrange(&[0.0, 2.0, 1.0], &v, 2, 1.0),
        Err(InterpError::UnorderedTimes(2))
    );
    assert_eq!(
        interp::lagrange(&[0.0, 1.0, 1.0], &v, 2, 1.0),
        Err(InterpError::UnorderedTimes(2))
    );
    assert_eq!(
        interp::lagrange(&[0.0, f64::NAN, 2.0], &v, 2, 1.0),
        Err(InterpError::UnorderedTimes(1))
    );
    assert_eq!(
        interp::lagrange(&t, &v[..2], 1, 1.0),
        Err(InterpError::LengthMismatch)
    );
    assert_eq!(
        interp::hermite(&t, &v, &v[..2], 1, 1.0),
        Err(InterpError::LengthMismatch)
    );
    assert_eq!(
        interp::hermite(&t, &v, &v, 7, 1.0),
        Err(InterpError::NotEnoughSamples(4, 3))
    );
    // degree 0 is a nearest-left lookup of one node: exactly the node's value
    let vals = [[1.0, 2.0, 3.0], [4.0, 5.0, 6.0], [7.0, 8.0, 9.0]];
    assert_eq!(
        interp::lagrange(&t, &vals, 0, 1.4).unwrap().0,
        [4.0, 5.0, 6.0]
    );
    // at a node, any degree returns the node's value
    for d in 0..3 {
        assert_eq!(
            interp::lagrange(&t, &vals, d, 1.0).unwrap().0,
            [4.0, 5.0, 6.0]
        );
    }
}
