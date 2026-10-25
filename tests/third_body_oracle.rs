//! Oracle: the third-body perturbation (`gravity::third_body`).
//!
//! `tests/golden/third_body.csv` (scripts/gen_third_body_golden.py): the direct
//! difference a = −μ₃ ((r − s)/|r − s|³ + s/|s|³) at 60 digits for satellites
//! from LEO to 300,000 km against Sun- and Moon-like bodies. In f64 that
//! difference cancels when |r| ≪ |s|; the test shows the cancellation is real
//! (the naive f64 formula misses the tolerance on Sun rows) and that the
//! crate's formulation (Battin's f(q)) meets it on every row.
//!
//! Author: Moroya Sakamoto

use alice_det_math::{max64, sqrt64};
use alice_space::gravity::third_body;

const EPS: f64 = f64::EPSILON;

fn norm(v: [f64; 3]) -> f64 {
    sqrt64(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
}

/// The direct difference in f64, as a careless implementation writes it.
fn naive(r: [f64; 3], s: [f64; 3], mu: f64) -> [f64; 3] {
    let d = [r[0] - s[0], r[1] - s[1], r[2] - s[2]];
    let (nd, ns) = (norm(d), norm(s));
    let (kd, ks) = (nd * nd * nd, ns * ns * ns);
    [
        -mu * (d[0] / kd + s[0] / ks),
        -mu * (d[1] / kd + s[1] / ks),
        -mu * (d[2] / kd + s[2] / ks),
    ]
}

/// q, f(q), the scaled sum and |r − s|³ each round in a few operations: the
/// result is within 32 ε of |a| (component errors measured against the
/// vector's magnitude).
#[test]
fn third_body_matches_the_60_digit_difference() {
    let text = std::fs::read_to_string(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/golden/third_body.csv"
    ))
    .unwrap();
    let mut compared = 0;
    let mut naive_missed = 0;
    let mut worst: f64 = 0.0;
    for line in text
        .lines()
        .filter(|l| l.starts_with("sun,") || l.starts_with("moon,"))
    {
        let v: Vec<f64> = line
            .split(',')
            .skip(1)
            .map(|x| x.parse().unwrap())
            .collect();
        let (mu, r, s, want) = (
            v[0],
            [v[1], v[2], v[3]],
            [v[4], v[5], v[6]],
            [v[7], v[8], v[9]],
        );
        let tol = 32.0 * EPS * norm(want);
        let got = third_body(r, s, mu);
        let n = naive(r, s, mu);
        let mut naive_err: f64 = 0.0;
        for (c, ((g, w), m)) in got.iter().zip(&want).zip(&n).enumerate() {
            let e = (g - w).abs();
            assert!(e <= tol, "{line}: {c}: {g} vs {w} ({e:e} > {tol:e})");
            worst = max64(worst, e / tol);
            naive_err = max64(naive_err, (m - w).abs());
        }
        naive_missed += usize::from(naive_err > tol);
        compared += 1;
    }
    assert_eq!(compared, 144);
    // the table does catch the cancellation of the direct difference
    assert!(naive_missed > 50, "naive missed only {naive_missed}");
    assert!(worst > 0.0);
}

/// The body at the satellite's position, at the centre, or a NaN input give
/// the canonical NaN (the same bits on every target: an arithmetic NaN has the
/// sign bit set on x86 and not on ARM) rather than a finite value; the
/// satellite at the centre feels no tidal acceleration
/// (r = 0: a = −μ₃(−s/|s|³ + s/|s|³) = 0).
#[test]
fn degenerate_positions() {
    let s = [149_597_870.7, 0.0, 0.0];
    let mu = 1.327_124_400_412_794e11;
    assert_eq!(third_body([0.0; 3], s, mu), [0.0; 3]);
    assert!(third_body(s, s, mu)
        .iter()
        .all(|v| v.to_bits() == f64::NAN.to_bits()));
    assert!(third_body([7000.0, 0.0, 0.0], [0.0; 3], mu)
        .iter()
        .all(|v| v.to_bits() == f64::NAN.to_bits()));
    assert!(third_body([f64::NAN, 0.0, 0.0], s, mu)
        .iter()
        .all(|v| v.to_bits() == f64::NAN.to_bits()));
}

/// A position whose squares overflow (|r|² = ∞, so 1 + q and |r − s|² are ∞ and
/// f(q) is ∞/∞) or an infinite one gives the canonical NaN too: the NaN made
/// by ∞/∞ has the sign bit set on x86 and not on ARM, so without the explicit
/// value the output bits differ by target
#[test]
fn overflowing_positions_give_the_canonical_nan() {
    let s = [149_597_870.7, 0.0, 0.0];
    let mu = 1.327_124_400_412_794e11;
    for r in [
        [1e160, 0.0, 0.0],
        [0.0, 1e200, 0.0],
        [f64::INFINITY, 0.0, 0.0],
        [0.0, 0.0, f64::NEG_INFINITY],
    ] {
        let a = third_body(r, s, mu);
        assert!(
            a.iter().all(|v| v.to_bits() == f64::NAN.to_bits()),
            "{r:?}: {:?}",
            a.map(|v| format!("{:#018x}", v.to_bits()))
        );
    }
    // a NaN μ with a payload of its own: the result is the canonical NaN, not that payload
    let payload = f64::from_bits(0x7ff8_0000_0000_0001);
    let a = third_body([7000.0, 0.0, 0.0], s, payload);
    assert!(
        a.iter().all(|v| v.to_bits() == f64::NAN.to_bits()),
        "{:?}",
        a.map(|v| format!("{:#018x}", v.to_bits()))
    );
}
