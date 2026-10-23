//! Oracle: an SP3 file tabulated from a closed-form orbit interpolates back to
//! that orbit.
//!
//! A GPS-like circular orbit (R = 26 560 km, period 43 082 s) is written as an
//! SP3 file every 15 minutes, positions printed to 1e-6 km as SP3 does, with one
//! epoch marked missing (all-zero position). Lagrange of degree 9 (the degree
//! used for precise orbits) at times between the epochs must reproduce the
//! closed form within the error bound of the centred window plus the printing
//! rounding (5e-6 km), also across the missing epoch.
//!
//! Author: Moroya Sakamoto

// the index pairs the same axis of two 3-vectors
#![allow(clippy::needless_range_loop)]

use alice_space::sp3;

const R: f64 = 26_560.0;
const PERIOD: f64 = 43_082.0;
/// 2020-05-24 00:00:00 UTC.
const START_UNIX: u64 = 1_590_278_400;

fn circle(t: f64) -> [f64; 3] {
    let w = 2.0 * core::f64::consts::PI / PERIOD;
    let (s, c) = (w * t).sin_cos();
    [R * c, R * s, 0.5 * R * s]
}

fn sp3_text(epochs: usize, missing: usize) -> String {
    let mut s = String::from(
        "#cP2020  5 24  0  0  0.00000000      96 ORBIT IGb14 HLM  IGS\n+    1   G01\n",
    );
    for i in 0..epochs {
        let t = 900 * i;
        let (h, m) = (t / 3600, t / 60 % 60);
        let day = 24 + h / 24;
        s.push_str(&format!(
            "*  2020  5 {day:2} {:2} {m:2}  0.00000000\n",
            h % 24
        ));
        #[allow(clippy::cast_precision_loss)]
        let p = if i == missing {
            [0.0; 3]
        } else {
            circle(t as f64)
        };
        s.push_str(&format!(
            "PG01{:14.6}{:14.6}{:14.6}{:14.6}\n",
            p[0], p[1], p[2], 1.5
        ));
    }
    s.push_str("EOF\n");
    s
}

#[test]
fn an_sp3_orbit_interpolates_back_to_the_closed_form() {
    let missing = 20;
    let file = sp3::parse(&sp3_text(48, missing)).expect("parses");
    assert_eq!(file.epochs.len(), 48);
    assert_eq!(file.epochs[0].unix_seconds, START_UNIX);
    let w = 2.0 * core::f64::consts::PI / PERIOD;
    let mut compared = 0;
    for k in 5..42 {
        let t = 900.0 * f64::from(k) + 333.0;
        #[allow(clippy::cast_precision_loss)]
        let p = file
            .interpolate("G01", START_UNIX as f64 + t, 9)
            .expect("interpolates");
        let exact = circle(t);
        // centred 10-node window (9 h^10 spread), less than 1e-6 km of truncation
        // for this orbit; the printed rows round to 5e-7 km
        let bound = R * (w * 900.0).powi(10) * 1000.0 / 3_628_800.0 + 5e-6;
        for c in 0..3 {
            assert!(
                (p[c] - exact[c]).abs() <= bound,
                "t={t} c={c}: {}",
                (p[c] - exact[c]).abs()
            );
            compared += 1;
        }
    }
    assert_eq!(compared, 37 * 3);
    // the missing epoch is not a node: interpolating at it gives the orbit, not 0
    #[allow(clippy::cast_precision_loss)]
    let at_missing = file
        .interpolate("G01", START_UNIX as f64 + 900.0 * missing as f64, 9)
        .unwrap();
    #[allow(clippy::cast_precision_loss)]
    let exact = circle(900.0 * missing as f64);
    assert!((at_missing[0] - exact[0]).abs() < 1e-5);
    // an unknown satellite and a time outside the file are errors
    assert!(file
        .interpolate("G02", START_UNIX as f64 + 1000.0, 9)
        .is_err());
    assert!(file.interpolate("G01", START_UNIX as f64 - 1.0, 9).is_err());
}
