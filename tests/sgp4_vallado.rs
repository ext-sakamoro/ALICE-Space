//! Oracle: Vallado's SGP4 verification set (AIAA 2006-6753).
//!
//! `tests/data/sgp4/SGP4-VER.TLE` holds the 33 verification element sets, each
//! line 2 followed by a start / stop / step (minutes); `tcppver.out` is the
//! output of Vallado's C++ reference for them. Both files come unchanged from
//! python-sgp4 (MIT, commit 8126f773, see `tests/data/sgp4/LICENSE-python-sgp4`)
//! and are pinned here by SHA-256.
//!
//! The walk is that of the reference driver (`testcpp`, as python-sgp4's
//! `generate_test_output` reproduces it): for each set, the state at 0, then from
//! start to stop by step (0 is not repeated), then the stop time itself when the
//! step does not land on it; the first error ends the set. Each state (TEME km,
//! km/s) is compared with the reference line field by field within 2e-7 (the
//! tolerance python-sgp4 uses against the same file; the file prints positions
//! to 1e-8 km and velocities to 1e-9 km/s). The other columns of the reference
//! (osculating elements and the calendar date) are not compared.
//!
//! The sequence of errors must be the reference's: codes 1, 1, 6, 6, 4, 3, 6.
//! When a set fails at the epoch, the reference prints the previous set's state
//! with time 0 (it has no new state); that line is skipped, not compared.
//!
//! Three error-case sets (33333, 33334, 33335) carry check digits that do not
//! match their lines in the published file; the TLE parser verifies check
//! digits, so the harness recomputes the digit of every line before parsing
//! (the elements are unchanged).
//!
//! Author: Moroya Sakamoto

use alice_space::sgp4::{GravityModel, MeanElements, OpsMode, Sgp4};
use alice_space::tle;
use sha2::{Digest, Sha256};

const TLE_FILE: &str = include_str!("data/sgp4/SGP4-VER.TLE");
const REFERENCE_FILE: &str = include_str!("data/sgp4/tcppver.out");

const TLE_SHA256: &str = "d246d1d9d768ace445a38a965713fa9ba52d80fd8a41a0502ff83d7acffe2881";
const REFERENCE_SHA256: &str = "687bf28dbe52df86e8e60ab5cb4a08d1aa3dbcaf4e63b1f7ab95f044fbe3833b";

const TOLERANCE: f64 = 2e-7;

fn sha256_hex(bytes: &[u8]) -> String {
    let digest: [u8; 32] = Sha256::digest(bytes).into();
    digest.iter().map(|b| format!("{b:02x}")).collect()
}

/// `line[..68]` followed by its recomputed check digit.
fn with_check_digit(line: &str) -> String {
    let payload = &line[..68];
    let sum: u32 = payload
        .chars()
        .map(|c| c.to_digit(10).unwrap_or(u32::from(c == '-')))
        .sum();
    format!("{payload}{}", sum % 10)
}

/// The verification sets: (line 1, line 2 cut to 69 columns, start, stop, step).
fn verification_sets() -> Vec<(String, String, f64, f64, f64)> {
    let mut out = Vec::new();
    let mut lines = TLE_FILE.lines().map(|l| l.trim_end_matches('\r'));
    while let Some(l1) = lines.next() {
        if !l1.starts_with("1 ") {
            continue;
        }
        let l2 = lines.next().expect("line 2 follows line 1");
        let times: Vec<f64> = l2[69..]
            .split_whitespace()
            .map(|f| f.parse().expect("start / stop / step"))
            .collect();
        assert_eq!(times.len(), 3, "start / stop / step after line 2");
        out.push((
            with_check_digit(&l1[..69]),
            with_check_digit(&l2[..69]),
            times[0],
            times[1],
            times[2],
        ));
    }
    out
}

/// One state as the reference prints it: tsince, position, velocity.
type Row = [f64; 7];

/// Our walk over the verification set: per set, the satellite number, the rows,
/// and the first error code (if any) with whether it happened at the epoch.
/// Per set: the satellite number and its rows.
type Sets = Vec<(u32, Vec<Row>)>;

fn our_walk() -> (Sets, Vec<(u8, bool)>, Vec<u8>) {
    let mut sets = Vec::new();
    let mut errors = Vec::new();
    let mut bits = Vec::new();
    for (l1, l2, tstart, tend, tstep) in verification_sets() {
        let parsed = tle::parse(&format!("{l1}\n{l2}")).expect("verification TLE parses");
        let sat = Sgp4::new(
            &MeanElements::from_tle(&parsed),
            GravityModel::Wgs72,
            OpsMode::Improved,
        )
        .expect("verification set initialises");
        let mut rows = Vec::new();
        let mut row_at = |t: f64, rows: &mut Vec<Row>| match sat.propagate(t) {
            Ok(s) => {
                let r = s.position_km;
                let v = s.velocity_km_s;
                for x in [t, r[0], r[1], r[2], v[0], v[1], v[2]] {
                    bits.extend_from_slice(&x.to_bits().to_le_bytes());
                }
                rows.push([t, r[0], r[1], r[2], v[0], v[1], v[2]]);
                None
            }
            Err(e) => Some(e.code()),
        };
        if let Some(code) = row_at(0.0, &mut rows) {
            errors.push((code, true));
            sets.push((parsed.satellite_number, rows));
            continue;
        }
        let mut tsince = tstart;
        let mut failed = false;
        while tsince <= tend {
            if tsince == tstart && tstart == 0.0 {
                tsince += tstep;
                continue;
            }
            if let Some(code) = row_at(tsince, &mut rows) {
                errors.push((code, false));
                failed = true;
                break;
            }
            tsince += tstep;
        }
        if !failed && tsince - tend < tstep - 1e-6 {
            if let Some(code) = row_at(tend, &mut rows) {
                errors.push((code, false));
            }
        }
        sets.push((parsed.satellite_number, rows));
    }
    (sets, errors, bits)
}

/// The reference file: per set, the satellite number and its data rows (the
/// first 7 fields of each line).
fn reference() -> Sets {
    let mut out: Sets = Vec::new();
    for line in REFERENCE_FILE.lines().map(|l| l.trim_end_matches('\r')) {
        if line.trim().is_empty() {
            continue;
        }
        if line.trim_end().ends_with("xx") {
            let n = line
                .split_whitespace()
                .next()
                .expect("number")
                .parse()
                .expect("number");
            out.push((n, Vec::new()));
            continue;
        }
        let fields: Vec<f64> = line
            .split_whitespace()
            .take(7)
            .map(|f| f.parse().expect("number"))
            .collect();
        let row: Row = fields.try_into().expect("7 fields");
        out.last_mut().expect("a set header first").1.push(row);
    }
    out
}

#[test]
fn the_verification_files_are_the_published_ones() {
    assert_eq!(sha256_hex(TLE_FILE.as_bytes()), TLE_SHA256);
    assert_eq!(sha256_hex(REFERENCE_FILE.as_bytes()), REFERENCE_SHA256);
}

#[test]
fn every_state_matches_vallados_reference_output() {
    let (ours, errors, _) = our_walk();
    let theirs = reference();
    assert_eq!(ours.len(), 33, "33 verification sets");
    assert_eq!(theirs.len(), ours.len(), "the reference has the same sets");

    let mut compared = 0usize;
    let mut worst = 0.0f64;
    let mut differing: Vec<String> = Vec::new();
    let mut epoch_failures = errors.iter().filter(|e| e.1).count();
    for ((sat, our_rows), (ref_sat, ref_rows)) in ours.iter().zip(&theirs) {
        assert_eq!(sat, ref_sat, "the sets come in the same order");
        let ref_rows: &[Row] = if our_rows.is_empty() && !ref_rows.is_empty() && epoch_failures > 0
        {
            // failed at the epoch: the reference reprints the previous state
            epoch_failures -= 1;
            assert_eq!(ref_rows.len(), 1, "set {sat}: only the reprinted line");
            &[]
        } else {
            ref_rows
        };
        assert_eq!(
            our_rows.len(),
            ref_rows.len(),
            "set {sat}: number of states"
        );
        for (a, e) in our_rows.iter().zip(ref_rows) {
            compared += 1;
            for (k, (x, y)) in a.iter().zip(e).enumerate() {
                let d = (x - y).abs();
                worst = worst.max(d);
                if d >= TOLERANCE {
                    differing.extend([format!("set {sat} t={} field {k}: {x} vs {y}", a[0])]);
                }
            }
        }
    }
    let codes: Vec<u8> = errors.iter().map(|e| e.0).collect();
    assert_eq!(
        codes,
        vec![1, 1, 6, 6, 4, 3, 6],
        "the reference's error sequence"
    );
    eprintln!("compared {compared} states, worst field difference {worst:e}");
    assert!(compared >= 600, "compared {compared} states");
    assert!(
        differing.is_empty(),
        "{} fields differ beyond {TOLERANCE} (worst {worst}):\n{}",
        differing.len(),
        differing
            .iter()
            .take(20)
            .cloned()
            .collect::<Vec<_>>()
            .join("\n")
    );
}

/// SHA-256 of the bits of every state of the walk (t, position, velocity as
/// little-endian `f64`). Every transcendental goes through alice-det-math, so
/// the walk gives these bits on every target.
const WALK_BITS_SHA256: &str = "0a767b6dd442e8eedb0a0923bcfa38a71f709c010223ba46bd4868f2ffb980ff";

#[test]
fn the_walk_gives_the_same_bits_on_every_target() {
    let (_, _, bits) = our_walk();
    assert!(!bits.is_empty(), "the walk produced states");
    assert_eq!(sha256_hex(&bits), WALK_BITS_SHA256);
}

#[test]
fn propagate_is_a_pure_function_of_the_time() {
    // a deep-space resonant set (12 h, e ≈ 0.7): the integrator restarts from the
    // epoch on every call, so the order of calls does not matter
    for (l1, l2, _, _, _) in verification_sets() {
        let parsed = tle::parse(&format!("{l1}\n{l2}")).expect("parses");
        let sat = Sgp4::new(
            &MeanElements::from_tle(&parsed),
            GravityModel::Wgs72,
            OpsMode::Improved,
        )
        .expect("verification set initialises");
        let forward: Vec<_> = [0.0, 720.0, 1440.0, -1440.0]
            .iter()
            .map(|&t| sat.propagate(t))
            .collect();
        let backward: Vec<_> = [-1440.0, 1440.0, 720.0, 0.0]
            .iter()
            .map(|&t| sat.propagate(t))
            .rev()
            .collect();
        assert_eq!(forward, backward, "set {}", parsed.satellite_number);
    }
}

#[test]
fn the_module_calls_no_platform_transcendental() {
    let src = include_str!("../src/sgp4.rs");
    assert!(src.len() > 1000, "the module source was read");
    for call in [
        ".sin(",
        ".cos(",
        ".tan(",
        ".atan2(",
        ".atan(",
        ".asin(",
        ".acos(",
        ".powf(",
        ".powi(",
        ".exp(",
        ".ln(",
        ".sqrt(",
        ".cbrt(",
        ".mul_add(",
        ".sin_cos(",
    ] {
        assert!(!src.contains(call), "src/sgp4.rs calls {call}");
    }
}
