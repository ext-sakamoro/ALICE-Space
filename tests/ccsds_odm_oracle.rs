//! Oracles for `ccsds_odm`: an OEM tabulated from a closed-form orbit is read
//! and interpolated back to that orbit.
//!
//! The trajectory is a circular orbit r(t) = R (cos ωt, sin ωt, 0) with R =
//! 7000 km and a 5800 s period, written as an OEM every 60 s (epoch strings
//! formed from whole seconds, so they are exact). Interpolating the read table
//! with the method it names must reproduce the closed form at times between the
//! rows within the interpolation error bound of that method (see
//! `tests/interp_oracle.rs`), and a table with every other row removed must
//! reproduce the removed rows. The message also carries what a reader must pass
//! over: comments, unit annotations, a covariance block and a second segment.
//!
//! Author: Moroya Sakamoto

// the index pairs the same axis of two 3-vectors
#![allow(clippy::needless_range_loop)]

use alice_space::ccsds_odm::{self, Interpolation, OdmEpoch, OdmError};

const R: f64 = 7000.0;

fn omega() -> f64 {
    2.0 * core::f64::consts::PI / 5800.0
}

fn circle(t: f64) -> ([f64; 3], [f64; 3]) {
    let w = omega();
    let (s, c) = (w * t).sin_cos();
    ([R * c, R * s, 0.0], [-R * w * s, R * w * c, 0.0])
}

/// `2026-10-23T00:00:00` plus `seconds` (whole seconds, past midnight into the
/// next days), as text.
fn epoch_text(seconds: i64) -> String {
    let (day, sod) = (23 + seconds / 86_400, seconds % 86_400);
    let (h, m, s) = (sod / 3600, sod / 60 % 60, sod % 60);
    format!("2026-10-{day:02}T{h:02}:{m:02}:{s:02}.000000")
}

/// One OEM segment from `t = 0` to `stop` s every `step` s, keeping the rows
/// for which `keep(i)` holds.
fn segment(
    name: &str,
    stop: i64,
    step: i64,
    method: &str,
    degree: usize,
    keep: impl Fn(i64) -> bool,
) -> String {
    segment_from(name, 0, stop, step, method, degree, keep)
}

/// [`segment`] starting at `t0` s instead of 0.
fn segment_from(
    name: &str,
    t0: i64,
    stop: i64,
    step: i64,
    method: &str,
    degree: usize,
    keep: impl Fn(i64) -> bool,
) -> String {
    let mut s = format!(
        "META_START\nOBJECT_NAME = {name}\nOBJECT_ID = 2026-001A\nCENTER_NAME = EARTH\n\
         REF_FRAME = EME2000\nTIME_SYSTEM = UTC\nSTART_TIME = {}\nUSEABLE_START_TIME = {}\n\
         USEABLE_STOP_TIME = {}\nSTOP_TIME = {}\nINTERPOLATION = {method}\n\
         INTERPOLATION_DEGREE = {degree}\nMETA_STOP\nCOMMENT the rows of a circular orbit\n",
        epoch_text(t0),
        epoch_text(t0),
        epoch_text(t0 + stop),
        epoch_text(t0 + stop),
    );
    for i in 0..=stop / step {
        if !keep(i) {
            continue;
        }
        let t = i * step;
        #[allow(clippy::cast_precision_loss)]
        let (p, v) = circle(t as f64);
        let t = t0 + t;
        s.push_str(&format!(
            "{} {:.9} {:.9} {:.9} {:.12} {:.12} {:.12}\n",
            epoch_text(t),
            p[0],
            p[1],
            p[2],
            v[0],
            v[1],
            v[2]
        ));
    }
    s
}

fn message(segments: &[String]) -> String {
    let mut s = String::from(
        "CCSDS_OEM_VERS = 2.0\nCOMMENT test message\nCREATION_DATE = 2026-10-23T00:00:00\nORIGINATOR = TEST\n",
    );
    for seg in segments {
        s.push_str(seg);
        s.push_str(
            "COVARIANCE_START\nEPOCH = 2026-10-23T00:00:00\nCOV_REF_FRAME = EME2000\n\
             1.0e-3\n1.0e-6 1.0e-3\nCOVARIANCE_STOP\n",
        );
    }
    s
}

#[test]
fn an_oem_of_a_circular_orbit_interpolates_back_to_it() {
    let text = message(&[
        segment("SAT-L", 3600, 60, "LAGRANGE", 7, |_| true),
        segment("SAT-H", 3600, 60, "HERMITE", 5, |_| true),
    ]);
    let oem = ccsds_odm::parse_oem_kvn(&text).expect("reads");
    assert_eq!(oem.version, "2.0");
    assert_eq!(oem.segments.len(), 2);
    assert_eq!(
        oem.segments[0].interpolation,
        Some((Interpolation::Lagrange, 7))
    );
    assert_eq!(
        oem.segments[1].interpolation,
        Some((Interpolation::Hermite, 5))
    );
    assert_eq!(oem.segments[0].epochs.len(), 61);

    let w = omega();
    // printed to 1e-9 km / 1e-12 km/s: the rounding of the rows adds ~1e-9 km
    let lag_bound = R * (w * 60.0).powi(8) / 4.0 / 8.0 + 2e-8;
    let her_bound = R * (w * 60.0).powi(6) * 4.0 / 720.0 + 2e-8;
    let mut compared = 0;
    for k in 0..59 {
        let t = 60.0 * f64::from(k) + 23.7;
        let e = OdmEpoch {
            day: oem.segments[0].epochs[0].day,
            second_of_day: t,
        };
        let exact = circle(t).0;
        for (seg, bound) in [(&oem.segments[0], lag_bound), (&oem.segments[1], her_bound)] {
            let (p, _) = seg.state_at(&e).expect("interpolates");
            for c in 0..3 {
                assert!(
                    (p[c] - exact[c]).abs() <= bound,
                    "{}: t={t} c={c}: {}",
                    seg.object_name,
                    (p[c] - exact[c]).abs()
                );
                compared += 1;
            }
        }
    }
    assert_eq!(compared, 59 * 2 * 3);
}

#[test]
fn a_decimated_oem_recovers_the_removed_rows() {
    let full =
        ccsds_odm::parse_oem_kvn(&message(&[segment("SAT", 3600, 60, "LAGRANGE", 7, |_| {
            true
        })]))
        .unwrap();
    let half =
        ccsds_odm::parse_oem_kvn(&message(&[segment("SAT", 3600, 60, "LAGRANGE", 7, |i| {
            i % 2 == 0
        })]))
        .unwrap();
    let (f, h) = (&full.segments[0], &half.segments[0]);
    assert_eq!(h.epochs.len(), 31);
    let w = omega();
    // spacing 120 s in the decimated table
    let bound = R * (w * 120.0).powi(8) / 4.0 / 8.0 + 2e-8;
    let mut compared = 0;
    for i in (1..60).step_by(2) {
        let (p, _) = h.state_at(&f.epochs[i]).expect("interpolates");
        for c in 0..3 {
            assert!(
                (p[c] - f.positions_km[i][c]).abs() <= bound,
                "row {i} c={c}: {}",
                (p[c] - f.positions_km[i][c]).abs()
            );
            compared += 1;
        }
    }
    assert_eq!(compared, 30 * 3);
}

#[test]
fn an_opm_reads_its_state_and_elements() {
    let text = "\
CCSDS_OPM_VERS = 2.0
CREATION_DATE = 2026-10-23T00:00:00
ORIGINATOR = TEST
OBJECT_NAME = SAT
OBJECT_ID = 2026-001A
CENTER_NAME = EARTH
REF_FRAME = EME2000
TIME_SYSTEM = UTC
COMMENT state vector
EPOCH = 2026-296T12:00:00.5
X = 7000.0 [km]
Y = -1.5 [km]
Z = 0.25 [km]
X_DOT = 0.001 [km/s]
Y_DOT = 7.546 [km/s]
Z_DOT = -0.002 [km/s]
SEMI_MAJOR_AXIS = 7000.1 [km]
ECCENTRICITY = 0.001
INCLINATION = 51.6 [deg]
RA_OF_ASC_NODE = 10.0 [deg]
ARG_OF_PERICENTER = 20.0 [deg]
TRUE_ANOMALY = 30.0 [deg]
GM = 398600.4415 [km**3/s**2]
MASS = 420000 [kg]
";
    let opm = ccsds_odm::parse_opm_kvn(text).expect("reads");
    assert_eq!(opm.position_km, [7000.0, -1.5, 0.25]);
    assert_eq!(opm.velocity_km_s, [0.001, 7.546, -0.002]);
    // day 296 of 2026 is October 23
    assert_eq!(opm.epoch, OdmEpoch::parse("2026-10-23T12:00:00.5").unwrap());
    let k = opm.keplerian.expect("elements");
    assert_eq!(k.anomaly_deg, (30.0, true));
    assert_eq!(k.gm_km3_s2, 398_600.441_5);
}

#[test]
fn malformed_and_out_of_window_inputs_are_errors() {
    let ok = message(&[segment("SAT", 600, 60, "LAGRANGE", 3, |_| true)]);
    let oem = ccsds_odm::parse_oem_kvn(&ok).unwrap();
    let seg = &oem.segments[0];
    let day = seg.epochs[0].day;
    // outside the useable window
    assert_eq!(
        seg.state_at(&OdmEpoch {
            day,
            second_of_day: 601.0
        }),
        Err(OdmError::OutsideUseableWindow)
    );
    assert_eq!(
        seg.state_at(&OdmEpoch {
            day: day - 1,
            second_of_day: 0.0
        }),
        Err(OdmError::OutsideUseableWindow)
    );
    // a degree the table cannot support
    assert!(matches!(
        seg.state_at_with(&seg.epochs[1], Interpolation::Lagrange, 20),
        Err(OdmError::Interp(_))
    ));
    // no interpolation named
    let none = ok.replace("INTERPOLATION = LAGRANGE\nINTERPOLATION_DEGREE = 3\n", "");
    assert_eq!(
        ccsds_odm::parse_oem_kvn(&none).unwrap().segments[0].state_at(&seg.epochs[1]),
        Err(OdmError::NoInterpolation)
    );
    // missing version, missing metadata key, unsupported interpolation, missing degree
    assert_eq!(
        ccsds_odm::parse_oem_kvn(&ok.replace("CCSDS_OEM_VERS = 2.0\n", "")),
        Err(OdmError::Missing("CCSDS_OEM_VERS"))
    );
    assert_eq!(
        ccsds_odm::parse_oem_kvn(&ok.replace("REF_FRAME = EME2000\n", "")),
        Err(OdmError::Missing("REF_FRAME"))
    );
    assert!(matches!(
        ccsds_odm::parse_oem_kvn(&ok.replace("= LAGRANGE", "= SPLINE")),
        Err(OdmError::Unsupported(..))
    ));
    assert_eq!(
        ccsds_odm::parse_oem_kvn(&ok.replace("INTERPOLATION_DEGREE = 3\n", "")),
        Err(OdmError::Missing("INTERPOLATION_DEGREE"))
    );
    // a short row, an unordered table, a non-number, a leap second in UTC
    let rows: Vec<&str> = ok
        .lines()
        .filter(|l| l.starts_with("2026-10-23T00:0"))
        .collect();
    let short = ok.replacen(rows[2], &rows[2][..rows[2].len() - 20], 1);
    assert!(matches!(
        ccsds_odm::parse_oem_kvn(&short),
        Err(OdmError::Syntax(..))
    ));
    let swapped = ok
        .replacen(rows[2], "TMP", 1)
        .replacen(rows[3], rows[2], 1)
        .replacen("TMP", rows[3], 1);
    assert!(matches!(
        ccsds_odm::parse_oem_kvn(&swapped),
        Err(OdmError::Syntax(..))
    ));
    let nan = ok.replacen(rows[1], &rows[1].replacen(" 0.000000000 ", " nan ", 1), 1);
    assert!(
        matches!(ccsds_odm::parse_oem_kvn(&nan), Err(OdmError::Invalid(..))),
        "{nan}"
    );
    let leap = "CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = S\nOBJECT_ID = X\nCENTER_NAME = EARTH\nREF_FRAME = EME2000\n\
         TIME_SYSTEM = UTC\nSTART_TIME = 2016-12-31T23:59:59\nSTOP_TIME = 2017-01-01T00:00:01\nMETA_STOP\n\
         2016-12-31T23:59:59 1 2 3 4 5 6\n2016-12-31T23:59:60 1 2 3 4 5 6\n2017-01-01T00:00:01 1 2 3 4 5 6\n";
    assert!(matches!(
        ccsds_odm::parse_oem_kvn(leap),
        Err(OdmError::Unsupported(..))
    ));
    // OdmEpoch::parse reads second 60 as a label only; whether it exists is
    // decided by the time system (tests/ccsds_odm_leap_second.rs)
    assert!(OdmEpoch::parse("2016-12-31T23:59:60").is_some());
    assert!(
        OdmEpoch::parse("2016-12-31T23:58:60").is_none(),
        "second 60 only at 23:59"
    );
    // data outside START_TIME..STOP_TIME
    let outside = ok.replace(
        &format!("STOP_TIME = {}", epoch_text(600)),
        &format!("STOP_TIME = {}", epoch_text(300)),
    );
    assert!(matches!(
        ccsds_odm::parse_oem_kvn(&outside),
        Err(OdmError::Invalid(..))
    ));
    // META_START without META_STOP, a key among the data, an empty message
    assert!(matches!(
        ccsds_odm::parse_oem_kvn("CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = S\n"),
        Err(OdmError::Syntax(..))
    ));
    assert!(matches!(
        ccsds_odm::parse_oem_kvn(&ok.replacen(rows[1], "FOO = 1", 1)),
        Err(OdmError::Syntax(..))
    ));
    assert!(ccsds_odm::parse_oem_kvn("").is_err());
    // OPM: a line that is not KEY = value, a missing coordinate
    assert!(matches!(
        ccsds_odm::parse_opm_kvn("CCSDS_OPM_VERS = 2.0\nnonsense\n"),
        Err(OdmError::Syntax(..))
    ));
    assert!(matches!(
        ccsds_odm::parse_opm_kvn(
            "CCSDS_OPM_VERS = 2.0\nEPOCH = 2026-10-23T00:00:00\nTIME_SYSTEM = UTC\n"
        ),
        Err(OdmError::Missing(_))
    ));
}

/// A segment from 23:30 to 01:00 the next day: times relative to the segment
/// count the day change as 86400 s (the closed form is in seconds from the
/// segment start).
#[test]
fn a_segment_across_midnight_interpolates_back_to_the_orbit() {
    let t0 = 84_600; // 23:30:00
    let oem = ccsds_odm::parse_oem_kvn(&message(&[segment_from(
        "SAT",
        t0,
        5400,
        60,
        "LAGRANGE",
        7,
        |_| true,
    )]))
    .unwrap();
    let seg = &oem.segments[0];
    assert_eq!(
        seg.epochs[0].day + 1,
        seg.epochs[seg.epochs.len() - 1].day,
        "the segment spans two days"
    );
    let w = omega();
    let bound = R * (w * 60.0).powi(8) / 4.0 / 8.0 + 2e-8;
    let mut compared = 0;
    for k in 0..89 {
        let t = 60.0 * f64::from(k) + 23.7; // seconds from the segment start
        let sod = 84_600.0 + t;
        let e = if sod < 86_400.0 {
            OdmEpoch {
                day: seg.epochs[0].day,
                second_of_day: sod,
            }
        } else {
            OdmEpoch {
                day: seg.epochs[0].day + 1,
                second_of_day: sod - 86_400.0,
            }
        };
        let (p, _) = seg.state_at(&e).expect("interpolates");
        let exact = circle(t).0;
        for c in 0..3 {
            assert!(
                (p[c] - exact[c]).abs() <= bound,
                "t={t} c={c}: {}",
                (p[c] - exact[c]).abs()
            );
            compared += 1;
        }
    }
    assert_eq!(compared, 89 * 3);
}
