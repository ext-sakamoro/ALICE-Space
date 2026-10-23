//! Oracle: second 60 in Orbit Data Messages.
//!
//! Only UTC has a 61st second (23:59:60 on the last day of a leap-second
//! month). In a UTC message it is rejected as unsupported wherever it occurs
//! (uniform seconds across it need the time-scale tables): an OPM epoch, an
//! OEM row, START / STOP / USEABLE times and the epoch asked of
//! `state_at`. In TAI, TT, GPS (and every other uniform time system) second
//! 60 does not exist and is invalid. Before this, a query or a window time of
//! 23:59:60 was read as 00:00:00 of the next day.
//!
//! Author: Moroya Sakamoto

use alice_space::ccsds_odm::{self, OdmEpoch, OdmError};

const OPM: &str = "\
CCSDS_OPM_VERS = 2.0
OBJECT_NAME = SAT
OBJECT_ID = 2026-001A
CENTER_NAME = EARTH
REF_FRAME = EME2000
TIME_SYSTEM = {ts}
EPOCH = {epoch}
X = 7000.0 [km]
Y = -1.5 [km]
Z = 0.25 [km]
X_DOT = 0.001 [km/s]
Y_DOT = 7.546 [km/s]
Z_DOT = -0.002 [km/s]
";

fn opm(ts: &str, epoch: &str) -> Result<ccsds_odm::Opm, OdmError> {
    ccsds_odm::parse_opm_kvn(&OPM.replace("{ts}", ts).replace("{epoch}", epoch))
}

/// A 3-row segment around the 2016-12-31 leap second; `meta` overrides keys.
fn oem(ts: &str, meta: &[(&str, &str)], rows: &[&str]) -> Result<ccsds_odm::Oem, OdmError> {
    let mut keys = vec![
        ("START_TIME", "2016-12-31T23:59:00"),
        ("STOP_TIME", "2017-01-01T00:01:00"),
        ("INTERPOLATION", "LAGRANGE"),
        ("INTERPOLATION_DEGREE", "2"),
    ];
    for (k, v) in meta {
        match keys.iter_mut().find(|(key, _)| key == k) {
            Some(slot) => slot.1 = v,
            None => keys.push((k, v)),
        }
    }
    let mut s = format!(
        "CCSDS_OEM_VERS = 2.0\nMETA_START\nOBJECT_NAME = S\nOBJECT_ID = X\nCENTER_NAME = EARTH\n\
         REF_FRAME = EME2000\nTIME_SYSTEM = {ts}\n"
    );
    for (k, v) in keys {
        s.push_str(&format!("{k} = {v}\n"));
    }
    s.push_str("META_STOP\n");
    for r in rows {
        s.push_str(&format!("{r} 1 2 3 4 5 6\n"));
    }
    ccsds_odm::parse_oem_kvn(&s)
}

const ROWS: [&str; 3] = [
    "2016-12-31T23:59:00",
    "2016-12-31T23:59:30",
    "2017-01-01T00:01:00",
];
const LEAP: &str = "2016-12-31T23:59:60";

#[test]
fn a_utc_opm_epoch_on_a_leap_second_is_unsupported() {
    assert!(matches!(opm("UTC", LEAP), Err(OdmError::Unsupported(k, _)) if k == "EPOCH"));
    assert!(opm("UTC", "2016-12-31T23:59:59.5").is_ok());
}

#[test]
fn a_utc_leap_second_in_a_window_time_is_unsupported() {
    assert!(oem("UTC", &[], &ROWS).is_ok(), "the base message is valid");
    for key in [
        "START_TIME",
        "STOP_TIME",
        "USEABLE_START_TIME",
        "USEABLE_STOP_TIME",
    ] {
        let got = oem("UTC", &[(key, LEAP)], &ROWS);
        assert!(
            matches!(&got, Err(OdmError::Unsupported(k, _)) if k == key),
            "{key}: {got:?}"
        );
    }
}

#[test]
fn a_utc_leap_second_query_is_unsupported_not_the_next_midnight() {
    let msg = oem("UTC", &[], &ROWS).unwrap();
    let seg = &msg.segments[0];
    let leap = OdmEpoch::parse(LEAP).unwrap();
    let got = seg.state_at(&leap);
    assert!(matches!(got, Err(OdmError::Unsupported(..))), "{got:?}");
    // the next midnight itself is an ordinary query
    let midnight = OdmEpoch::parse("2017-01-01T00:00:00").unwrap();
    assert!(seg.state_at(&midnight).is_ok());
}

#[test]
fn second_60_is_invalid_in_uniform_time_systems() {
    for ts in ["TAI", "TT", "GPS", "TDB"] {
        assert!(
            matches!(opm(ts, LEAP), Err(OdmError::Invalid(k, _)) if k == "EPOCH"),
            "{ts} OPM"
        );
        let rows = ["2016-12-31T23:59:00", LEAP, "2017-01-01T00:01:00"];
        assert!(
            matches!(oem(ts, &[], &rows), Err(OdmError::Invalid(..))),
            "{ts} row"
        );
        for key in ["START_TIME", "USEABLE_STOP_TIME"] {
            let got = oem(ts, &[(key, LEAP)], &ROWS);
            assert!(
                matches!(&got, Err(OdmError::Invalid(k, _)) if k == key),
                "{ts} {key}: {got:?}"
            );
        }
        let msg = oem(ts, &[], &ROWS).unwrap();
        let got = msg.segments[0].state_at(&OdmEpoch::parse(LEAP).unwrap());
        assert!(
            matches!(got, Err(OdmError::Invalid(..))),
            "{ts} query: {got:?}"
        );
    }
}
