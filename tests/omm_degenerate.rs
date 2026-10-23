//! Degenerate and malformed OMM inputs: each case states the expected result
//! (the element set, or the error naming the key or the syntax) and asserts it.
//!
//! Author: Moroya Sakamoto

use alice_space::omm::{self, OmmError};

/// A complete KVN element set (MARIO) with `{EPOCH}` / `{EXTRA}` slots.
const KVN: &str = "\
OBJECT_NAME = MARIO
OBJECT_ID = 1998-067UQ
CENTER_NAME = EARTH
REF_FRAME = TEME
TIME_SYSTEM = UTC
MEAN_ELEMENT_THEORY = SGP4
EPOCH = {EPOCH}
MEAN_MOTION = 15.99081912
ECCENTRICITY = .0014649
INCLINATION = 51.6242
RA_OF_ASC_NODE = 216.2930
ARG_OF_PERICENTER = 331.8976
MEAN_ANOMALY = 28.1241
NORAD_CAT_ID = 55123
BSTAR = .1568E-2
MEAN_MOTION_DOT = .787702E-2
MEAN_MOTION_DDOT = .29408E-3
{EXTRA}";

fn kvn(epoch: &str, extra: &str) -> String {
    KVN.replace("{EPOCH}", epoch).replace("{EXTRA}", extra)
}

fn kvn_without(key: &str) -> String {
    kvn("2023-04-25T10:45:30.642912", "")
        .lines()
        .filter(|l| !l.starts_with(&format!("{key} ")))
        .collect::<Vec<_>>()
        .join("\n")
}

#[test]
fn the_base_message_reads() {
    let sets = omm::parse_kvn(&kvn("2023-04-25T10:45:30.642912", "")).expect("reads");
    assert_eq!(sets.len(), 1);
    assert_eq!(sets[0].norad_cat_id, 55123);
    assert_eq!(sets[0].classification, 'U', "absent classification is U");
}

/// The day-of-year form of the epoch (2023 day 115) is the same instant as the
/// calendar form (2023-04-25).
#[test]
fn a_day_of_year_epoch_equals_the_calendar_epoch() {
    let a = omm::parse_kvn(&kvn("2023-04-25T10:45:30.642912", "")).unwrap();
    let b = omm::parse_kvn(&kvn("2023-115T10:45:30.642912", "")).unwrap();
    let c = omm::parse_kvn(&kvn("2023-115T10:45:30.642912Z", "")).unwrap();
    assert_eq!(
        a[0].epoch_days_since_1950.to_bits(),
        b[0].epoch_days_since_1950.to_bits()
    );
    assert_eq!(
        a[0].epoch_days_since_1950.to_bits(),
        c[0].epoch_days_since_1950.to_bits()
    );
}

/// Each required key, removed: `Missing` naming it.
#[test]
fn each_required_key_missing_is_named() {
    let keys = [
        "EPOCH",
        "MEAN_MOTION",
        "ECCENTRICITY",
        "INCLINATION",
        "RA_OF_ASC_NODE",
        "ARG_OF_PERICENTER",
        "MEAN_ANOMALY",
        "NORAD_CAT_ID",
        "BSTAR",
        "MEAN_MOTION_DOT",
        "MEAN_MOTION_DDOT",
    ];
    for key in keys {
        let err = omm::parse_kvn(&kvn_without(key)).unwrap_err();
        assert_eq!(err, OmmError::Missing(key), "{key}");
    }
}

/// Epochs that are not a time: month 13, February 29 of a non-leap year, a leap
/// second, day 366 of a non-leap year, hour 24, no `T`.
#[test]
fn an_invalid_epoch_is_rejected() {
    for epoch in [
        "2023-13-01T00:00:00",
        "2023-02-29T00:00:00",
        "2016-12-31T23:59:60",
        "2023-366T00:00:00",
        "2023-04-25T24:00:00",
        "2023-04-25 10:45:30",
        "2023-04-25T10:45:3x",
        "2023-04-25T10:45:30.",
    ] {
        let err = omm::parse_kvn(&kvn(epoch, "")).unwrap_err();
        assert_eq!(
            err,
            OmmError::Invalid("EPOCH".into(), epoch.into()),
            "{epoch}"
        );
    }
    // February 29 of a leap year and day 366 of a leap year are valid
    assert!(omm::parse_kvn(&kvn("2024-02-29T00:00:00", "")).is_ok());
    assert!(omm::parse_kvn(&kvn("2024-366T00:00:00", "")).is_ok());
}

/// Metadata SGP4 elements cannot have: `Unsupported` naming the key.
#[test]
fn unsupported_metadata_is_rejected() {
    for (key, value) in [
        ("CENTER_NAME", "MOON"),
        ("REF_FRAME", "GCRF"),
        ("TIME_SYSTEM", "TAI"),
        ("MEAN_ELEMENT_THEORY", "DSST"),
    ] {
        let text = kvn("2023-04-25T10:45:30.642912", "")
            .lines()
            .map(|l| {
                if l.starts_with(&format!("{key} ")) {
                    format!("{key} = {value}")
                } else {
                    l.to_owned()
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let err = omm::parse_kvn(&text).unwrap_err();
        assert_eq!(
            err,
            OmmError::Unsupported(key.into(), value.into()),
            "{key}"
        );
    }
}

/// Values out of range or not numbers: `Invalid` naming the key.
#[test]
fn invalid_values_are_rejected() {
    for (key, value) in [
        ("ECCENTRICITY", "1.0"),
        ("ECCENTRICITY", "-0.1"),
        ("MEAN_MOTION", "0"),
        ("MEAN_MOTION", "nan"),
        ("INCLINATION", "inf"),
        ("BSTAR", "1e999"),
        ("NORAD_CAT_ID", "-5"),
    ] {
        let text = kvn("2023-04-25T10:45:30.642912", "")
            .lines()
            .map(|l| {
                if l.starts_with(&format!("{key} ")) {
                    format!("{key} = {value}")
                } else {
                    l.to_owned()
                }
            })
            .collect::<Vec<_>>()
            .join("\n");
        let err = omm::parse_kvn(&text).unwrap_err();
        assert_eq!(
            err,
            OmmError::Invalid(key.into(), value.into()),
            "{key} = {value}"
        );
    }
}

/// `OBJECT_ID = UNKNOWN` (CCSDS 502.0-B-3 Table 4-2) and an empty id read as an
/// empty designator.
#[test]
fn an_unknown_or_empty_object_id_is_empty() {
    for id in ["UNKNOWN", ""] {
        let text = kvn("2023-04-25T10:45:30.642912", "")
            .replace("OBJECT_ID = 1998-067UQ", &format!("OBJECT_ID = {id}"));
        assert_eq!(omm::parse_kvn(&text).unwrap()[0].object_id, "");
    }
}

/// Two element sets in one KVN message (a new `OBJECT_NAME` after an `EPOCH`).
#[test]
fn a_kvn_message_with_two_sets_reads_both() {
    let one = kvn("2023-04-25T10:45:30.642912", "");
    let two = format!(
        "{one}\n{}",
        one.replace("NORAD_CAT_ID = 55123", "NORAD_CAT_ID = 55124")
    );
    let sets = omm::parse_kvn(&two).unwrap();
    assert_eq!(
        sets.iter().map(|s| s.norad_cat_id).collect::<Vec<_>>(),
        vec![55123, 55124]
    );
}

/// Malformed text in each encoding: a syntax error, never a panic.
#[test]
fn malformed_text_is_a_syntax_error() {
    assert!(matches!(omm::parse_kvn(""), Err(OmmError::Syntax(_))));
    assert!(matches!(
        omm::parse_kvn("EPOCH 2023"),
        Err(OmmError::Syntax(_))
    ));
    for json in [
        "",
        "[",
        "[{",
        "{\"EPOCH\":}",
        "[1, 2]",
        "{\"A\":\"\\u12\"}",
        "[]",
        "{} x",
    ] {
        assert!(
            matches!(omm::parse_json(json), Err(OmmError::Syntax(_))),
            "JSON {json:?}"
        );
    }
    assert!(matches!(
        omm::parse_xml("<omm></omm>"),
        Err(OmmError::Syntax(_))
    ));
    assert!(matches!(
        omm::parse_xml("<segment><EPOCH>x"),
        Err(OmmError::Syntax(_))
    ));
    assert!(matches!(omm::parse_csv(""), Err(OmmError::Syntax(_))));
    assert!(matches!(omm::parse_csv("A,B\n1"), Err(OmmError::Syntax(_))));
    assert!(matches!(omm::parse_csv("A,B\n"), Err(OmmError::Syntax(_))));
}
