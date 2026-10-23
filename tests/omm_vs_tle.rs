//! Oracle: an OMM and the TLE of the same element set give the same SGP4 input.
//!
//! The pairs are CelesTrak's: the MARIO element set (catalogue 55123, every field
//! non-zero, so a wrong unit conversion cannot hide behind a zero) as TLE, as the
//! NDM XML and as the CSV that CelesTrak serves for it, and Vanguard 1
//! (catalogue 5) as TLE and as the OMM fields of the same element set; both are
//! the test data of python-sgp4 (MIT, commit 8126f773). The KVN and JSON forms
//! here carry the same values as the XML, written in those encodings.
//!
//! Every mean element, the epoch included, must have the same bits as the TLE's:
//! both are decimals converted with the same expressions, and the OMM epoch (a
//! calendar time to the microsecond) goes through the same day-of-year epoch
//! arithmetic as the TLE's (both print the same instant here). The SGP4 states
//! of both, at the epoch and a day later, must then agree within 8 mm.
//!
//! Author: Moroya Sakamoto

use alice_space::omm::{self, Omm};
use alice_space::sgp4::{GravityModel, MeanElements, OpsMode, Sgp4};
use alice_space::tle;

const MARIO_TLE: &str = "\
1 55123U 98067UQ  23115.44827133  .00787702  29408-3  15680-2 0  9999
2 55123  51.6242 216.2930 0014649 331.8976  28.1241 15.99081912 18396";

const MARIO_XML: &str = r#"<?xml version="1.0" encoding="UTF-8"?>
<ndm xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:noNamespaceSchemaLocation="https://sanaregistry.org/r/ndmxml_unqualified/ndmxml-2.0.0-master-2.0.xsd">
<omm id="CCSDS_OMM_VERS" version="2.0">
<header><CREATION_DATE/><ORIGINATOR/></header><body><segment><metadata><OBJECT_NAME>MARIO</OBJECT_NAME><OBJECT_ID>1998-067UQ</OBJECT_ID><CENTER_NAME>EARTH</CENTER_NAME><REF_FRAME>TEME</REF_FRAME><TIME_SYSTEM>UTC</TIME_SYSTEM><MEAN_ELEMENT_THEORY>SGP4</MEAN_ELEMENT_THEORY></metadata><data><meanElements><EPOCH>2023-04-25T10:45:30.642912</EPOCH><MEAN_MOTION>15.99081912</MEAN_MOTION><ECCENTRICITY>.0014649</ECCENTRICITY><INCLINATION>51.6242</INCLINATION><RA_OF_ASC_NODE>216.2930</RA_OF_ASC_NODE><ARG_OF_PERICENTER>331.8976</ARG_OF_PERICENTER><MEAN_ANOMALY>28.1241</MEAN_ANOMALY></meanElements><tleParameters><EPHEMERIS_TYPE>0</EPHEMERIS_TYPE><CLASSIFICATION_TYPE>U</CLASSIFICATION_TYPE><NORAD_CAT_ID>55123</NORAD_CAT_ID><ELEMENT_SET_NO>999</ELEMENT_SET_NO><REV_AT_EPOCH>1839</REV_AT_EPOCH><BSTAR>.1568E-2</BSTAR><MEAN_MOTION_DOT>.787702E-2</MEAN_MOTION_DOT><MEAN_MOTION_DDOT>.29408E-3</MEAN_MOTION_DDOT></tleParameters></data></segment></body></omm>
</ndm>
"#;

const MARIO_CSV: &str = "\
OBJECT_NAME,OBJECT_ID,EPOCH,MEAN_MOTION,ECCENTRICITY,INCLINATION,RA_OF_ASC_NODE,ARG_OF_PERICENTER,MEAN_ANOMALY,EPHEMERIS_TYPE,CLASSIFICATION_TYPE,NORAD_CAT_ID,ELEMENT_SET_NO,REV_AT_EPOCH,BSTAR,MEAN_MOTION_DOT,MEAN_MOTION_DDOT
MARIO,1998-067UQ,2023-04-25T10:45:30.642912,15.99081912,.0014649,51.6242,216.2930,331.8976,28.1241,0,U,55123,999,1839,.1568E-2,.787702E-2,.29408E-3
";

const MARIO_KVN: &str = "\
CCSDS_OMM_VERS = 2.0
COMMENT the MARIO element set of the XML, in KVN
CREATION_DATE = 2023-04-25T12:00:00
ORIGINATOR = CELESTRAK
OBJECT_NAME = MARIO
OBJECT_ID = 1998-067UQ
CENTER_NAME = EARTH
REF_FRAME = TEME
TIME_SYSTEM = UTC
MEAN_ELEMENT_THEORY = SGP4
EPOCH = 2023-04-25T10:45:30.642912
MEAN_MOTION = 15.99081912 [rev/day]
ECCENTRICITY = .0014649
INCLINATION = 51.6242 [deg]
RA_OF_ASC_NODE = 216.2930 [deg]
ARG_OF_PERICENTER = 331.8976 [deg]
MEAN_ANOMALY = 28.1241 [deg]
EPHEMERIS_TYPE = 0
CLASSIFICATION_TYPE = U
NORAD_CAT_ID = 55123
ELEMENT_SET_NO = 999
REV_AT_EPOCH = 1839
BSTAR = .1568E-2 [1/ER]
MEAN_MOTION_DOT = .787702E-2 [rev/day**2]
MEAN_MOTION_DDOT = .29408E-3 [rev/day**3]
";

const MARIO_JSON: &str = r#"[{"OBJECT_NAME":"MARIO","OBJECT_ID":"1998-067UQ","EPOCH":"2023-04-25T10:45:30.642912","MEAN_MOTION":15.99081912,"ECCENTRICITY":0.0014649,"INCLINATION":51.6242,"RA_OF_ASC_NODE":216.293,"ARG_OF_PERICENTER":331.8976,"MEAN_ANOMALY":28.1241,"EPHEMERIS_TYPE":0,"CLASSIFICATION_TYPE":"U","NORAD_CAT_ID":55123,"ELEMENT_SET_NO":999,"REV_AT_EPOCH":1839,"BSTAR":0.001568,"MEAN_MOTION_DOT":0.00787702,"MEAN_MOTION_DDOT":0.00029408}]"#;

const VANGUARD_TLE: &str = "\
1 00005U 58002B   20287.20333880 -.00000016  00000-0 -22483-4 0  9998
2 00005  34.2443 225.5254 1845686 162.2516 205.2356 10.84869164218149";

/// The OMM fields python-sgp4 exports for the Vanguard element set.
const VANGUARD_JSON: &str = r#"{"ARG_OF_PERICENTER": 162.2516, "BSTAR": -2.2483e-05, "CENTER_NAME": "EARTH",
"CLASSIFICATION_TYPE": "U", "ECCENTRICITY": 0.1845686, "ELEMENT_SET_NO": 999, "EPHEMERIS_TYPE": 0,
"EPOCH": "2020-10-13T04:52:48.472320", "INCLINATION": 34.2443, "MEAN_ANOMALY": 205.2356,
"MEAN_ELEMENT_THEORY": "SGP4", "MEAN_MOTION": 10.84869164, "MEAN_MOTION_DDOT": 0.0,
"MEAN_MOTION_DOT": -1.6e-07, "NORAD_CAT_ID": 5, "OBJECT_ID": "1958-002B", "OBJECT_NAME": "VANGUARD 1",
"RA_OF_ASC_NODE": 225.5254, "REF_FRAME": "TEME", "REV_AT_EPOCH": 21814, "TIME_SYSTEM": "UTC"}"#;

fn from_tle(text: &str) -> (tle::Tle, MeanElements) {
    let t = tle::parse(text).expect("TLE parses");
    let m = MeanElements::from_tle(&t);
    (t, m)
}

/// The fields of `omm` against the TLE; returns the number of fields compared.
fn assert_same_input(what: &str, omm: &Omm, tle_text: &str) -> usize {
    let (t, a) = from_tle(tle_text);
    let b = omm.to_mean_elements();
    let bits = [
        ("bstar", a.bstar, b.bstar),
        ("ndot", a.ndot, b.ndot),
        ("nddot", a.nddot, b.nddot),
        ("eccentricity", a.eccentricity, b.eccentricity),
        ("arg_perigee", a.arg_perigee, b.arg_perigee),
        ("inclination", a.inclination, b.inclination),
        ("mean_anomaly", a.mean_anomaly, b.mean_anomaly),
        (
            "mean_motion_kozai",
            a.mean_motion_kozai,
            b.mean_motion_kozai,
        ),
        ("raan", a.raan, b.raan),
    ];
    for (name, x, y) in bits {
        assert_eq!(x.to_bits(), y.to_bits(), "{what}: {name} {x:e} vs {y:e}");
    }
    assert_eq!(
        a.epoch_days_since_1950.to_bits(),
        b.epoch_days_since_1950.to_bits(),
        "{what}: epoch {} vs {} days",
        a.epoch_days_since_1950,
        b.epoch_days_since_1950
    );
    assert_eq!(
        omm.norad_cat_id, t.satellite_number,
        "{what}: catalogue number"
    );
    assert_eq!(
        omm.classification, t.classification,
        "{what}: classification"
    );

    // SGP4 from both: within the distance of one microsecond (8 mm) at 0 and 1 day
    let sa = Sgp4::new(&a, GravityModel::Wgs72, OpsMode::Improved).expect("TLE set");
    let sb = Sgp4::new(&b, GravityModel::Wgs72, OpsMode::Improved).expect("OMM set");
    for t_min in [0.0, 1440.0] {
        let (pa, pb) = (sa.propagate(t_min).unwrap(), sb.propagate(t_min).unwrap());
        for k in 0..3 {
            let d = (pa.position_km[k] - pb.position_km[k]).abs();
            assert!(
                d < 8.0e-6,
                "{what}: t={t_min} position[{k}] differs by {d} km"
            );
        }
    }
    bits.len() + 3
}

#[test]
fn every_encoding_of_mario_gives_the_tles_sgp4_input() {
    let mut compared = 0;
    for (what, sets) in [
        ("XML", omm::parse_xml(MARIO_XML)),
        ("CSV", omm::parse_csv(MARIO_CSV)),
        ("KVN", omm::parse_kvn(MARIO_KVN)),
        ("JSON", omm::parse_json(MARIO_JSON)),
    ] {
        let sets = sets.unwrap_or_else(|e| panic!("{what}: {e}"));
        assert_eq!(sets.len(), 1, "{what}: one element set");
        assert_eq!(sets[0].object_name, "MARIO");
        assert_eq!(sets[0].object_id, "1998-067UQ");
        assert_eq!(sets[0].element_set_no, 999);
        assert_eq!(sets[0].rev_at_epoch, 1839);
        compared += assert_same_input(what, &sets[0], MARIO_TLE);
    }
    assert_eq!(compared, 4 * 12, "4 encodings x 12 comparisons");
}

#[test]
fn vanguards_exported_fields_give_the_tles_sgp4_input() {
    let sets = omm::parse_json(VANGUARD_JSON).expect("JSON");
    assert_eq!(sets.len(), 1);
    assert_eq!(sets[0].object_name, "VANGUARD 1");
    assert_eq!(assert_same_input("Vanguard", &sets[0], VANGUARD_TLE), 12);
}

const VERIFICATION: &str = include_str!("data/sgp4/SGP4-VER.TLE");

/// `line[..68]` followed by its recomputed check digit.
fn with_check_digit(line: &str) -> String {
    let payload = &line[..68];
    let sum: u32 = payload
        .chars()
        .map(|c| c.to_digit(10).unwrap_or(u32::from(c == '-')))
        .sum();
    format!("{payload}{}", sum % 10)
}

/// An assumed-decimal TLE field written out: `" -22483-4"` -> `"-0.22483e-4"`.
fn written_out(field: &str) -> String {
    let f = field.trim();
    if f.is_empty() {
        return "0".to_owned();
    }
    let (sign, body) = match f.as_bytes()[0] {
        b'-' => ("-", &f[1..]),
        b'+' => ("", &f[1..]),
        _ => ("", f),
    };
    let cut = body
        .rfind(['+', '-'])
        .filter(|&i| i > 0)
        .unwrap_or(body.len());
    let exponent = if cut == body.len() { "0" } else { &body[cut..] };
    format!("{sign}0.{}e{exponent}", &body[..cut])
}

/// The TLE `l1` / `l2` as a KVN OMM, field by field from the TLE's columns.
fn kvn_of_tle(l1: &str, l2: &str) -> String {
    let yy: u32 = l1[18..20].parse().unwrap();
    let year = if yy < 57 { 2000 + yy } else { 1900 + yy };
    let ddd = &l1[20..23];
    assert_eq!(&l1[23..24], ".", "epoch format");
    let frac: u64 = l1[24..32].parse().unwrap();
    let micros = frac * 864; // 1e-8 day = 864 µs
    let (secs, us) = (micros / 1_000_000, micros % 1_000_000);
    let (h, m, sec) = (secs / 3600, secs / 60 % 60, secs % 60);
    format!(
        "OBJECT_NAME = VALLADO {sat}\n\
         CENTER_NAME = EARTH\nREF_FRAME = TEME\nTIME_SYSTEM = UTC\nMEAN_ELEMENT_THEORY = SGP4\n\
         EPOCH = {year}-{ddd}T{h:02}:{m:02}:{sec:02}.{us:06}\n\
         MEAN_MOTION = {n}\nECCENTRICITY = 0.{e}\nINCLINATION = {i}\nRA_OF_ASC_NODE = {o}\n\
         ARG_OF_PERICENTER = {w}\nMEAN_ANOMALY = {ma}\nCLASSIFICATION_TYPE = {c}\n\
         NORAD_CAT_ID = {sat}\nBSTAR = {b}\nMEAN_MOTION_DOT = {nd}\nMEAN_MOTION_DDOT = {ndd}\n",
        sat = l1[2..7].trim(),
        n = l2[52..63].trim(),
        e = l2[26..33].replace(' ', "0"),
        i = l2[8..16].trim(),
        o = l2[17..25].trim(),
        w = l2[34..42].trim(),
        ma = l2[43..51].trim(),
        c = &l1[7..8],
        b = written_out(&l1[53..61]),
        nd = l1[33..43].trim(),
        ndd = written_out(&l1[44..52]),
    )
}

#[test]
fn every_verification_tle_written_as_an_omm_gives_the_same_sgp4_input() {
    let mut lines = VERIFICATION.lines().map(|l| l.trim_end_matches('\r'));
    let mut compared = 0;
    let mut sets = 0;
    while let Some(l1) = lines.next() {
        if !l1.starts_with("1 ") {
            continue;
        }
        let l2 = lines.next().expect("line 2");
        let (l1, l2) = (with_check_digit(&l1[..69]), with_check_digit(&l2[..69]));
        let kvn = kvn_of_tle(&l1, &l2);
        let omm = omm::parse_kvn(&kvn).unwrap_or_else(|e| panic!("{kvn}\n{e}"));
        assert_eq!(omm.len(), 1);
        let (t, a) = from_tle(&format!("{l1}\n{l2}"));
        let b = omm[0].to_mean_elements();
        for (name, x, y) in [
            ("epoch", a.epoch_days_since_1950, b.epoch_days_since_1950),
            ("bstar", a.bstar, b.bstar),
            ("ndot", a.ndot, b.ndot),
            ("nddot", a.nddot, b.nddot),
            ("eccentricity", a.eccentricity, b.eccentricity),
            ("arg_perigee", a.arg_perigee, b.arg_perigee),
            ("inclination", a.inclination, b.inclination),
            ("mean_anomaly", a.mean_anomaly, b.mean_anomaly),
            (
                "mean_motion_kozai",
                a.mean_motion_kozai,
                b.mean_motion_kozai,
            ),
            ("raan", a.raan, b.raan),
        ] {
            assert_eq!(
                x.to_bits(),
                y.to_bits(),
                "set {}: {name} {x:e} vs {y:e}",
                t.satellite_number
            );
            compared += 1;
        }
        sets += 1;
    }
    assert_eq!(sets, 33, "33 verification sets");
    assert_eq!(compared, 330);
}
