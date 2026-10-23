//! Oracle: the assumed-decimal fields of a TLE (BSTAR and the second derivative
//! of the mean motion, e.g. ` -22483-4` = −0.22483 × 10⁻⁴) decode to the
//! correctly rounded `f64` of the decimal they denote.
//!
//! The reference value is Rust's own decimal parse of the same number written
//! out (`"-0.22483e-4".parse::<f64>()`), which IEEE 754 requires to be the
//! nearest double. A decode that multiplies by a power of ten computed at run
//! time (`powi`, whose precision Rust leaves unspecified) can land one ulp away,
//! and the field feeds SGP4.
//!
//! The inputs are every assumed-decimal field of Vallado's verification set and
//! two CelesTrak element sets.
//!
//! Author: Moroya Sakamoto

use alice_space::tle;

const TLE_FILE: &str = include_str!("data/sgp4/SGP4-VER.TLE");

/// Two CelesTrak element sets (Vanguard 1 and MARIO, every field non-zero).
const CELESTRAK: [&str; 2] = [
    "1 00005U 58002B   20287.20333880 -.00000016  00000-0 -22483-4 0  9998\n\
     2 00005  34.2443 225.5254 1845686 162.2516 205.2356 10.84869164218149",
    "1 55123U 98067UQ  23115.44827133  .00787702  29408-3  15680-2 0  9999\n\
     2 55123  51.6242 216.2930 0014649 331.8976  28.1241 15.99081912 18396",
];

/// `line[..68]` followed by its recomputed check digit (three verification
/// error-case sets carry digits that do not match their lines).
fn with_check_digit(line: &str) -> String {
    let payload = &line[..68];
    let sum: u32 = payload
        .chars()
        .map(|c| c.to_digit(10).unwrap_or(u32::from(c == '-')))
        .sum();
    format!("{payload}{}", sum % 10)
}

/// The decimal an assumed-decimal field denotes, written out: `" -22483-4"` →
/// `"-0.22483e-4"` (an empty field is 0).
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
    let cut = body.rfind(['+', '-']).filter(|&i| i > 0).unwrap_or(body.len());
    let exponent = if cut == body.len() { "0" } else { &body[cut..] };
    format!("{sign}0.{}e{exponent}", &body[..cut])
}

/// (line 1, line 2) of every element set used here.
fn element_sets() -> Vec<(String, String)> {
    let mut out = Vec::new();
    let mut lines = TLE_FILE.lines().map(|l| l.trim_end_matches('\r'));
    while let Some(l1) = lines.next() {
        if l1.starts_with("1 ") {
            let l2 = lines.next().expect("line 2");
            out.push((with_check_digit(&l1[..69]), with_check_digit(&l2[..69])));
        }
    }
    for set in CELESTRAK {
        let mut it = set.lines().map(str::trim);
        out.push((it.next().unwrap().to_owned(), it.next().unwrap().to_owned()));
    }
    out
}

#[test]
fn assumed_decimal_fields_decode_to_the_correctly_rounded_double() {
    let mut compared = 0usize;
    let mut wrong = Vec::new();
    for (l1, l2) in element_sets() {
        let parsed = tle::parse(&format!("{l1}\n{l2}")).expect("parses");
        for (name, field, got) in [
            ("nddot", &l1[44..52], parsed.mean_motion_ddot),
            ("bstar", &l1[53..61], parsed.bstar),
        ] {
            let want: f64 = written_out(field).parse().expect("decimal");
            compared += 1;
            if got.to_bits() != want.to_bits() {
                wrong.extend([format!(
                    "set {} {name} {field:?}: {got:e} (bits {:x}) vs {want:e} (bits {:x})",
                    parsed.satellite_number,
                    got.to_bits(),
                    want.to_bits()
                )]);
            }
        }
    }
    assert_eq!(compared, 70, "35 element sets x 2 fields");
    assert!(
        wrong.is_empty(),
        "{} of {compared} fields differ:\n{}",
        wrong.len(),
        wrong.join("\n")
    );
}
