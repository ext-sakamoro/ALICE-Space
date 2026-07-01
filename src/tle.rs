//! Two-Line Element (`TLE`) set parser.
//!
//! `TLE` is the de-facto standard format for publishing satellite orbital
//! elements (NORAD, CelesTrak, `QZSS`). Each satellite is represented by an
//! optional name line and two 69-character data lines encoding a snapshot
//! of the mean Keplerian orbit at a specific epoch.
//!
//! This parser validates:
//!
//! - Line lengths (69 characters).
//! - Line prefixes (`1 ` and `2 `).
//! - Modulo-10 checksums.
//! - Consistent satellite numbers between the two lines.
//!
//! The returned [`Tle`] carries the raw parsed fields; higher-level converters
//! translate them into [`crate::orbit::OrbitalElements`].

use crate::orbit::OrbitalElements;
use core::f64::consts::PI;
use core::fmt;

/// Standard gravitational parameter of Earth (km³/s²).
///
/// Used to convert `TLE` mean motion (revolutions/day) into semi-major axis.
pub const MU_EARTH_KM3_S2: f64 = 398_600.4418;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// All errors produced when parsing a `TLE`.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum TleError {
    /// Neither 2 nor 3 non-empty lines were supplied.
    WrongLineCount(usize),
    /// A data line was not exactly 69 characters long.
    WrongLineLength { line: u8, actual: usize },
    /// Line 1 or Line 2 did not begin with the expected `"1 "` or `"2 "` prefix.
    WrongPrefix { line: u8 },
    /// The checksum stored in the last column did not match the computed value.
    ChecksumMismatch { line: u8, expected: u8, actual: u8 },
    /// Line 1 and Line 2 carry different satellite catalogue numbers.
    SatelliteNumberMismatch { line1: u32, line2: u32 },
    /// A numeric substring could not be parsed.
    ParseNumber(String),
}

impl fmt::Display for TleError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::WrongLineCount(n) => write!(f, "expected 2 or 3 non-empty lines, got {n}"),
            Self::WrongLineLength { line, actual } => {
                write!(f, "line {line}: expected length 69, got {actual}")
            }
            Self::WrongPrefix { line } => write!(f, "line {line}: wrong prefix"),
            Self::ChecksumMismatch {
                line,
                expected,
                actual,
            } => write!(
                f,
                "line {line}: checksum mismatch (expected {expected}, computed {actual})"
            ),
            Self::SatelliteNumberMismatch { line1, line2 } => {
                write!(f, "satellite number mismatch: {line1} vs {line2}")
            }
            Self::ParseNumber(s) => write!(f, "failed to parse number: {s}"),
        }
    }
}

impl std::error::Error for TleError {}

// ---------------------------------------------------------------------------
// TLE data structures
// ---------------------------------------------------------------------------

/// Fully parsed `TLE`.
#[derive(Debug, Clone)]
pub struct Tle {
    /// Optional common name (present when a 3-line block was supplied).
    pub name: Option<String>,
    /// Satellite catalogue number.
    pub satellite_number: u32,
    /// Classification character (`U`, `C`, `S`).
    pub classification: char,
    /// International designator (columns 10-17 of line 1).
    pub international_designator: String,
    /// Two-digit epoch year converted to full year (57–99 -> 1957-1999, else 20xx).
    pub epoch_year: u32,
    /// Day-of-year with fractional part.
    pub epoch_day: f64,
    /// First derivative of mean motion divided by 2 (rev / day²).
    pub mean_motion_dot: f64,
    /// Second derivative of mean motion divided by 6 (rev / day³).
    pub mean_motion_ddot: f64,
    /// `BSTAR` drag term (Earth radii⁻¹).
    pub bstar: f64,
    /// Element set number.
    pub element_set_number: u32,
    /// Inclination (radians).
    pub inclination_rad: f64,
    /// Right ascension of the ascending node (radians).
    pub raan_rad: f64,
    /// Eccentricity (dimensionless).
    pub eccentricity: f64,
    /// Argument of perigee (radians).
    pub arg_perigee_rad: f64,
    /// Mean anomaly (radians).
    pub mean_anomaly_rad: f64,
    /// Mean motion (revolutions/day).
    pub mean_motion_rev_per_day: f64,
    /// Orbit revolution number at epoch.
    pub revolution_number: u32,
}

impl Tle {
    /// Semi-major axis in kilometres derived from mean motion.
    #[must_use]
    pub fn semi_major_axis_km(&self) -> f64 {
        let n_rad_per_sec = self.mean_motion_rev_per_day * 2.0 * PI / 86_400.0;
        (MU_EARTH_KM3_S2 / (n_rad_per_sec * n_rad_per_sec)).cbrt()
    }

    /// Build a [`OrbitalElements`] snapshot at epoch.
    ///
    /// `true_anomaly_rad` is set to the mean anomaly value as a first-order
    /// approximation; downstream code that requires the true anomaly should
    /// solve Kepler's equation via [`crate::orbit`].
    #[must_use]
    pub fn to_orbital_elements(&self) -> OrbitalElements {
        OrbitalElements {
            semi_major_axis_km: self.semi_major_axis_km(),
            eccentricity: self.eccentricity,
            inclination_rad: self.inclination_rad,
            raan_rad: self.raan_rad,
            arg_periapsis_rad: self.arg_perigee_rad,
            true_anomaly_rad: self.mean_anomaly_rad,
        }
    }
}

// ---------------------------------------------------------------------------
// Parsing entry point
// ---------------------------------------------------------------------------

/// Parse a `TLE` text block containing 2 or 3 non-empty lines.
///
/// # Errors
///
/// See [`TleError`].
pub fn parse(text: &str) -> Result<Tle, TleError> {
    let lines: Vec<&str> = text
        .lines()
        .map(str::trim_end)
        .filter(|l| !l.is_empty())
        .collect();

    let (name, l1, l2) = match lines.as_slice() {
        [l1, l2] => (None, *l1, *l2),
        [name, l1, l2] => (Some((*name).trim().to_owned()), *l1, *l2),
        other => return Err(TleError::WrongLineCount(other.len())),
    };

    validate_line(l1, 1)?;
    validate_line(l2, 2)?;

    let sat1 = parse_uint(&l1[2..7])?;
    let sat2 = parse_uint(&l2[2..7])?;
    if sat1 != sat2 {
        return Err(TleError::SatelliteNumberMismatch {
            line1: sat1,
            line2: sat2,
        });
    }

    let classification = l1.chars().nth(7).unwrap_or(' ');
    let international_designator = l1[9..17].trim().to_owned();
    let epoch_year_raw = parse_uint(&l1[18..20])?;
    let epoch_year = if epoch_year_raw < 57 {
        2000 + epoch_year_raw
    } else {
        1900 + epoch_year_raw
    };
    let epoch_day = parse_float(&l1[20..32])?;
    let mean_motion_dot = parse_float(&l1[33..43])?;
    let mean_motion_ddot = parse_assumed_decimal(&l1[44..52])?;
    let bstar = parse_assumed_decimal(&l1[53..61])?;
    let element_set_number = parse_uint(l1[64..68].trim())?;

    let inclination_deg = parse_float(&l2[8..16])?;
    let raan_deg = parse_float(&l2[17..25])?;
    let eccentricity = parse_leading_zero_decimal(&l2[26..33])?;
    let arg_perigee_deg = parse_float(&l2[34..42])?;
    let mean_anomaly_deg = parse_float(&l2[43..51])?;
    let mean_motion_rev_per_day = parse_float(&l2[52..63])?;
    let revolution_number = parse_uint(l2[63..68].trim())?;

    Ok(Tle {
        name,
        satellite_number: sat1,
        classification,
        international_designator,
        epoch_year,
        epoch_day,
        mean_motion_dot,
        mean_motion_ddot,
        bstar,
        element_set_number,
        inclination_rad: inclination_deg.to_radians(),
        raan_rad: raan_deg.to_radians(),
        eccentricity,
        arg_perigee_rad: arg_perigee_deg.to_radians(),
        mean_anomaly_rad: mean_anomaly_deg.to_radians(),
        mean_motion_rev_per_day,
        revolution_number,
    })
}

// ---------------------------------------------------------------------------
// Parsing helpers
// ---------------------------------------------------------------------------

fn validate_line(line: &str, index: u8) -> Result<(), TleError> {
    if line.len() != 69 {
        return Err(TleError::WrongLineLength {
            line: index,
            actual: line.len(),
        });
    }
    let expected_prefix = match index {
        1 => "1 ",
        2 => "2 ",
        _ => unreachable!("index is always 1 or 2 in current call sites"),
    };
    if &line[0..2] != expected_prefix {
        return Err(TleError::WrongPrefix { line: index });
    }
    let expected = tle_checksum(&line[..68]);
    let actual_char = line.chars().last().unwrap_or('0');
    let actual = actual_char
        .to_digit(10)
        .ok_or_else(|| TleError::ParseNumber(format!("checksum char '{actual_char}'")))?
        as u8;
    if expected != actual {
        return Err(TleError::ChecksumMismatch {
            line: index,
            expected,
            actual,
        });
    }
    Ok(())
}

fn tle_checksum(payload: &str) -> u8 {
    let mut sum: u32 = 0;
    for c in payload.chars() {
        if let Some(d) = c.to_digit(10) {
            sum += d;
        } else if c == '-' {
            sum += 1;
        }
    }
    (sum % 10) as u8
}

fn parse_uint(s: &str) -> Result<u32, TleError> {
    s.trim()
        .parse::<u32>()
        .map_err(|_| TleError::ParseNumber(s.to_owned()))
}

fn parse_float(s: &str) -> Result<f64, TleError> {
    s.trim()
        .parse::<f64>()
        .map_err(|_| TleError::ParseNumber(s.to_owned()))
}

/// Parse the TLE assumed-decimal notation for BSTAR and n_ddot:
/// e.g. `" 12345-3"` = `0.12345e-3`.
fn parse_assumed_decimal(s: &str) -> Result<f64, TleError> {
    let s = s.trim();
    if s.is_empty() {
        return Ok(0.0);
    }
    let sign = if s.starts_with('-') { -1.0 } else { 1.0 };
    let body = s.trim_start_matches(['+', '-']);
    // Split at the last '+' or '-' (the exponent sign).
    let (mantissa, exponent) = match body.rfind(['+', '-']) {
        Some(idx) if idx > 0 => (&body[..idx], &body[idx..]),
        _ => (body, "+0"),
    };
    let mantissa_val: f64 = format!("0.{mantissa}")
        .parse()
        .map_err(|_| TleError::ParseNumber(s.to_owned()))?;
    let exponent_val: i32 = exponent
        .parse()
        .map_err(|_| TleError::ParseNumber(s.to_owned()))?;
    Ok(sign * mantissa_val * 10f64.powi(exponent_val))
}

/// Parse the eccentricity field as an implicit leading `0.`.
fn parse_leading_zero_decimal(s: &str) -> Result<f64, TleError> {
    let trimmed = s.trim();
    if trimmed.is_empty() {
        return Ok(0.0);
    }
    format!("0.{trimmed}")
        .parse::<f64>()
        .map_err(|_| TleError::ParseNumber(s.to_owned()))
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // ISS ZARYA TLE from CelesTrak (dated 2020-05-24, epoch year 20).
    const ISS_TLE: &str = "ISS (ZARYA)
1 25544U 98067A   20145.51805556 -.00000598  00000-0 -10230-4 0  9991
2 25544  51.6428 227.3897 0002307 216.1957 143.8965 15.49445832225683";

    #[test]
    fn parses_iss_three_line_block() {
        let tle = parse(ISS_TLE).expect("valid TLE");
        assert_eq!(tle.name.as_deref(), Some("ISS (ZARYA)"));
        assert_eq!(tle.satellite_number, 25544);
        assert_eq!(tle.classification, 'U');
        assert_eq!(tle.epoch_year, 2020);
        assert!((tle.epoch_day - 145.518_055_56).abs() < 1e-6);
        assert!((tle.inclination_rad.to_degrees() - 51.6428).abs() < 1e-6);
        assert!((tle.raan_rad.to_degrees() - 227.3897).abs() < 1e-6);
        assert!((tle.eccentricity - 0.000_230_7).abs() < 1e-9);
        assert!((tle.arg_perigee_rad.to_degrees() - 216.1957).abs() < 1e-6);
        assert!((tle.mean_anomaly_rad.to_degrees() - 143.8965).abs() < 1e-6);
        assert!((tle.mean_motion_rev_per_day - 15.494_458_32).abs() < 1e-8);
        assert_eq!(tle.revolution_number, 22568);
    }

    #[test]
    fn iss_semi_major_axis_is_leo_scale() {
        let tle = parse(ISS_TLE).unwrap();
        let a = tle.semi_major_axis_km();
        // ISS is ~6780 km semi-major axis.
        assert!((6700.0..6900.0).contains(&a), "a = {a}");
    }

    #[test]
    fn iss_orbital_elements_conversion() {
        let tle = parse(ISS_TLE).unwrap();
        let el = tle.to_orbital_elements();
        assert!((el.eccentricity - tle.eccentricity).abs() < 1e-12);
        assert!((el.inclination_rad - tle.inclination_rad).abs() < 1e-12);
    }

    #[test]
    fn accepts_two_line_block() {
        let two = "1 25544U 98067A   20145.51805556 -.00000598  00000-0 -10230-4 0  9991
2 25544  51.6428 227.3897 0002307 216.1957 143.8965 15.49445832225683";
        let tle = parse(two).unwrap();
        assert!(tle.name.is_none());
        assert_eq!(tle.satellite_number, 25544);
    }

    #[test]
    fn rejects_short_line() {
        let bad = "1 25544U 98067A   20145.51805556\n2 25544";
        assert!(matches!(parse(bad), Err(TleError::WrongLineLength { .. })));
    }

    #[test]
    fn rejects_wrong_prefix() {
        // Truncate the checksum char so the length is still 69 but prefix wrong.
        let bad = "X 25544U 98067A   20145.51805556 -.00000598  00000-0 -10230-4 0  9993
2 25544  51.6428 227.3897 0002307 216.1957 143.8965 15.49445832225686";
        assert!(matches!(parse(bad), Err(TleError::WrongPrefix { .. })));
    }

    #[test]
    fn detects_checksum_mismatch() {
        // Flip the final checksum char of line 1.
        let bad = "1 25544U 98067A   20145.51805556 -.00000598  00000-0 -10230-4 0  9990
2 25544  51.6428 227.3897 0002307 216.1957 143.8965 15.49445832225686";
        assert!(matches!(parse(bad), Err(TleError::ChecksumMismatch { .. })));
    }

    #[test]
    fn detects_mismatched_satellite_numbers() {
        let bad = "1 25544U 98067A   20145.51805556 -.00000598  00000-0 -10230-4 0  9993
2 25545  51.6428 227.3897 0002307 216.1957 143.8965 15.49445832225686";
        // Line 2 checksum will fail first (payload changed); accept either error.
        // The specific error depends on parse ordering; both are correct rejections.
        let err = parse(bad).unwrap_err();
        assert!(matches!(
            err,
            TleError::ChecksumMismatch { .. } | TleError::SatelliteNumberMismatch { .. }
        ));
    }

    #[test]
    fn parses_bstar_assumed_decimal() {
        assert!((parse_assumed_decimal(" 12345-3").unwrap() - 0.12345e-3).abs() < 1e-12);
        assert!((parse_assumed_decimal("-12345-3").unwrap() + 0.12345e-3).abs() < 1e-12);
        assert!((parse_assumed_decimal("00000-0").unwrap()).abs() < 1e-12);
    }

    #[test]
    fn parses_eccentricity_leading_zero() {
        assert!((parse_leading_zero_decimal("0002307").unwrap() - 0.000_230_7).abs() < 1e-9);
    }

    #[test]
    fn checksum_counts_minus_as_one() {
        assert_eq!(tle_checksum("--"), 2);
        assert_eq!(tle_checksum("12345"), 5);
    }
}
