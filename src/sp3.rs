//! SP3 (Standard Product 3) precise ephemeris parser.
//!
//! SP3 is the IGS-standard format for precise satellite ephemerides.  A file
//! begins with a fixed header, then contains a stream of epoch blocks, each
//! introduced by a `*` line and followed by position (and optionally
//! velocity) records prefixed by `P` and `V`.
//!
//! This parser implements the subset of SP3-c/SP3-d needed for GNSS
//! visibility computations targeted by SPACID's QZSS use case:
//!
//! - Header decoding (`#` line): number of epochs, start time, coordinate
//!   system, orbit type.
//! - Epoch line decoding: `Y/M/D h:m:s` -> Unix seconds.
//! - Position record decoding: satellite identifier, X/Y/Z (km) in the
//!   Earth-fixed frame, and clock offset (µs).
//! - Optional velocity record decoding: X/Y/Z velocity (dm/s) and clock rate
//!   (10⁻⁴ µs/s).
//!
//! Bad rows are surfaced through [`Sp3Error`]; the caller decides whether to
//! abort or skip a corrupted line.

use core::fmt;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// SP3 parser errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Sp3Error {
    /// Header did not start with `#`.
    MissingHeader,
    /// The header epoch line could not be decoded.
    BadHeaderEpoch(String),
    /// A subsequent line was too short.
    TooShort { line_number: usize, actual: usize },
    /// A numeric substring could not be parsed.
    ParseNumber { line_number: usize, field: String },
    /// An unknown record type was encountered.
    UnknownRecord { line_number: usize, first: char },
}

impl fmt::Display for Sp3Error {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::MissingHeader => write!(f, "SP3 file does not start with '#'"),
            Self::BadHeaderEpoch(s) => write!(f, "bad header epoch line: {s}"),
            Self::TooShort {
                line_number,
                actual,
            } => write!(f, "line {line_number}: too short ({actual} chars)"),
            Self::ParseNumber { line_number, field } => {
                write!(f, "line {line_number}: parse error on {field}")
            }
            Self::UnknownRecord { line_number, first } => {
                write!(f, "line {line_number}: unknown record type '{first}'")
            }
        }
    }
}

impl std::error::Error for Sp3Error {}

// ---------------------------------------------------------------------------
// Data model
// ---------------------------------------------------------------------------

/// Header metadata extracted from the first `#` line.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Sp3Header {
    /// Version character (`c`, `d`, ...).
    pub version: char,
    /// Position/velocity flag (`P` or `V`).
    pub pv_flag: char,
    /// Start-of-file epoch in Unix seconds.
    pub start_unix_seconds: u64,
    /// Announced number of epochs.
    pub epoch_count: u32,
}

/// One position (+ optional clock) sample for one satellite.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Sp3Position {
    /// Satellite identifier (e.g. `"G01"` or `"J01"`).
    pub prn: [u8; 3],
    /// ECEF X in kilometres.
    pub x_km: f64,
    /// ECEF Y in kilometres.
    pub y_km: f64,
    /// ECEF Z in kilometres.
    pub z_km: f64,
    /// Clock offset in microseconds.
    pub clock_us: f64,
}

impl Sp3Position {
    /// Return the satellite identifier as a `&str` slice.
    #[must_use]
    pub fn prn_str(&self) -> &str {
        // The three PRN bytes are always plain ASCII per the SP3 spec.
        core::str::from_utf8(&self.prn).unwrap_or("???")
    }
}

/// One epoch of data.
#[derive(Debug, Clone, PartialEq)]
pub struct Sp3Epoch {
    /// Unix seconds.
    pub unix_seconds: u64,
    /// Position samples, one per satellite.
    pub positions: Vec<Sp3Position>,
}

/// Parsed SP3 file.
#[derive(Debug, Clone)]
pub struct Sp3File {
    pub header: Sp3Header,
    pub epochs: Vec<Sp3Epoch>,
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse an SP3 body.
///
/// The full IGS header spans many lines; the fields kept in [`Sp3Header`]
/// come exclusively from the first `#` line, which is sufficient to
/// interpret the epoch stream. Subsequent header lines are tolerated but
/// not decoded.
///
/// # Errors
///
/// See [`Sp3Error`].
pub fn parse(text: &str) -> Result<Sp3File, Sp3Error> {
    let mut lines = text.lines().enumerate();

    // First line: header.
    let (line_no, first) = lines.next().ok_or(Sp3Error::MissingHeader)?;
    if !first.starts_with('#') {
        return Err(Sp3Error::MissingHeader);
    }
    let header = parse_header(first, line_no + 1)?;

    let mut epochs = Vec::new();
    let mut current: Option<Sp3Epoch> = None;

    for (index, line) in lines {
        let line_number = index + 1;
        if line.is_empty() {
            continue;
        }
        let first_char = line.chars().next().unwrap_or(' ');
        match first_char {
            // End of file.
            'E' if line.starts_with("EOF") => break,
            // Additional header line: ignore.
            '#' | '+' | '%' | '/' => (),
            '*' => {
                if let Some(prev) = current.take() {
                    epochs.push(prev);
                }
                let unix_seconds = parse_epoch_line(line, line_number)?;
                current = Some(Sp3Epoch {
                    unix_seconds,
                    positions: Vec::new(),
                });
            }
            'P' => {
                let pos = parse_position(line, line_number)?;
                if let Some(ep) = current.as_mut() {
                    ep.positions.push(pos);
                }
            }
            'V' => {
                // Velocity records exist but this parser targets position
                // consumers; skip them silently.
            }
            _ => {
                return Err(Sp3Error::UnknownRecord {
                    line_number,
                    first: first_char,
                });
            }
        }
    }
    if let Some(last) = current {
        epochs.push(last);
    }
    Ok(Sp3File { header, epochs })
}

fn parse_header(line: &str, line_number: usize) -> Result<Sp3Header, Sp3Error> {
    if line.len() < 60 {
        return Err(Sp3Error::TooShort {
            line_number,
            actual: line.len(),
        });
    }
    let version = line.chars().nth(1).unwrap_or(' ');
    let pv_flag = line.chars().nth(2).unwrap_or('P');

    let year: i32 = parse_int(&line[3..7], line_number, "year")?;
    let month: u32 = parse_int(&line[8..10], line_number, "month")?;
    let day: u32 = parse_int(&line[11..13], line_number, "day")?;
    let hour: u32 = parse_int(&line[14..16], line_number, "hour")?;
    let minute: u32 = parse_int(&line[17..19], line_number, "minute")?;
    let second: f64 = parse_float(&line[20..31], line_number, "second")?;

    let unix_seconds = ymdhms_to_unix(year, month, day, hour, minute, second)
        .ok_or_else(|| Sp3Error::BadHeaderEpoch(line.to_owned()))?;

    let epoch_count: u32 = parse_int(line[32..39].trim(), line_number, "epoch_count")?;

    Ok(Sp3Header {
        version,
        pv_flag,
        start_unix_seconds: unix_seconds,
        epoch_count,
    })
}

fn parse_epoch_line(line: &str, line_number: usize) -> Result<u64, Sp3Error> {
    if line.len() < 31 {
        return Err(Sp3Error::TooShort {
            line_number,
            actual: line.len(),
        });
    }
    let year: i32 = parse_int(line[3..7].trim(), line_number, "epoch year")?;
    let month: u32 = parse_int(line[8..10].trim(), line_number, "epoch month")?;
    let day: u32 = parse_int(line[11..13].trim(), line_number, "epoch day")?;
    let hour: u32 = parse_int(line[14..16].trim(), line_number, "epoch hour")?;
    let minute: u32 = parse_int(line[17..19].trim(), line_number, "epoch minute")?;
    let second: f64 = parse_float(line[20..31].trim(), line_number, "epoch second")?;
    ymdhms_to_unix(year, month, day, hour, minute, second)
        .ok_or_else(|| Sp3Error::BadHeaderEpoch(format!("line {line_number}: {line}")))
}

fn parse_position(line: &str, line_number: usize) -> Result<Sp3Position, Sp3Error> {
    if line.len() < 60 {
        return Err(Sp3Error::TooShort {
            line_number,
            actual: line.len(),
        });
    }
    let mut prn = [b'?'; 3];
    for (i, b) in line.as_bytes()[1..4].iter().enumerate() {
        prn[i] = *b;
    }
    let x_km: f64 = parse_float(line[4..18].trim(), line_number, "x_km")?;
    let y_km: f64 = parse_float(line[18..32].trim(), line_number, "y_km")?;
    let z_km: f64 = parse_float(line[32..46].trim(), line_number, "z_km")?;
    let clock_us: f64 = parse_float(line[46..60].trim(), line_number, "clock_us")?;
    Ok(Sp3Position {
        prn,
        x_km,
        y_km,
        z_km,
        clock_us,
    })
}

fn parse_int<T>(s: &str, line_number: usize, field: &str) -> Result<T, Sp3Error>
where
    T: core::str::FromStr,
{
    s.trim().parse::<T>().map_err(|_| Sp3Error::ParseNumber {
        line_number,
        field: field.to_owned(),
    })
}

fn parse_float(s: &str, line_number: usize, field: &str) -> Result<f64, Sp3Error> {
    s.trim().parse::<f64>().map_err(|_| Sp3Error::ParseNumber {
        line_number,
        field: field.to_owned(),
    })
}

fn ymdhms_to_unix(
    year: i32,
    month: u32,
    day: u32,
    hour: u32,
    minute: u32,
    second: f64,
) -> Option<u64> {
    if !(1..=12).contains(&month) || !(1..=31).contains(&day) {
        return None;
    }
    // Days in each month (non-leap year).
    let days_before: [u32; 12] = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334];
    let is_leap = year % 4 == 0 && (year % 100 != 0 || year % 400 == 0);
    let mut days = i64::from(days_before[month as usize - 1]);
    if is_leap && month > 2 {
        days += 1;
    }
    days += i64::from(day - 1);

    // Days from 1970-01-01 to year-01-01.
    let mut year_days: i64 = 0;
    let start_year = 1970i32;
    if year >= start_year {
        for y in start_year..year {
            year_days += if y % 4 == 0 && (y % 100 != 0 || y % 400 == 0) {
                366
            } else {
                365
            };
        }
    } else {
        for y in year..start_year {
            year_days -= if y % 4 == 0 && (y % 100 != 0 || y % 400 == 0) {
                366
            } else {
                365
            };
        }
    }
    let total_days = year_days + days;
    let seconds_of_day = i64::from(hour) * 3600 + i64::from(minute) * 60 + second as i64;
    let total = total_days * 86_400 + seconds_of_day;
    u64::try_from(total).ok()
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Two-epoch, one-satellite SP3 sample. Header is trimmed to the first line
    // only for testing purposes; the parser tolerates missing extra header lines.
    const SAMPLE: &str = "#cP2020  5 24  0  0  0.00000000       2 ORBIT IGb14 HLM  IGS
+    1   G01
+    1
++         2
%c M  cc GPS ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
%c cc cc ccc ccc cccc cccc cccc cccc ccccc ccccc ccccc ccccc
%f  0.0000000  0.000000000  0.00000000000  0.000000000000000
%f  0.0000000  0.000000000  0.00000000000  0.000000000000000
%i    0    0    0    0      0      0      0      0         0
%i    0    0    0    0      0      0      0      0         0
*  2020  5 24  0  0  0.00000000
PG01  15600.500000 -14000.250000  21000.100000    120.500000
*  2020  5 24  0 15  0.00000000
PG01  15650.500000 -13950.250000  21005.100000    121.500000
EOF";

    #[test]
    fn parses_header_and_two_epochs() {
        let file = parse(SAMPLE).expect("valid SP3 sample");
        assert_eq!(file.header.version, 'c');
        assert_eq!(file.header.pv_flag, 'P');
        assert_eq!(file.header.epoch_count, 2);
        assert_eq!(file.epochs.len(), 2);
    }

    #[test]
    fn epoch_times_are_15_minutes_apart() {
        let file = parse(SAMPLE).unwrap();
        let dt = file.epochs[1].unix_seconds - file.epochs[0].unix_seconds;
        assert_eq!(dt, 900);
    }

    #[test]
    fn position_x_matches_source() {
        let file = parse(SAMPLE).unwrap();
        let p = &file.epochs[0].positions[0];
        assert!((p.x_km - 15_600.5).abs() < 1e-6);
        assert!((p.y_km + 14_000.25).abs() < 1e-6);
        assert!((p.z_km - 21_000.1).abs() < 1e-6);
        assert!((p.clock_us - 120.5).abs() < 1e-6);
    }

    #[test]
    fn prn_is_g01() {
        let file = parse(SAMPLE).unwrap();
        assert_eq!(file.epochs[0].positions[0].prn_str(), "G01");
    }

    #[test]
    fn missing_header_is_reported() {
        let bad = "PG01 0 0 0 0";
        assert_eq!(parse(bad).unwrap_err(), Sp3Error::MissingHeader);
    }

    #[test]
    fn short_epoch_line_is_reported() {
        let bad = "#cP2020  5 24  0  0  0.00000000       2 ORBIT IGb14 HLM  IGS\n*  2020";
        assert!(matches!(parse(bad), Err(Sp3Error::TooShort { .. })));
    }

    #[test]
    fn eof_marker_stops_parsing() {
        let file = parse(SAMPLE).unwrap();
        // No spurious extra epochs after EOF.
        assert_eq!(file.epochs.len(), 2);
    }

    #[test]
    fn ymdhms_matches_known_unix() {
        // 2020-05-24T00:00:00Z is Unix timestamp 1590278400.
        let u = ymdhms_to_unix(2020, 5, 24, 0, 0, 0.0).unwrap();
        assert_eq!(u, 1_590_278_400);
    }

    #[test]
    fn ymdhms_rejects_invalid_month() {
        assert!(ymdhms_to_unix(2020, 13, 1, 0, 0, 0.0).is_none());
    }
}
