//! Minimal RINEX 3.x observation file parser.
//!
//! RINEX (Receiver Independent Exchange Format) is the standard container
//! for `GNSS` observation data. This parser implements a subset sufficient
//! for SPACID's authenticity workflow:
//!
//! - Header parsing (`RINEX VERSION / TYPE`, `SYS / # / OBS TYPES`,
//!   `END OF HEADER`).
//! - Epoch header lines starting with `>`.
//! - Per-satellite observation records (satellite id + observation values).
//!
//! Missing fields (e.g. LLI, signal strength indicators) are tolerated and
//! ignored — they do not affect visibility or spoofing checks.

use core::fmt;
use std::collections::BTreeMap;

// ---------------------------------------------------------------------------
// Errors
// ---------------------------------------------------------------------------

/// RINEX parser errors.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum RinexError {
    /// Header was truncated before `END OF HEADER`.
    UnterminatedHeader,
    /// Version identifier could not be decoded.
    BadVersion(String),
    /// A numeric substring could not be parsed.
    ParseNumber { line_number: usize, field: String },
    /// An epoch header line was malformed.
    BadEpoch { line_number: usize, actual: String },
    /// An observation line was too short.
    TooShort { line_number: usize, actual: usize },
}

impl fmt::Display for RinexError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::UnterminatedHeader => write!(f, "RINEX header not terminated"),
            Self::BadVersion(s) => write!(f, "bad version line: {s}"),
            Self::ParseNumber { line_number, field } => {
                write!(f, "line {line_number}: parse error on {field}")
            }
            Self::BadEpoch {
                line_number,
                actual,
            } => {
                write!(f, "line {line_number}: bad epoch header: {actual}")
            }
            Self::TooShort {
                line_number,
                actual,
            } => write!(f, "line {line_number}: too short ({actual} chars)"),
        }
    }
}

impl std::error::Error for RinexError {}

// ---------------------------------------------------------------------------
// Data model
// ---------------------------------------------------------------------------

/// Header metadata retained after parsing.
#[derive(Debug, Clone, PartialEq)]
pub struct RinexHeader {
    /// RINEX version (e.g. `3.05`).
    pub version: f64,
    /// File type code (e.g. `O` for observation).
    pub file_type: char,
    /// Observation type strings per constellation letter (e.g. `G` -> `["C1C", "L1C"]`).
    pub obs_types: BTreeMap<char, Vec<String>>,
}

/// One epoch header + list of satellite observations.
#[derive(Debug, Clone, PartialEq)]
pub struct RinexEpoch {
    /// Year (4-digit).
    pub year: i32,
    pub month: u32,
    pub day: u32,
    pub hour: u32,
    pub minute: u32,
    pub second: f64,
    /// Epoch flag (typically `0` for normal).
    pub flag: u8,
    /// Per-satellite observations.
    pub observations: Vec<RinexObs>,
}

/// Per-satellite observation record.
#[derive(Debug, Clone, PartialEq)]
pub struct RinexObs {
    /// Satellite identifier, e.g. `"G01"` or `"J07"`.
    pub prn: String,
    /// Observation values, one per header-declared type; `None` when missing.
    pub values: Vec<Option<f64>>,
}

/// Fully parsed RINEX observation file.
#[derive(Debug, Clone)]
pub struct RinexFile {
    pub header: RinexHeader,
    pub epochs: Vec<RinexEpoch>,
}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a RINEX 3 observation file body.
///
/// # Errors
///
/// See [`RinexError`].
pub fn parse(text: &str) -> Result<RinexFile, RinexError> {
    let mut lines = text.lines().enumerate();
    let (header, header_last_line) = parse_header(&mut lines)?;

    let mut epochs: Vec<RinexEpoch> = Vec::new();
    let mut current_epoch: Option<(RinexEpoch, u32)> = None;

    for (index, line) in lines {
        let line_number = index + 1;
        let _ = header_last_line;
        if line.is_empty() {
            continue;
        }
        if line.starts_with('>') {
            // Push previously accumulated epoch (if any).
            if let Some((ep, _)) = current_epoch.take() {
                epochs.push(ep);
            }
            let (ep, expected_count) = parse_epoch_header(line, line_number)?;
            current_epoch = Some((ep, expected_count));
        } else if let Some((ep, remaining)) = current_epoch.as_mut() {
            if *remaining == 0 {
                continue;
            }
            let letter = line.chars().next().unwrap_or('?');
            let obs_types_len = header.obs_types.get(&letter).map_or(0, std::vec::Vec::len);
            let obs = parse_obs_line(line, obs_types_len, line_number)?;
            ep.observations.push(obs);
            *remaining -= 1;
        }
    }
    if let Some((ep, _)) = current_epoch {
        epochs.push(ep);
    }
    Ok(RinexFile { header, epochs })
}

fn parse_header(
    lines: &mut core::iter::Enumerate<core::str::Lines<'_>>,
) -> Result<(RinexHeader, usize), RinexError> {
    let mut version = 0.0_f64;
    let mut file_type = 'O';
    let mut obs_types: BTreeMap<char, Vec<String>> = BTreeMap::new();
    let mut pending_letter: Option<char> = None;
    let mut pending_remaining: usize = 0;
    let mut last_line_no = 0;

    for (index, line) in lines.by_ref() {
        let line_number = index + 1;
        last_line_no = line_number.max(last_line_no);
        if line.contains("RINEX VERSION / TYPE") {
            let version_str = line.get(0..9).unwrap_or("").trim();
            version = version_str
                .parse()
                .map_err(|_| RinexError::BadVersion(line.to_owned()))?;
            file_type = line.chars().nth(20).unwrap_or('O');
        } else if line.contains("SYS / # / OBS TYPES") {
            let letter = line.chars().next().unwrap_or('?');
            let count_str = line.get(3..6).unwrap_or("").trim();
            let count: usize = count_str.parse().map_err(|_| RinexError::ParseNumber {
                line_number,
                field: "obs type count".into(),
            })?;
            let types = collect_obs_types(&line[7..]);
            if types.len() >= count {
                obs_types.insert(letter, types.into_iter().take(count).collect());
                pending_letter = None;
                pending_remaining = 0;
            } else {
                pending_letter = Some(letter);
                pending_remaining = count - types.len();
                obs_types.insert(letter, types);
            }
        } else if let Some(letter) = pending_letter {
            // Continuation line for a previously started SYS entry.
            let types = collect_obs_types(line);
            let bucket = obs_types.entry(letter).or_default();
            for t in types.iter().take(pending_remaining) {
                bucket.push(t.clone());
            }
            if types.len() >= pending_remaining {
                pending_letter = None;
                pending_remaining = 0;
            } else {
                pending_remaining -= types.len();
            }
        }
        if line.contains("END OF HEADER") {
            return Ok((
                RinexHeader {
                    version,
                    file_type,
                    obs_types,
                },
                last_line_no,
            ));
        }
    }
    Err(RinexError::UnterminatedHeader)
}

fn collect_obs_types(s: &str) -> Vec<String> {
    s.split_ascii_whitespace()
        .filter(|t| !t.is_empty())
        .map(str::to_owned)
        .collect()
}

fn parse_epoch_header(line: &str, line_number: usize) -> Result<(RinexEpoch, u32), RinexError> {
    // Format: > YYYY MM DD HH MM SS.SSSSSSS  FLAG  N_SATS
    let parts: Vec<&str> = line[1..].split_ascii_whitespace().collect();
    if parts.len() < 8 {
        return Err(RinexError::BadEpoch {
            line_number,
            actual: line.to_owned(),
        });
    }
    let year: i32 = parts[0].parse().map_err(|_| RinexError::ParseNumber {
        line_number,
        field: "epoch year".into(),
    })?;
    let month: u32 = parts[1].parse().map_err(|_| RinexError::ParseNumber {
        line_number,
        field: "epoch month".into(),
    })?;
    let day: u32 = parts[2].parse().map_err(|_| RinexError::ParseNumber {
        line_number,
        field: "epoch day".into(),
    })?;
    let hour: u32 = parts[3].parse().map_err(|_| RinexError::ParseNumber {
        line_number,
        field: "epoch hour".into(),
    })?;
    let minute: u32 = parts[4].parse().map_err(|_| RinexError::ParseNumber {
        line_number,
        field: "epoch minute".into(),
    })?;
    let second: f64 = parts[5].parse().map_err(|_| RinexError::ParseNumber {
        line_number,
        field: "epoch second".into(),
    })?;
    let flag: u8 = parts[6].parse().map_err(|_| RinexError::ParseNumber {
        line_number,
        field: "epoch flag".into(),
    })?;
    let count: u32 = parts[7].parse().map_err(|_| RinexError::ParseNumber {
        line_number,
        field: "epoch sat count".into(),
    })?;
    Ok((
        RinexEpoch {
            year,
            month,
            day,
            hour,
            minute,
            second,
            flag,
            observations: Vec::with_capacity(count as usize),
        },
        count,
    ))
}

fn parse_obs_line(
    line: &str,
    obs_types_len: usize,
    line_number: usize,
) -> Result<RinexObs, RinexError> {
    if line.len() < 3 {
        return Err(RinexError::TooShort {
            line_number,
            actual: line.len(),
        });
    }
    let prn = line[..3].to_owned();
    let rest = &line[3..];
    let mut values = Vec::with_capacity(obs_types_len);
    for i in 0..obs_types_len {
        let start = i * 16;
        let end = (start + 14).min(rest.len());
        let field = if start >= rest.len() {
            ""
        } else {
            rest[start..end].trim()
        };
        if field.is_empty() {
            values.push(None);
        } else {
            let v: f64 = field.parse().map_err(|_| RinexError::ParseNumber {
                line_number,
                field: format!("observation #{}", i + 1),
            })?;
            values.push(Some(v));
        }
    }
    Ok(RinexObs { prn, values })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    // Minimal but valid RINEX 3.05 observation snippet with one epoch and
    // two satellites (GPS and QZSS).
    const SAMPLE: &str =
        "     3.05           O                   M                   RINEX VERSION / TYPE
G    2 C1C L1C                                              SYS / # / OBS TYPES
J    2 C1C L1C                                              SYS / # / OBS TYPES
                                                            END OF HEADER
> 2025 06 24 12 00  0.0000000  0  2
G01  21005432.123    120345678.9
J01  38500123.456    215678901.2";

    #[test]
    fn header_captures_version_and_file_type() {
        let file = parse(SAMPLE).unwrap();
        assert!((file.header.version - 3.05).abs() < 1e-9);
        assert_eq!(file.header.file_type, 'O');
    }

    #[test]
    fn header_captures_obs_types_per_constellation() {
        let file = parse(SAMPLE).unwrap();
        let g = file.header.obs_types.get(&'G').unwrap();
        assert_eq!(g, &vec!["C1C".to_string(), "L1C".to_string()]);
        let j = file.header.obs_types.get(&'J').unwrap();
        assert_eq!(j, &vec!["C1C".to_string(), "L1C".to_string()]);
    }

    #[test]
    fn epoch_time_is_captured() {
        let file = parse(SAMPLE).unwrap();
        let ep = &file.epochs[0];
        assert_eq!((ep.year, ep.month, ep.day), (2025, 6, 24));
        assert_eq!((ep.hour, ep.minute), (12, 0));
    }

    #[test]
    fn two_satellites_recorded() {
        let file = parse(SAMPLE).unwrap();
        let ep = &file.epochs[0];
        assert_eq!(ep.observations.len(), 2);
        assert_eq!(ep.observations[0].prn, "G01");
        assert_eq!(ep.observations[1].prn, "J01");
    }

    #[test]
    fn observation_values_decode() {
        let file = parse(SAMPLE).unwrap();
        let obs = &file.epochs[0].observations[0];
        assert_eq!(obs.values.len(), 2);
        assert!((obs.values[0].unwrap() - 21_005_432.123).abs() < 1e-3);
        assert!((obs.values[1].unwrap() - 120_345_678.9).abs() < 1e-1);
    }

    #[test]
    fn missing_header_terminator_errors() {
        let bad = "     3.05           O                   M                   RINEX VERSION / TYPE\nno terminator";
        assert!(matches!(parse(bad), Err(RinexError::UnterminatedHeader)));
    }

    #[test]
    fn short_obs_line_reports_error() {
        // Craft a snippet with a bogus obs line.
        let bad = "     3.05           O                   M                   RINEX VERSION / TYPE
G    1 C1C                                                  SYS / # / OBS TYPES
                                                            END OF HEADER
> 2025 06 24 12 00  0.0000000  0  1
G";
        assert!(matches!(parse(bad), Err(RinexError::TooShort { .. })));
    }
}
