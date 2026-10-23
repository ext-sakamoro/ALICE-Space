//! CCSDS Orbit Data Messages (CCSDS 502.0-B-3) in KVN: the Orbit Ephemeris
//! Message (OEM, a table of states per segment) and the Orbit Parameter Message
//! (OPM, one state with optional osculating elements).
//!
//! An OEM segment is read with its metadata, its useable time window and the
//! interpolation it names (`INTERPOLATION` = `LAGRANGE` / `HERMITE` / `LINEAR`,
//! `INTERPOLATION_DEGREE`); [`OemSegment::state_at`] interpolates the table with
//! [`crate::interp`]. Covariance blocks are skipped. Epochs are held as whole
//! days and seconds of the day ([`OdmEpoch`]), so times relative to a segment
//! keep full precision whatever the epoch. Second 60 (`23:59:60`) is rejected
//! wherever an epoch occurs (rows, window times, the OPM epoch, the epoch asked
//! of [`OemSegment::state_at`]): in UTC as unsupported, since uniform seconds
//! across a leap second need the time-scale tables, and in any other time
//! system as invalid, since it has no such second.
//!
//! Author: Moroya Sakamoto

use std::collections::BTreeMap;

use crate::interp::{self, InterpError, State};
use crate::sgp4::days_from_1949_12_31;

/// Why a message could not be read or interpolated.
#[derive(Clone, Debug, PartialEq)]
pub enum OdmError {
    /// A required key is missing (the key).
    Missing(&'static str),
    /// A value is not of the form its key requires (key, value).
    Invalid(String, String),
    /// A value the reader does not support (key, value).
    Unsupported(String, String),
    /// The text does not follow the KVN structure (what was found, line number).
    Syntax(String, usize),
    /// The requested epoch is outside the segment's useable window.
    OutsideUseableWindow,
    /// The segment names no interpolation and none was given.
    NoInterpolation,
    /// The interpolation failed.
    Interp(InterpError),
}

impl core::fmt::Display for OdmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Missing(k) => write!(f, "ODM: required key {k} is missing"),
            Self::Invalid(k, v) => write!(f, "ODM: {k} = {v:?} is not a valid value"),
            Self::Unsupported(k, v) => write!(f, "ODM: {k} = {v:?} is not supported"),
            Self::Syntax(s, n) => write!(f, "ODM: line {n}: {s}"),
            Self::OutsideUseableWindow => write!(f, "ODM: epoch outside the useable window"),
            Self::NoInterpolation => write!(f, "ODM: the segment names no interpolation"),
            Self::Interp(e) => write!(f, "ODM: {e}"),
        }
    }
}

impl std::error::Error for OdmError {}

impl From<InterpError> for OdmError {
    fn from(e: InterpError) -> Self {
        Self::Interp(e)
    }
}

/// An epoch as whole days since 1949 December 31 and seconds of that day.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OdmEpoch {
    /// Whole days since 1949 December 31.
    pub day: i64,
    /// Seconds of the day (0 ≤ s < 86400, or < 86401 on a leap second).
    pub second_of_day: f64,
}

impl OdmEpoch {
    /// Parse `YYYY-MM-DDThh:mm:ss[.f…]` or `YYYY-DDDThh:mm:ss[.f…]` (an optional
    /// trailing `Z`); `None` when it is not such a time. Second 60 is accepted
    /// (the caller decides whether the time system has leap seconds).
    #[must_use]
    pub fn parse(text: &str) -> Option<Self> {
        let t = text.trim();
        let t = t.strip_suffix('Z').unwrap_or(t);
        let (date, time) = t.split_once('T')?;
        let digits = |s: &str, n: usize| -> Option<i64> {
            (s.len() == n && s.bytes().all(|b| b.is_ascii_digit())).then(|| s.parse().ok())?
        };
        let leap = |y: i64| (y % 4 == 0 && y % 100 != 0) || y % 400 == 0;
        let parts: Vec<&str> = date.split('-').collect();
        let day = match parts.as_slice() {
            [y, m, d] => {
                let (y, m, d) = (digits(y, 4)?, digits(m, 2)?, digits(d, 2)?);
                let len = [
                    31,
                    if leap(y) { 29 } else { 28 },
                    31,
                    30,
                    31,
                    30,
                    31,
                    31,
                    30,
                    31,
                    30,
                    31,
                ];
                if !(1..=12).contains(&m) || d < 1 || d > len[usize::try_from(m - 1).ok()?] {
                    return None;
                }
                days_from_1949_12_31(y, m, d)
            }
            [y, doy] => {
                let (y, doy) = (digits(y, 4)?, digits(doy, 3)?);
                if doy < 1 || doy > if leap(y) { 366 } else { 365 } {
                    return None;
                }
                days_from_1949_12_31(y, 1, 1) + doy - 1
            }
            _ => return None,
        };
        let tp: Vec<&str> = time.split(':').collect();
        let [h, mi, s] = tp.as_slice() else {
            return None;
        };
        let (h, mi) = (digits(h, 2)?, digits(mi, 2)?);
        let (si, frac) = match s.split_once('.') {
            Some((i, f)) if !f.is_empty() && f.bytes().all(|b| b.is_ascii_digit()) => {
                (digits(i, 2)?, format!("0.{f}").parse::<f64>().ok()?)
            }
            Some(_) => return None,
            None => (digits(s, 2)?, 0.0),
        };
        if h > 23 || mi > 59 || si > 60 || (si == 60 && (h != 23 || mi != 59)) {
            return None;
        }
        #[allow(clippy::cast_precision_loss)] // below 86401
        let second_of_day = (h * 3600 + mi * 60 + si) as f64 + frac;
        Some(Self { day, second_of_day })
    }

    /// Seconds from `origin` to `self`, counting every day as 86400 s.
    #[must_use]
    pub fn seconds_since(&self, origin: &Self) -> f64 {
        #[allow(clippy::cast_precision_loss)] // a day difference far below 2^53 / 86400
        let days = (self.day - origin.day) as f64;
        days * 86_400.0 + (self.second_of_day - origin.second_of_day)
    }

    fn is_leap_second(&self) -> bool {
        self.second_of_day >= 86_400.0
    }
}

/// How an OEM segment is to be interpolated.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Interpolation {
    /// Lagrange of the positions (the velocity is the derivative of the
    /// interpolating polynomial).
    Lagrange,
    /// Hermite of positions and velocities.
    Hermite,
    /// Linear (Lagrange of degree 1).
    Linear,
}

/// One OEM segment: metadata and the table of states.
#[derive(Clone, Debug, PartialEq)]
pub struct OemSegment {
    /// `OBJECT_NAME`.
    pub object_name: String,
    /// `OBJECT_ID`.
    pub object_id: String,
    /// `CENTER_NAME`.
    pub center_name: String,
    /// `REF_FRAME`.
    pub ref_frame: String,
    /// `TIME_SYSTEM`.
    pub time_system: String,
    /// `START_TIME`.
    pub start_time: OdmEpoch,
    /// `STOP_TIME`.
    pub stop_time: OdmEpoch,
    /// `USEABLE_START_TIME` (the start time when absent).
    pub useable_start_time: OdmEpoch,
    /// `USEABLE_STOP_TIME` (the stop time when absent).
    pub useable_stop_time: OdmEpoch,
    /// `INTERPOLATION` and `INTERPOLATION_DEGREE`, when given.
    pub interpolation: Option<(Interpolation, usize)>,
    /// Epochs of the rows.
    pub epochs: Vec<OdmEpoch>,
    /// Positions (km).
    pub positions_km: Vec<[f64; 3]>,
    /// Velocities (km/s).
    pub velocities_km_s: Vec<[f64; 3]>,
}

impl OemSegment {
    /// The interpolated state (km, km/s) at `epoch` with the segment's own
    /// interpolation.
    ///
    /// # Errors
    ///
    /// [`OdmError::NoInterpolation`] when the segment names none, and the errors
    /// of [`Self::state_at_with`].
    pub fn state_at(&self, epoch: &OdmEpoch) -> Result<State, OdmError> {
        let (method, degree) = self.interpolation.ok_or(OdmError::NoInterpolation)?;
        self.state_at_with(epoch, method, degree)
    }

    /// The interpolated state at `epoch` with `method` of `degree`.
    ///
    /// # Errors
    ///
    /// [`OdmError::OutsideUseableWindow`] outside the useable window,
    /// [`OdmError::Interp`] when the table is too short for the degree, and
    /// for second 60 [`OdmError::Unsupported`] (UTC) or [`OdmError::Invalid`]
    /// (other time systems).
    pub fn state_at_with(
        &self,
        epoch: &OdmEpoch,
        method: Interpolation,
        degree: usize,
    ) -> Result<State, OdmError> {
        reject_second_60(epoch, &self.time_system, "query epoch", "second 60")?;
        let origin = self.epochs[0];
        let t = epoch.seconds_since(&origin);
        let lo = self.useable_start_time.seconds_since(&origin);
        let hi = self.useable_stop_time.seconds_since(&origin);
        if !(lo..=hi).contains(&t) {
            return Err(OdmError::OutsideUseableWindow);
        }
        let times: Vec<f64> = self
            .epochs
            .iter()
            .map(|e| e.seconds_since(&origin))
            .collect();
        Ok(match method {
            Interpolation::Lagrange => interp::lagrange(&times, &self.positions_km, degree, t)?,
            Interpolation::Linear => interp::lagrange(&times, &self.positions_km, 1, t)?,
            Interpolation::Hermite => {
                interp::hermite(&times, &self.positions_km, &self.velocities_km_s, degree, t)?
            }
        })
    }
}

/// An OEM: its version and segments.
#[derive(Clone, Debug, PartialEq)]
pub struct Oem {
    /// `CCSDS_OEM_VERS`.
    pub version: String,
    /// The segments, in order.
    pub segments: Vec<OemSegment>,
}

/// A KVN line with its line number: comments and blank lines removed, unit
/// annotations (`[km]`) stripped from values.
enum Line<'a> {
    Kv(String, String),
    Raw(&'a str),
}

fn kvn_lines(text: &str) -> Vec<(usize, Line<'_>)> {
    let mut out = Vec::new();
    for (n, raw) in text.lines().enumerate() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with("COMMENT") {
            continue;
        }
        let entry = match line.split_once('=') {
            Some((k, v)) => Line::Kv(
                k.trim().to_ascii_uppercase(),
                v.split('[').next().unwrap_or("").trim().to_owned(),
            ),
            None => Line::Raw(line),
        };
        out.extend([(n + 1, entry)]);
    }
    out
}

fn epoch_of(map: &BTreeMap<String, String>, key: &'static str) -> Result<OdmEpoch, OdmError> {
    let v = map.get(key).ok_or(OdmError::Missing(key))?;
    OdmEpoch::parse(v).ok_or_else(|| OdmError::Invalid(key.to_owned(), v.clone()))
}

/// Second 60 in `time_system`: a UTC leap second is unsupported (uniform
/// seconds across it need the time-scale tables), and in any other time
/// system the second does not exist.
fn reject_second_60(
    e: &OdmEpoch,
    time_system: &str,
    key: &str,
    text: &str,
) -> Result<(), OdmError> {
    if !e.is_leap_second() {
        Ok(())
    } else if time_system.eq_ignore_ascii_case("UTC") {
        Err(OdmError::Unsupported(key.to_owned(), text.to_owned()))
    } else {
        Err(OdmError::Invalid(key.to_owned(), text.to_owned()))
    }
}

fn text_of(map: &BTreeMap<String, String>, key: &'static str) -> Result<String, OdmError> {
    map.get(key).cloned().ok_or(OdmError::Missing(key))
}

fn number(key: &str, v: &str) -> Result<f64, OdmError> {
    let ok = !v.is_empty()
        && v.bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'+' | b'-' | b'.' | b'e' | b'E'))
        && v.bytes().any(|b| b.is_ascii_digit());
    ok.then(|| v.parse::<f64>().ok())
        .flatten()
        .filter(|x| x.is_finite())
        .ok_or_else(|| OdmError::Invalid(key.to_owned(), v.to_owned()))
}

fn segment_from(
    meta: &BTreeMap<String, String>,
    rows: Vec<(usize, &str)>,
) -> Result<OemSegment, OdmError> {
    let time_system = text_of(meta, "TIME_SYSTEM")?;
    let interpolation = match meta.get("INTERPOLATION") {
        None => None,
        Some(m) => {
            let method = match m.to_ascii_uppercase().as_str() {
                "LAGRANGE" => Interpolation::Lagrange,
                "HERMITE" => Interpolation::Hermite,
                "LINEAR" => Interpolation::Linear,
                _ => return Err(OdmError::Unsupported("INTERPOLATION".to_owned(), m.clone())),
            };
            let degree = match (method, meta.get("INTERPOLATION_DEGREE")) {
                (Interpolation::Linear, None) => 1,
                (_, None) => return Err(OdmError::Missing("INTERPOLATION_DEGREE")),
                (_, Some(d)) => d
                    .parse()
                    .map_err(|_| OdmError::Invalid("INTERPOLATION_DEGREE".to_owned(), d.clone()))?,
            };
            Some((method, degree))
        }
    };
    let start_time = epoch_of(meta, "START_TIME")?;
    let stop_time = epoch_of(meta, "STOP_TIME")?;
    let useable_start_time = if meta.contains_key("USEABLE_START_TIME") {
        epoch_of(meta, "USEABLE_START_TIME")?
    } else {
        start_time
    };
    let useable_stop_time = if meta.contains_key("USEABLE_STOP_TIME") {
        epoch_of(meta, "USEABLE_STOP_TIME")?
    } else {
        stop_time
    };

    for (key, e) in [
        ("START_TIME", &start_time),
        ("STOP_TIME", &stop_time),
        ("USEABLE_START_TIME", &useable_start_time),
        ("USEABLE_STOP_TIME", &useable_stop_time),
    ] {
        if let Some(v) = meta.get(key) {
            reject_second_60(e, &time_system, key, v)?;
        }
    }

    let mut epochs = Vec::new();
    let mut positions_km = Vec::new();
    let mut velocities_km_s = Vec::new();
    for (n, row) in rows {
        let cols: Vec<&str> = row.split_whitespace().collect();
        if cols.len() != 7 && cols.len() != 10 {
            return Err(OdmError::Syntax(
                format!("an ephemeris row has 7 or 10 columns, found {}", cols.len()),
                n,
            ));
        }
        let e = OdmEpoch::parse(cols[0])
            .ok_or_else(|| OdmError::Invalid("EPOCH".to_owned(), cols[0].to_owned()))?;
        reject_second_60(&e, &time_system, "EPOCH", cols[0])?;
        if let Some(prev) = epochs.last() {
            if e.seconds_since(prev) <= 0.0 {
                return Err(OdmError::Syntax("epochs are not increasing".to_owned(), n));
            }
        }
        let v: Vec<f64> = cols[1..7]
            .iter()
            .map(|c| number("ephemeris value", c))
            .collect::<Result<_, _>>()?;
        epochs.extend([e]);
        positions_km.extend([[v[0], v[1], v[2]]]);
        velocities_km_s.extend([[v[3], v[4], v[5]]]);
    }
    if epochs.is_empty() {
        return Err(OdmError::Missing("ephemeris data"));
    }
    let first = epochs[0];
    let last = epochs[epochs.len() - 1];
    let within =
        |e: &OdmEpoch| e.seconds_since(&start_time) >= 0.0 && stop_time.seconds_since(e) >= 0.0;
    if !within(&first)
        || !within(&last)
        || !within(&useable_start_time)
        || !within(&useable_stop_time)
    {
        return Err(OdmError::Invalid(
            "START_TIME / STOP_TIME".to_owned(),
            "the data or the useable window is outside the segment".to_owned(),
        ));
    }
    Ok(OemSegment {
        object_name: text_of(meta, "OBJECT_NAME")?,
        object_id: text_of(meta, "OBJECT_ID")?,
        center_name: text_of(meta, "CENTER_NAME")?,
        ref_frame: text_of(meta, "REF_FRAME")?,
        time_system,
        start_time,
        stop_time,
        useable_start_time,
        useable_stop_time,
        interpolation,
        epochs,
        positions_km,
        velocities_km_s,
    })
}

/// Read a KVN OEM.
///
/// # Errors
///
/// [`OdmError`] for a missing version or metadata key, a malformed row or value,
/// rows outside the segment, an unsupported interpolation, or second 60
/// (unsupported in UTC, invalid in other time systems).
pub fn parse_oem_kvn(text: &str) -> Result<Oem, OdmError> {
    let lines = kvn_lines(text);
    let mut version = None;
    let mut segments = Vec::new();
    let mut i = 0;
    while i < lines.len() {
        let (n, line) = &lines[i];
        match line {
            Line::Kv(k, v) if k == "CCSDS_OEM_VERS" => version = Some(v.clone()),
            Line::Kv(k, _) if k == "CREATION_DATE" || k == "ORIGINATOR" || k == "MESSAGE_ID" => {}
            Line::Raw("META_START") => {
                let mut meta = BTreeMap::new();
                i += 1;
                loop {
                    match lines.get(i) {
                        Some((_, Line::Raw("META_STOP"))) => break,
                        Some((_, Line::Kv(k, v))) => {
                            meta.insert(k.clone(), v.clone());
                        }
                        Some((m, Line::Raw(r))) => {
                            return Err(OdmError::Syntax(format!("{r:?} inside metadata"), *m))
                        }
                        None => {
                            return Err(OdmError::Syntax(
                                "META_START without META_STOP".to_owned(),
                                *n,
                            ))
                        }
                    }
                    i += 1;
                }
                i += 1;
                let mut rows = Vec::new();
                while let Some((m, l)) = lines.get(i) {
                    match l {
                        Line::Raw("META_START") => break,
                        Line::Raw("COVARIANCE_START") => {
                            // skip the covariance block
                            while !matches!(
                                lines.get(i),
                                Some((_, Line::Raw("COVARIANCE_STOP"))) | None
                            ) {
                                i += 1;
                            }
                        }
                        Line::Raw(r) => rows.extend([(*m, *r)]),
                        Line::Kv(k, _) => {
                            return Err(OdmError::Syntax(format!("key {k} among the data"), *m))
                        }
                    }
                    i += 1;
                }
                segments.extend([segment_from(&meta, rows)?]);
                continue;
            }
            Line::Kv(k, _) => return Err(OdmError::Syntax(format!("unexpected key {k}"), *n)),
            Line::Raw(r) => return Err(OdmError::Syntax(format!("unexpected line {r:?}"), *n)),
        }
        i += 1;
    }
    let version = version.ok_or(OdmError::Missing("CCSDS_OEM_VERS"))?;
    if segments.is_empty() {
        return Err(OdmError::Missing("META_START"));
    }
    Ok(Oem { version, segments })
}

/// Osculating Keplerian elements of an OPM.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OpmKeplerian {
    /// `SEMI_MAJOR_AXIS` (km).
    pub semi_major_axis_km: f64,
    /// `ECCENTRICITY`.
    pub eccentricity: f64,
    /// `INCLINATION` (deg).
    pub inclination_deg: f64,
    /// `RA_OF_ASC_NODE` (deg).
    pub raan_deg: f64,
    /// `ARG_OF_PERICENTER` (deg).
    pub arg_pericenter_deg: f64,
    /// `TRUE_ANOMALY` or `MEAN_ANOMALY` (deg), and whether it is the true one.
    pub anomaly_deg: (f64, bool),
    /// `GM` (km³/s²).
    pub gm_km3_s2: f64,
}

/// An OPM: metadata, the state vector, and optional Keplerian elements.
#[derive(Clone, Debug, PartialEq)]
pub struct Opm {
    /// `CCSDS_OPM_VERS`.
    pub version: String,
    /// `OBJECT_NAME`.
    pub object_name: String,
    /// `OBJECT_ID`.
    pub object_id: String,
    /// `CENTER_NAME`.
    pub center_name: String,
    /// `REF_FRAME`.
    pub ref_frame: String,
    /// `TIME_SYSTEM`.
    pub time_system: String,
    /// `EPOCH`.
    pub epoch: OdmEpoch,
    /// `X`, `Y`, `Z` (km).
    pub position_km: [f64; 3],
    /// `X_DOT`, `Y_DOT`, `Z_DOT` (km/s).
    pub velocity_km_s: [f64; 3],
    /// The Keplerian elements block, when present.
    pub keplerian: Option<OpmKeplerian>,
}

/// Read a KVN OPM (the state vector and the optional Keplerian elements; the
/// spacecraft parameters, covariance and maneuvers are read past).
///
/// # Errors
///
/// [`OdmError`] for a missing or malformed required key or a line that is not
/// `KEY = value`.
pub fn parse_opm_kvn(text: &str) -> Result<Opm, OdmError> {
    let mut map = BTreeMap::new();
    for (n, line) in kvn_lines(text) {
        match line {
            Line::Kv(k, v) => {
                map.insert(k, v);
            }
            Line::Raw(r) => return Err(OdmError::Syntax(format!("{r:?} is not KEY = value"), n)),
        }
    }
    let num = |k: &'static str| -> Result<f64, OdmError> {
        number(k, map.get(k).ok_or(OdmError::Missing(k))?)
    };
    let keplerian = if map.contains_key("SEMI_MAJOR_AXIS") {
        let anomaly_deg = match (map.get("TRUE_ANOMALY"), map.get("MEAN_ANOMALY")) {
            (Some(v), _) => (number("TRUE_ANOMALY", v)?, true),
            (None, Some(v)) => (number("MEAN_ANOMALY", v)?, false),
            (None, None) => return Err(OdmError::Missing("TRUE_ANOMALY")),
        };
        Some(OpmKeplerian {
            semi_major_axis_km: num("SEMI_MAJOR_AXIS")?,
            eccentricity: num("ECCENTRICITY")?,
            inclination_deg: num("INCLINATION")?,
            raan_deg: num("RA_OF_ASC_NODE")?,
            arg_pericenter_deg: num("ARG_OF_PERICENTER")?,
            anomaly_deg,
            gm_km3_s2: num("GM")?,
        })
    } else {
        None
    };
    let epoch = epoch_of(&map, "EPOCH")?;
    let time_system = text_of(&map, "TIME_SYSTEM")?;
    reject_second_60(&epoch, &time_system, "EPOCH", &text_of(&map, "EPOCH")?)?;
    Ok(Opm {
        version: text_of(&map, "CCSDS_OPM_VERS")?,
        object_name: text_of(&map, "OBJECT_NAME")?,
        object_id: text_of(&map, "OBJECT_ID")?,
        center_name: text_of(&map, "CENTER_NAME")?,
        ref_frame: text_of(&map, "REF_FRAME")?,
        time_system,
        epoch,
        position_km: [num("X")?, num("Y")?, num("Z")?],
        velocity_km_s: [num("X_DOT")?, num("Y_DOT")?, num("Z_DOT")?],
        keplerian,
    })
}
