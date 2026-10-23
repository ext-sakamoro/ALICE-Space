//! CCSDS Orbit Mean-Elements Message (OMM, CCSDS 502.0-B-3) reader for SGP4
//! element sets.
//!
//! An OMM carries the same mean elements as a TLE, with a four-digit year, a
//! calendar epoch and catalogue numbers past 99999. This module reads the four
//! encodings in use for SGP4 element sets — KVN (`KEY = value` lines), the NDM
//! XML of CCSDS 505.0, and the JSON and CSV that CelesTrak serves for its GP data
//! (one object or row per element set, keys named as in the standard) — into
//! [`Omm`], and converts that to the [`MeanElements`] of [`crate::sgp4`] with the
//! same unit conversions as [`MeanElements::from_tle`], so an OMM and the TLE of
//! the same element set give the same SGP4 input (the epoch to the precision the
//! OMM prints it).
//!
//! Metadata that SGP4 depends on is checked when present: `CENTER_NAME` must be
//! `EARTH`, `REF_FRAME` `TEME`, `TIME_SYSTEM` `UTC` and `MEAN_ELEMENT_THEORY`
//! `SGP4` (CelesTrak's CSV has no metadata columns).
//!
//! Author: Moroya Sakamoto

use std::collections::BTreeMap;

use crate::sgp4::{days_from_1949_12_31, epoch_days_since_1950, MeanElements};

/// Why a message could not be read.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum OmmError {
    /// A required key is missing (the key).
    Missing(&'static str),
    /// A value is not of the form its key requires (key, value).
    Invalid(String, String),
    /// A metadata value names something SGP4 elements cannot be (key, value).
    Unsupported(String, String),
    /// The text is not in the encoding being read (what was expected).
    Syntax(String),
}

impl core::fmt::Display for OmmError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Missing(k) => write!(f, "OMM: required key {k} is missing"),
            Self::Invalid(k, v) => write!(f, "OMM: {k} = {v:?} is not a valid value"),
            Self::Unsupported(k, v) => {
                write!(f, "OMM: {k} = {v:?} is not supported for SGP4 elements")
            }
            Self::Syntax(s) => write!(f, "OMM: {s}"),
        }
    }
}

impl std::error::Error for OmmError {}

/// One element set of an OMM (the values as the message prints them).
#[derive(Clone, Debug, PartialEq)]
pub struct Omm {
    /// `OBJECT_NAME` (empty when absent).
    pub object_name: String,
    /// `OBJECT_ID`, the international designator `YYYY-NNNP{PP}` (empty when
    /// absent or `UNKNOWN`).
    pub object_id: String,
    /// `EPOCH` as days since 1949 December 31 00:00 UTC.
    pub epoch_days_since_1950: f64,
    /// `MEAN_MOTION` (rev/day).
    pub mean_motion_rev_per_day: f64,
    /// `ECCENTRICITY`.
    pub eccentricity: f64,
    /// `INCLINATION` (deg).
    pub inclination_deg: f64,
    /// `RA_OF_ASC_NODE` (deg).
    pub raan_deg: f64,
    /// `ARG_OF_PERICENTER` (deg).
    pub arg_pericenter_deg: f64,
    /// `MEAN_ANOMALY` (deg).
    pub mean_anomaly_deg: f64,
    /// `EPHEMERIS_TYPE` (0 when absent).
    pub ephemeris_type: u32,
    /// `CLASSIFICATION_TYPE` (`U` when absent).
    pub classification: char,
    /// `NORAD_CAT_ID`.
    pub norad_cat_id: u32,
    /// `ELEMENT_SET_NO` (0 when absent).
    pub element_set_no: u32,
    /// `REV_AT_EPOCH` (0 when absent).
    pub rev_at_epoch: u32,
    /// `BSTAR` (1 / Earth radii).
    pub bstar: f64,
    /// `MEAN_MOTION_DOT` (rev/day²), the first derivative divided by two as in
    /// a TLE.
    pub mean_motion_dot: f64,
    /// `MEAN_MOTION_DDOT` (rev/day³), the second derivative divided by six as in
    /// a TLE.
    pub mean_motion_ddot: f64,
}

impl Omm {
    /// The SGP4 mean elements of this element set.
    #[must_use]
    pub fn to_mean_elements(&self) -> MeanElements {
        // the same expressions as MeanElements::from_tle (via the TLE's radians)
        let xpdotp = 1440.0 / (2.0 * core::f64::consts::PI);
        MeanElements {
            epoch_days_since_1950: self.epoch_days_since_1950,
            bstar: self.bstar,
            ndot: self.mean_motion_dot / (xpdotp * 1440.0),
            nddot: self.mean_motion_ddot / (xpdotp * 1440.0 * 1440.0),
            eccentricity: self.eccentricity,
            arg_perigee: self.arg_pericenter_deg.to_radians(),
            inclination: self.inclination_deg.to_radians(),
            mean_anomaly: self.mean_anomaly_deg.to_radians(),
            mean_motion_kozai: self.mean_motion_rev_per_day / xpdotp,
            raan: self.raan_deg.to_radians(),
        }
    }

    /// Read one element set from its fields (key → value as text).
    ///
    /// # Errors
    ///
    /// [`OmmError`] for a missing required key, a malformed value, or metadata
    /// that SGP4 elements cannot have.
    pub fn from_fields(fields: &BTreeMap<String, String>) -> Result<Self, OmmError> {
        let get = |k: &'static str| fields.get(k).map(|v| v.trim()).filter(|v| !v.is_empty());
        let req = |k: &'static str| get(k).ok_or(OmmError::Missing(k));
        let num = |k: &'static str| -> Result<f64, OmmError> {
            let v = req(k)?;
            parse_f64(v).ok_or_else(|| OmmError::Invalid(k.to_owned(), v.to_owned()))
        };
        let int = |k: &'static str, default: u32| -> Result<u32, OmmError> {
            match get(k) {
                None => Ok(default),
                Some(v) => v
                    .parse()
                    .map_err(|_| OmmError::Invalid(k.to_owned(), v.to_owned())),
            }
        };

        for (key, want) in [
            ("CENTER_NAME", "EARTH"),
            ("REF_FRAME", "TEME"),
            ("TIME_SYSTEM", "UTC"),
            ("MEAN_ELEMENT_THEORY", "SGP4"),
        ] {
            if let Some(v) = get(key) {
                if !v.eq_ignore_ascii_case(want) {
                    return Err(OmmError::Unsupported(key.to_owned(), v.to_owned()));
                }
            }
        }

        let epoch_text = req("EPOCH")?;
        let epoch_days_since_1950 = parse_epoch(epoch_text)
            .ok_or_else(|| OmmError::Invalid("EPOCH".to_owned(), epoch_text.to_owned()))?;

        let eccentricity = num("ECCENTRICITY")?;
        if !(0.0..1.0).contains(&eccentricity) {
            return Err(OmmError::Invalid(
                "ECCENTRICITY".to_owned(),
                req("ECCENTRICITY")?.to_owned(),
            ));
        }
        let mean_motion_rev_per_day = num("MEAN_MOTION")?;
        if mean_motion_rev_per_day <= 0.0 {
            return Err(OmmError::Invalid(
                "MEAN_MOTION".to_owned(),
                req("MEAN_MOTION")?.to_owned(),
            ));
        }

        let classification = match get("CLASSIFICATION_TYPE") {
            None => 'U',
            Some(v) => {
                let mut c = v.chars();
                match (c.next(), c.next()) {
                    (Some(ch), None) => ch,
                    _ => {
                        return Err(OmmError::Invalid(
                            "CLASSIFICATION_TYPE".to_owned(),
                            v.to_owned(),
                        ))
                    }
                }
            }
        };
        let object_id = get("OBJECT_ID")
            .filter(|v| !v.eq_ignore_ascii_case("UNKNOWN"))
            .unwrap_or("")
            .to_owned();

        Ok(Self {
            object_name: get("OBJECT_NAME").unwrap_or("").to_owned(),
            object_id,
            epoch_days_since_1950,
            mean_motion_rev_per_day,
            eccentricity,
            inclination_deg: num("INCLINATION")?,
            raan_deg: num("RA_OF_ASC_NODE")?,
            arg_pericenter_deg: num("ARG_OF_PERICENTER")?,
            mean_anomaly_deg: num("MEAN_ANOMALY")?,
            ephemeris_type: int("EPHEMERIS_TYPE", 0)?,
            classification,
            norad_cat_id: {
                let v = req("NORAD_CAT_ID")?;
                v.parse()
                    .map_err(|_| OmmError::Invalid("NORAD_CAT_ID".to_owned(), v.to_owned()))?
            },
            element_set_no: int("ELEMENT_SET_NO", 0)?,
            rev_at_epoch: int("REV_AT_EPOCH", 0)?,
            bstar: num("BSTAR")?,
            mean_motion_dot: num("MEAN_MOTION_DOT")?,
            mean_motion_ddot: num("MEAN_MOTION_DDOT")?,
        })
    }
}

/// A decimal number as an OMM prints it (`.1568E-2`, `-0.00000016`, `15.99`):
/// Rust's correctly rounded parse after checking the form (no `inf` / `nan`).
fn parse_f64(v: &str) -> Option<f64> {
    let ok = !v.is_empty()
        && v.bytes()
            .all(|b| b.is_ascii_digit() || matches!(b, b'+' | b'-' | b'.' | b'e' | b'E'))
        && v.bytes().any(|b| b.is_ascii_digit());
    if !ok {
        return None;
    }
    v.parse().ok().filter(|x: &f64| x.is_finite())
}

fn is_leap(year: i64) -> bool {
    (year % 4 == 0 && year % 100 != 0) || year % 400 == 0
}

/// The OMM epoch (`YYYY-MM-DDThh:mm:ss[.f…]` or `YYYY-DDDThh:mm:ss[.f…]`, an
/// optional trailing `Z`) as days since 1949 December 31 00:00, computed as a
/// TLE epoch of the same year and day of the year
/// ([`crate::sgp4::epoch_days_since_1950`]). A leap second (`60`) is rejected:
/// SGP4 time has no place for it.
fn parse_epoch(text: &str) -> Option<f64> {
    let t = text.strip_suffix('Z').unwrap_or(text);
    let (date, time) = t.split_once('T')?;
    let digits = |s: &str, n: usize| -> Option<i64> {
        (s.len() == n && s.bytes().all(|b| b.is_ascii_digit())).then(|| s.parse().ok())?
    };
    let parts: Vec<&str> = date.split('-').collect();
    let day_number = match parts.as_slice() {
        [y, m, d] => {
            let (y, m, d) = (digits(y, 4)?, digits(m, 2)?, digits(d, 2)?);
            let month_len = [
                31,
                if is_leap(y) { 29 } else { 28 },
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
            if !(1..=12).contains(&m) || d < 1 || d > month_len[usize::try_from(m - 1).ok()?] {
                return None;
            }
            days_from_1949_12_31(y, m, d)
        }
        [y, doy] => {
            let (y, doy) = (digits(y, 4)?, digits(doy, 3)?);
            if doy < 1 || doy > if is_leap(y) { 366 } else { 365 } {
                return None;
            }
            days_from_1949_12_31(y, 1, 1) + doy - 1
        }
        _ => return None,
    };
    let tparts: Vec<&str> = time.split(':').collect();
    let [h, mi, s] = tparts.as_slice() else {
        return None;
    };
    let (h, mi) = (digits(h, 2)?, digits(mi, 2)?);
    let (s_int, s_frac) = match s.split_once('.') {
        Some((i, f)) => (digits(i, 2)?, Some(f)),
        None => (digits(s, 2)?, None),
    };
    if h > 23 || mi > 59 || s_int > 59 {
        return None;
    }
    let frac = match s_frac {
        None => 0.0,
        Some(f) if !f.is_empty() && f.bytes().all(|b| b.is_ascii_digit()) => {
            format!("0.{f}").parse::<f64>().ok()?
        }
        Some(_) => return None,
    };
    // the day of the year (1.0 = January 1 00:00), then the reference's epoch
    // arithmetic shared with the TLE path
    let year = day_number_year(day_number);
    let jan1 = days_from_1949_12_31(year, 1, 1);
    #[allow(clippy::cast_precision_loss)] // a day count and seconds far below 2^53
    let day_of_year =
        (day_number - jan1 + 1) as f64 + ((h * 3600 + mi * 60 + s_int) as f64 + frac) / 86_400.0;
    Some(epoch_days_since_1950(
        u32::try_from(year).ok()?,
        day_of_year,
    ))
}

/// The calendar year of a day counted from 1949 December 31.
fn day_number_year(day_number: i64) -> i64 {
    let mut year = 1950 + day_number.div_euclid(366);
    while days_from_1949_12_31(year + 1, 1, 1) <= day_number {
        year += 1;
    }
    while days_from_1949_12_31(year, 1, 1) > day_number {
        year -= 1;
    }
    year
}

/// Read every element set of a KVN OMM (`KEY = value` lines; `COMMENT` lines,
/// blank lines and the `CCSDS_OMM_VERS` header are skipped; a new `EPOCH` or
/// `OBJECT_NAME` after an `EPOCH` starts the next element set).
///
/// # Errors
///
/// [`OmmError`] for a line that is not `KEY = value`, an empty message, or an
/// element set that [`Omm::from_fields`] rejects.
pub fn parse_kvn(text: &str) -> Result<Vec<Omm>, OmmError> {
    let mut sets: Vec<BTreeMap<String, String>> = Vec::new();
    let mut cur: BTreeMap<String, String> = BTreeMap::new();
    for raw in text.lines() {
        let line = raw.trim();
        if line.is_empty() || line.starts_with("COMMENT") {
            continue;
        }
        let (key, value) = line
            .split_once('=')
            .ok_or_else(|| OmmError::Syntax(format!("line {line:?} is not KEY = value")))?;
        let key = key.trim().to_ascii_uppercase();
        // the unit annotation of KVN values, e.g. `15.5 [rev/day]`
        let value = value.split('[').next().unwrap_or("").trim().to_owned();
        if key == "CCSDS_OMM_VERS" {
            continue;
        }
        if (key == "OBJECT_NAME" || key == "EPOCH") && cur.contains_key("EPOCH") {
            sets.extend([std::mem::take(&mut cur)]);
        }
        cur.insert(key, value);
    }
    if !cur.is_empty() {
        sets.extend([cur]);
    }
    if sets.is_empty() {
        return Err(OmmError::Syntax("no element set".to_owned()));
    }
    sets.iter().map(Omm::from_fields).collect()
}

/// Read every element set of an NDM XML OMM (each `<segment>` with its
/// `<metadata>`, `<meanElements>` and `<tleParameters>` leaf elements).
///
/// # Errors
///
/// [`OmmError`] for a document without segments or an element set that
/// [`Omm::from_fields`] rejects.
pub fn parse_xml(text: &str) -> Result<Vec<Omm>, OmmError> {
    let mut out = Vec::new();
    let mut rest = text;
    while let Some(start) = rest.find("<segment") {
        let after = &rest[start..];
        let end = after
            .find("</segment>")
            .ok_or_else(|| OmmError::Syntax("<segment> without </segment>".to_owned()))?;
        let segment = &after[..end];
        out.extend([Omm::from_fields(&xml_leaves(segment)?)?]);
        rest = &after[end + "</segment>".len()..];
    }
    if out.is_empty() {
        return Err(OmmError::Syntax("no <segment>".to_owned()));
    }
    Ok(out)
}

/// The leaf elements (`<KEY>value</KEY>` and empty `<KEY/>`) of an XML
/// fragment, by tag name, entity references decoded.
fn xml_leaves(fragment: &str) -> Result<BTreeMap<String, String>, OmmError> {
    let mut map = BTreeMap::new();
    let mut rest = fragment;
    while let Some(lt) = rest.find('<') {
        rest = &rest[lt + 1..];
        let gt = rest
            .find('>')
            .ok_or_else(|| OmmError::Syntax("unterminated tag".to_owned()))?;
        let tag = &rest[..gt];
        rest = &rest[gt + 1..];
        if tag.starts_with('/') || tag.starts_with('?') || tag.starts_with('!') {
            continue;
        }
        let name = tag
            .split(|c: char| c.is_whitespace() || c == '/')
            .next()
            .unwrap_or("");
        if tag.ends_with('/') {
            map.insert(name.to_owned(), String::new());
            continue;
        }
        let close = format!("</{name}>");
        // a leaf: text up to its own closing tag, no child element in between
        if let Some(pos) = rest.find(&close) {
            let body = &rest[..pos];
            if !body.contains('<') {
                map.insert(name.to_owned(), xml_unescape(body.trim()));
                rest = &rest[pos + close.len()..];
            }
        }
    }
    Ok(map)
}

fn xml_unescape(s: &str) -> String {
    s.replace("&lt;", "<")
        .replace("&gt;", ">")
        .replace("&quot;", "\"")
        .replace("&apos;", "'")
        .replace("&amp;", "&")
}

/// Read every element set of CelesTrak's GP JSON (an array of flat objects, or
/// one object; values are strings, numbers or `null`).
///
/// # Errors
///
/// [`OmmError::Syntax`] for text that is not such JSON, and the errors of
/// [`Omm::from_fields`].
pub fn parse_json(text: &str) -> Result<Vec<Omm>, OmmError> {
    let mut p = JsonFlat {
        s: text.as_bytes(),
        i: 0,
    };
    p.ws();
    let objects = if p.peek() == Some(b'[') {
        p.i += 1;
        let mut v = Vec::new();
        p.ws();
        if p.peek() == Some(b']') {
            p.i += 1;
        } else {
            loop {
                v.extend([p.object()?]);
                p.ws();
                match p.next() {
                    Some(b',') => {}
                    Some(b']') => break,
                    _ => return Err(p.err("',' or ']'")),
                }
            }
        }
        v
    } else {
        vec![p.object()?]
    };
    p.ws();
    if p.i != p.s.len() {
        return Err(p.err("end of text"));
    }
    if objects.is_empty() {
        return Err(OmmError::Syntax("no element set".to_owned()));
    }
    objects.iter().map(Omm::from_fields).collect()
}

/// A reader for JSON arrays of flat objects (the shape of GP data).
struct JsonFlat<'a> {
    s: &'a [u8],
    i: usize,
}

impl JsonFlat<'_> {
    fn peek(&self) -> Option<u8> {
        self.s.get(self.i).copied()
    }
    fn next(&mut self) -> Option<u8> {
        let c = self.peek();
        self.i += 1;
        c
    }
    fn ws(&mut self) {
        while matches!(self.peek(), Some(b' ' | b'\t' | b'\n' | b'\r')) {
            self.i += 1;
        }
    }
    fn err(&self, want: &str) -> OmmError {
        OmmError::Syntax(format!("JSON: expected {want} at byte {}", self.i))
    }
    fn object(&mut self) -> Result<BTreeMap<String, String>, OmmError> {
        self.ws();
        if self.next() != Some(b'{') {
            return Err(self.err("'{'"));
        }
        let mut map = BTreeMap::new();
        self.ws();
        if self.peek() == Some(b'}') {
            self.i += 1;
            return Ok(map);
        }
        loop {
            self.ws();
            let key = self.string()?;
            self.ws();
            if self.next() != Some(b':') {
                return Err(self.err("':'"));
            }
            self.ws();
            let value = match self.peek() {
                Some(b'"') => self.string()?,
                Some(b'n') => {
                    self.literal("null")?;
                    String::new()
                }
                Some(b't') => {
                    self.literal("true")?;
                    "true".to_owned()
                }
                Some(b'f') => {
                    self.literal("false")?;
                    "false".to_owned()
                }
                _ => self.number()?,
            };
            map.insert(key, value);
            self.ws();
            match self.next() {
                Some(b',') => {}
                Some(b'}') => return Ok(map),
                _ => return Err(self.err("',' or '}'")),
            }
        }
    }
    fn literal(&mut self, word: &str) -> Result<(), OmmError> {
        if self.s[self.i..].starts_with(word.as_bytes()) {
            self.i += word.len();
            Ok(())
        } else {
            Err(self.err(word))
        }
    }
    fn number(&mut self) -> Result<String, OmmError> {
        let start = self.i;
        while matches!(
            self.peek(),
            Some(b'0'..=b'9' | b'+' | b'-' | b'.' | b'e' | b'E')
        ) {
            self.i += 1;
        }
        if self.i == start {
            return Err(self.err("a value"));
        }
        Ok(String::from_utf8_lossy(&self.s[start..self.i]).into_owned())
    }
    fn string(&mut self) -> Result<String, OmmError> {
        if self.next() != Some(b'"') {
            return Err(self.err("'\"'"));
        }
        let mut out: Vec<u8> = Vec::new();
        loop {
            match self.next() {
                None => return Err(self.err("closing '\"'")),
                Some(b'"') => break,
                Some(b'\\') => match self.next() {
                    Some(b'"') => out.extend([b'"']),
                    Some(b'\\') => out.extend([b'\\']),
                    Some(b'/') => out.extend([b'/']),
                    Some(b'n') => out.extend([b'\n']),
                    Some(b't') => out.extend([b'\t']),
                    Some(b'r') => out.extend([b'\r']),
                    Some(b'b') => out.extend([0x08]),
                    Some(b'f') => out.extend([0x0c]),
                    Some(b'u') => {
                        let hex = self
                            .s
                            .get(self.i..self.i + 4)
                            .and_then(|h| core::str::from_utf8(h).ok())
                            .and_then(|h| u32::from_str_radix(h, 16).ok())
                            .ok_or_else(|| self.err("4 hex digits"))?;
                        self.i += 4;
                        let ch = char::from_u32(hex).ok_or_else(|| self.err("a scalar value"))?;
                        let mut buf = [0u8; 4];
                        out.extend(ch.encode_utf8(&mut buf).as_bytes());
                    }
                    _ => return Err(self.err("an escape")),
                },
                Some(c) => out.extend([c]),
            }
        }
        String::from_utf8(out).map_err(|_| self.err("UTF-8"))
    }
}

/// Read every element set of CelesTrak's GP CSV (a header row of keys, then
/// one row per element set; fields contain no commas or quotes).
///
/// # Errors
///
/// [`OmmError::Syntax`] for an empty text or a row whose length differs from
/// the header's, and the errors of [`Omm::from_fields`].
pub fn parse_csv(text: &str) -> Result<Vec<Omm>, OmmError> {
    let mut lines = text.lines().map(str::trim).filter(|l| !l.is_empty());
    let header: Vec<&str> = lines
        .next()
        .ok_or_else(|| OmmError::Syntax("CSV: no header row".to_owned()))?
        .split(',')
        .map(str::trim)
        .collect();
    let mut out = Vec::new();
    for (n, line) in lines.enumerate() {
        let cells: Vec<&str> = line.split(',').map(str::trim).collect();
        if cells.len() != header.len() {
            return Err(OmmError::Syntax(format!(
                "CSV: row {} has {} cells, the header {}",
                n + 1,
                cells.len(),
                header.len()
            )));
        }
        let fields = header
            .iter()
            .zip(&cells)
            .map(|(k, v)| ((*k).to_ascii_uppercase(), (*v).to_owned()))
            .collect();
        out.extend([Omm::from_fields(&fields)?]);
    }
    if out.is_empty() {
        return Err(OmmError::Syntax("CSV: no element set".to_owned()));
    }
    Ok(out)
}
