//! QZSS L1S SLAS message framing.
//!
//! The L1S signal (Sub-meter Level Augmentation Service, SLAS) is
//! transmitted on 1575.42 MHz (same as GPS L1) at 250 sps and carries
//! SBAS-compatible augmentation messages plus Japan-specific safety
//! notifications (`DC Report Service`, disaster / crisis management).
//!
//! Each L1S message is 250 bits (32 bytes, last 6 bits padding) with the
//! following structure:
//!
//! ```text
//! bit  0-7   preamble  (0x53, 0x9A, 0xC6 cycle)
//! bit  8-13  message type (6 bits, 0..=63)
//! bit 14-...  data body (up to 212 bits)
//! last 24 bits Reed-Solomon parity
//! ```
//!
//! This module decodes the **preamble + message type + body/parity split**
//! required to feed downstream SLAS or DC Report parsers.
//!
//! # References
//!
//! - Cabinet Office Japan, "Quasi-Zenith Satellite System Interface
//!   Specification — Sub-Meter Level Augmentation Service" (`IS-QZSS-L1S-005`).
//! - Cabinet Office Japan, "Quasi-Zenith Satellite System Interface
//!   Specification — DC Report Service" (`IS-QZSS-DCR-013`).

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// L1S message length in bytes (250 bits + 6 pad = 32 bytes).
pub const MESSAGE_BYTES: usize = 32;
/// L1S carrier frequency in Hz.
pub const L1S_HZ: f64 = 1_575_420_000.0;
/// The three allowed preamble bytes, transmitted in cyclic sequence.
pub const PREAMBLES: [u8; 3] = [0x53, 0x9A, 0xC6];

// ---------------------------------------------------------------------------
// L1SMessageKind
// ---------------------------------------------------------------------------

/// Decoded L1S message type.
///
/// Type codes 0..=48 are SBAS-shared (per RTCA DO-229 Table A-6). Type
/// codes 51-53 are QZSS-specific SLAS augmentation. Type codes 54-56 are
/// reserved for the DC Report Service.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum L1SMessageKind {
    /// Type 1 — PRN mask assignment.
    PrnMask,
    /// Type 2/3/4/5 — fast corrections.
    FastCorrection,
    /// Type 6 — integrity information.
    Integrity,
    /// Type 7 — fast correction degradation.
    FastCorrectionDegradation,
    /// Type 9 — GEO ranging function parameters.
    GeoNavigation,
    /// Type 10 — degradation parameters.
    DegradationParameters,
    /// Type 17 — GEO satellite almanacs.
    GeoAlmanac,
    /// Type 18 — ionospheric grid point mask.
    IonoGridMask,
    /// Type 24 — mixed fast / long-term correction.
    MixedCorrection,
    /// Type 25 — long-term corrections.
    LongTermCorrection,
    /// Type 26 — ionospheric delay corrections.
    IonoDelayCorrection,
    /// Types 51-53 — QZSS-specific SLAS augmentation.
    Slas,
    /// Types 54-56 — QZSS DC Report (safety alerts, disaster info).
    DcReport,
    /// Type 63 — null message.
    Null,
    /// Reserved / unknown type code.
    Unknown,
}

impl L1SMessageKind {
    /// Short code used in canonical serialization.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::PrnMask => "MASK",
            Self::FastCorrection => "FAST",
            Self::Integrity => "INTEG",
            Self::FastCorrectionDegradation => "FASTDEG",
            Self::GeoNavigation => "GEONAV",
            Self::DegradationParameters => "DEGPARAM",
            Self::GeoAlmanac => "GEOALM",
            Self::IonoGridMask => "IGMASK",
            Self::MixedCorrection => "MIXED",
            Self::LongTermCorrection => "LTC",
            Self::IonoDelayCorrection => "IDC",
            Self::Slas => "SLAS",
            Self::DcReport => "DCR",
            Self::Null => "NULL",
            Self::Unknown => "UNK",
        }
    }

    /// Decode a 6-bit L1S / SBAS message type code.
    #[must_use]
    pub const fn from_type_code(code: u8) -> Self {
        match code & 0x3F {
            1 => Self::PrnMask,
            2..=5 => Self::FastCorrection,
            6 => Self::Integrity,
            7 => Self::FastCorrectionDegradation,
            9 => Self::GeoNavigation,
            10 => Self::DegradationParameters,
            17 => Self::GeoAlmanac,
            18 => Self::IonoGridMask,
            24 => Self::MixedCorrection,
            25 => Self::LongTermCorrection,
            26 => Self::IonoDelayCorrection,
            51..=53 => Self::Slas,
            54..=56 => Self::DcReport,
            63 => Self::Null,
            _ => Self::Unknown,
        }
    }
}

// ---------------------------------------------------------------------------
// L1SMessage
// ---------------------------------------------------------------------------

/// A parsed L1S message.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct L1SMessage {
    /// The preamble byte transmitted for this message.
    pub preamble: u8,
    /// Raw 6-bit message type code.
    pub type_code: u8,
    /// Decoded message kind.
    pub kind: L1SMessageKind,
    /// Data body bits (up to 212 bits, packed MSB-first into bytes).
    pub body: Vec<u8>,
}

// ---------------------------------------------------------------------------
// L1SParseError
// ---------------------------------------------------------------------------

/// Errors returned by [`parse_message`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum L1SParseError {
    /// The buffer was not exactly [`MESSAGE_BYTES`] bytes long.
    WrongMessageLength { got: usize, expected: usize },
    /// The preamble did not match any of [`PREAMBLES`].
    UnknownPreamble { got: u8 },
}

impl core::fmt::Display for L1SParseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::WrongMessageLength { got, expected } => {
                write!(f, "L1S message length {got} != expected {expected}")
            }
            Self::UnknownPreamble { got } => write!(f, "L1S unknown preamble 0x{got:02X}"),
        }
    }
}

impl std::error::Error for L1SParseError {}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a raw L1S message frame.
///
/// The buffer must be exactly [`MESSAGE_BYTES`] bytes; the first byte is
/// the preamble, the next 6 bits are the type code (in the top of the
/// second byte), and the remainder is the body + RS parity.
///
/// # Errors
///
/// Returns [`L1SParseError`] variants if the length or preamble is invalid.
pub fn parse_message(buf: &[u8]) -> Result<L1SMessage, L1SParseError> {
    if buf.len() != MESSAGE_BYTES {
        return Err(L1SParseError::WrongMessageLength {
            got: buf.len(),
            expected: MESSAGE_BYTES,
        });
    }
    let preamble = buf[0];
    if !PREAMBLES.contains(&preamble) {
        return Err(L1SParseError::UnknownPreamble { got: preamble });
    }
    let type_code = (buf[1] >> 2) & 0x3F;
    let kind = L1SMessageKind::from_type_code(type_code);
    let body = buf[1..].to_vec();
    Ok(L1SMessage {
        preamble,
        type_code,
        kind,
        body,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_msg(preamble: u8, type_code: u8) -> Vec<u8> {
        let mut buf = vec![0u8; MESSAGE_BYTES];
        buf[0] = preamble;
        buf[1] = (type_code & 0x3F) << 2;
        buf
    }

    #[test]
    fn slas_type_51_is_recognized() {
        let m = parse_message(&make_msg(PREAMBLES[0], 51)).unwrap();
        assert_eq!(m.kind, L1SMessageKind::Slas);
    }

    #[test]
    fn dc_report_type_54_is_recognized() {
        let m = parse_message(&make_msg(PREAMBLES[0], 54)).unwrap();
        assert_eq!(m.kind, L1SMessageKind::DcReport);
    }

    #[test]
    fn fast_correction_types_2_to_5_are_all_fast() {
        for t in 2..=5u8 {
            let m = parse_message(&make_msg(PREAMBLES[0], t)).unwrap();
            assert_eq!(m.kind, L1SMessageKind::FastCorrection);
        }
    }

    #[test]
    fn iono_delay_correction_type_26_is_recognized() {
        let m = parse_message(&make_msg(PREAMBLES[0], 26)).unwrap();
        assert_eq!(m.kind, L1SMessageKind::IonoDelayCorrection);
    }

    #[test]
    fn null_type_63_is_recognized() {
        let m = parse_message(&make_msg(PREAMBLES[0], 63)).unwrap();
        assert_eq!(m.kind, L1SMessageKind::Null);
    }

    #[test]
    fn unknown_type_code_returns_unknown() {
        let m = parse_message(&make_msg(PREAMBLES[0], 40)).unwrap();
        assert_eq!(m.kind, L1SMessageKind::Unknown);
    }

    #[test]
    fn all_three_preambles_are_accepted() {
        for &p in &PREAMBLES {
            let m = parse_message(&make_msg(p, 1)).unwrap();
            assert_eq!(m.preamble, p);
        }
    }

    #[test]
    fn unknown_preamble_is_rejected() {
        let err = parse_message(&make_msg(0x00, 1)).unwrap_err();
        assert_eq!(err, L1SParseError::UnknownPreamble { got: 0x00 });
    }

    #[test]
    fn wrong_length_is_rejected() {
        let err = parse_message(&[0x53u8; 5]).unwrap_err();
        assert_eq!(
            err,
            L1SParseError::WrongMessageLength {
                got: 5,
                expected: MESSAGE_BYTES,
            }
        );
    }

    #[test]
    fn code_labels_are_stable() {
        assert_eq!(L1SMessageKind::Slas.code(), "SLAS");
        assert_eq!(L1SMessageKind::DcReport.code(), "DCR");
        assert_eq!(L1SMessageKind::FastCorrection.code(), "FAST");
        assert_eq!(L1SMessageKind::PrnMask.code(), "MASK");
        assert_eq!(L1SMessageKind::Null.code(), "NULL");
    }

    #[test]
    fn error_display_shows_length() {
        let err = L1SParseError::WrongMessageLength {
            got: 12,
            expected: 32,
        };
        let s = err.to_string();
        assert!(s.contains("12"));
    }

    #[test]
    fn error_display_shows_preamble() {
        let err = L1SParseError::UnknownPreamble { got: 0xFF };
        let s = err.to_string();
        assert!(s.contains("FF"));
    }
}
