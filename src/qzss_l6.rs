//! QZSS L6 signal message framing.
//!
//! QZSS transmits four families of augmentation messages on the L6
//! frequency (1278.75 MHz):
//!
//! - **L6D MADOCA-PPP** — State Space Representation (SSR) corrections
//!   for global PPP.
//! - **L6E MADOCALIB** — an extension of MADOCA with additional
//!   augmentation fields (test-signal since 2022).
//! - **L6D CLAS** — Centimeter Level Augmentation Service, network-RTK
//!   style corrections for the Japanese archipelago.
//! - **L6E CLAS Test** — pre-operational CLAS variant.
//!
//! Every L6 frame is 2000 symbols (BPSK-modulated, 250 bytes) at
//! 250 sps. The first 4 bytes are the **PRN, MessageType, Alert Flag,
//! Vendor ID** header; the remainder is the message body plus a
//! Reed-Solomon (255,223) code.
//!
//! This module implements the **frame header + type dispatch** — the
//! entry point required by any SPACID-facing SDK. Full SSR / CLAS
//! decoders can be layered on top.
//!
//! # References
//!
//! - Cabinet Office Japan, "Quasi-Zenith Satellite System Interface
//!   Specification — Centimeter Level Augmentation Service" (`IS-QZSS-L6-005`).
//! - Cabinet Office Japan, "Quasi-Zenith Satellite System Interface
//!   Specification — MADOCA-PPP" (`IS-QZSS-MDC-002`).
//! - Cabinet Office Japan, "IS-QZSS-PNT-005" §4.1.2.5 (L6 signal).

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// L6 frame length in bytes (250 bytes = 2000 symbols BPSK).
pub const FRAME_BYTES: usize = 250;
/// L6 header length in bytes.
pub const HEADER_BYTES: usize = 4;
/// L6 carrier frequency in Hz.
pub const L6_HZ: f64 = 1_278_750_000.0;

// ---------------------------------------------------------------------------
// L6MessageKind
// ---------------------------------------------------------------------------

/// The kind of augmentation carried by an L6 frame.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum L6MessageKind {
    /// MADOCA-PPP State Space Representation corrections (`L6D`).
    MadocaPpp,
    /// MADOCALIB extension (`L6E` production).
    Madocalib,
    /// Centimeter Level Augmentation Service (`L6D CLAS`).
    Clas,
    /// CLAS test-signal (`L6E CLAS`).
    ClasTest,
    /// Unknown / reserved message type.
    Unknown,
}

impl L6MessageKind {
    /// Short code used in canonical serialization.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::MadocaPpp => "MADOCA",
            Self::Madocalib => "MADOCALIB",
            Self::Clas => "CLAS",
            Self::ClasTest => "CLAS-TEST",
            Self::Unknown => "UNK",
        }
    }

    /// Decode the 5-bit MessageType field of the L6 header.
    ///
    /// Per `IS-QZSS-L6-005` Table 4.1.2-1, the current type-code
    /// allocation is:
    ///
    /// | Code | Meaning |
    /// |------|---------|
    /// | 0x1F | Reserved / null |
    /// | 0x1D | CLAS (operational L6D) |
    /// | 0x1E | CLAS test (L6E test) |
    /// | 0x1C | MADOCA-PPP (operational L6D) |
    /// | 0x1B | MADOCALIB (L6E production) |
    #[must_use]
    pub const fn from_type_code(code: u8) -> Self {
        match code & 0x1F {
            0x1D => Self::Clas,
            0x1E => Self::ClasTest,
            0x1C => Self::MadocaPpp,
            0x1B => Self::Madocalib,
            _ => Self::Unknown,
        }
    }
}

// ---------------------------------------------------------------------------
// L6Header
// ---------------------------------------------------------------------------

/// Parsed L6 frame header.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct L6Header {
    /// Transmitting satellite PRN.
    pub prn: u8,
    /// Message type (raw 5-bit code as broadcast).
    pub message_type_code: u8,
    /// Decoded message kind.
    pub kind: L6MessageKind,
    /// Alert flag (bit 15 of header byte 2).
    pub alert_flag: bool,
    /// Vendor identifier (bits 8-15 of header byte 3).
    pub vendor_id: u8,
}

// ---------------------------------------------------------------------------
// L6Frame
// ---------------------------------------------------------------------------

/// A fully-decoded L6 frame.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct L6Frame {
    /// Parsed header.
    pub header: L6Header,
    /// Message payload (body bytes after the 4-byte header, before RS).
    pub payload: Vec<u8>,
}

// ---------------------------------------------------------------------------
// L6ParseError
// ---------------------------------------------------------------------------

/// Errors returned by [`parse_frame`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum L6ParseError {
    /// The buffer was not exactly [`FRAME_BYTES`] bytes long.
    WrongFrameLength { got: usize, expected: usize },
}

impl core::fmt::Display for L6ParseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::WrongFrameLength { got, expected } => {
                write!(f, "L6 frame length {got} != expected {expected}")
            }
        }
    }
}

impl std::error::Error for L6ParseError {}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a raw L6 frame into a header + body.
///
/// The frame is expected to be exactly [`FRAME_BYTES`] long. This function
/// does **not** verify the trailing Reed-Solomon (255,223) parity block;
/// callers should feed frames from a receiver that has already validated
/// (or repaired) the RS block.
///
/// # Errors
///
/// Returns [`L6ParseError::WrongFrameLength`] if the buffer is not
/// [`FRAME_BYTES`] bytes.
pub fn parse_frame(buf: &[u8]) -> Result<L6Frame, L6ParseError> {
    if buf.len() != FRAME_BYTES {
        return Err(L6ParseError::WrongFrameLength {
            got: buf.len(),
            expected: FRAME_BYTES,
        });
    }
    let prn = buf[0];
    let raw_type = buf[1] & 0x1F;
    let kind = L6MessageKind::from_type_code(raw_type);
    let alert_flag = (buf[2] & 0x80) != 0;
    let vendor_id = buf[3];
    let payload = buf[HEADER_BYTES..].to_vec();
    Ok(L6Frame {
        header: L6Header {
            prn,
            message_type_code: raw_type,
            kind,
            alert_flag,
            vendor_id,
        },
        payload,
    })
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_frame(prn: u8, type_code: u8, alert: bool) -> Vec<u8> {
        let mut buf = vec![0u8; FRAME_BYTES];
        buf[0] = prn;
        buf[1] = type_code & 0x1F;
        buf[2] = if alert { 0x80 } else { 0x00 };
        buf[3] = 0x42;
        buf
    }

    #[test]
    fn madoca_message_type_is_recognized() {
        let f = parse_frame(&make_frame(193, 0x1C, false)).unwrap();
        assert_eq!(f.header.kind, L6MessageKind::MadocaPpp);
    }

    #[test]
    fn clas_message_type_is_recognized() {
        let f = parse_frame(&make_frame(194, 0x1D, false)).unwrap();
        assert_eq!(f.header.kind, L6MessageKind::Clas);
    }

    #[test]
    fn madocalib_message_type_is_recognized() {
        let f = parse_frame(&make_frame(195, 0x1B, false)).unwrap();
        assert_eq!(f.header.kind, L6MessageKind::Madocalib);
    }

    #[test]
    fn clas_test_message_type_is_recognized() {
        let f = parse_frame(&make_frame(196, 0x1E, false)).unwrap();
        assert_eq!(f.header.kind, L6MessageKind::ClasTest);
    }

    #[test]
    fn unknown_type_returns_unknown() {
        let f = parse_frame(&make_frame(197, 0x01, false)).unwrap();
        assert_eq!(f.header.kind, L6MessageKind::Unknown);
    }

    #[test]
    fn alert_flag_is_captured() {
        let f = parse_frame(&make_frame(193, 0x1C, true)).unwrap();
        assert!(f.header.alert_flag);
    }

    #[test]
    fn vendor_id_is_captured() {
        let f = parse_frame(&make_frame(193, 0x1C, false)).unwrap();
        assert_eq!(f.header.vendor_id, 0x42);
    }

    #[test]
    fn payload_length_is_frame_minus_header() {
        let f = parse_frame(&make_frame(193, 0x1C, false)).unwrap();
        assert_eq!(f.payload.len(), FRAME_BYTES - HEADER_BYTES);
    }

    #[test]
    fn short_frame_is_rejected() {
        let err = parse_frame(&[0u8; 10]).unwrap_err();
        assert_eq!(
            err,
            L6ParseError::WrongFrameLength {
                got: 10,
                expected: FRAME_BYTES
            }
        );
    }

    #[test]
    fn code_labels_are_stable() {
        assert_eq!(L6MessageKind::MadocaPpp.code(), "MADOCA");
        assert_eq!(L6MessageKind::Madocalib.code(), "MADOCALIB");
        assert_eq!(L6MessageKind::Clas.code(), "CLAS");
        assert_eq!(L6MessageKind::ClasTest.code(), "CLAS-TEST");
        assert_eq!(L6MessageKind::Unknown.code(), "UNK");
    }

    #[test]
    fn error_display_is_readable() {
        let err = L6ParseError::WrongFrameLength {
            got: 100,
            expected: 250,
        };
        let s = err.to_string();
        assert!(s.contains("100"));
        assert!(s.contains("250"));
    }
}
