//! RTCM 3.3 message frame parser.
//!
//! RTCM 3.3 (Radio Technical Commission for Maritime Services Standard
//! for Differential Global Navigation Satellite Systems) is the industry
//! transport for network-RTK and PPP-RTK correction streams. Every frame
//! has the form:
//!
//! ```text
//! byte 0        preamble  = 0xD3
//! byte 1  bits 0-5 reserved (all 0)
//! byte 1-2 bits 6-15  message length (10 bits, big-endian)
//! byte 3..N  message payload (length bytes)
//! byte N..N+3  CRC-24Q parity
//! ```
//!
//! Common message types recognised here are the observation groups
//! (1001-1013), ephemeris broadcast (1019-1020, 1042, 1044-1046), MSM
//! (Multiple Signal Messages, 1071-1131) and SSR (State Space
//! Representation) types.
//!
//! # References
//!
//! - Radio Technical Commission for Maritime Services (2016),
//!   "RTCM Standard 10403.3 — Differential GNSS Services Version 3".
//! - Kaplan, E. D. & Hegarty, C. J. (2017), §8.6.

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// RTCM 3 frame preamble.
pub const RTCM3_PREAMBLE: u8 = 0xD3;
/// RTCM 3 header length (preamble + reserved + 10-bit length).
pub const RTCM3_HEADER_BYTES: usize = 3;
/// RTCM 3 CRC-24Q trailer length.
pub const RTCM3_CRC_BYTES: usize = 3;

/// Maximum RTCM 3 payload length (10 bits = 1023 bytes).
pub const MAX_PAYLOAD_BYTES: usize = 1023;

// ---------------------------------------------------------------------------
// Rtcm3MessageKind
// ---------------------------------------------------------------------------

/// Rough classification of an RTCM 3 message by its message-number.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Rtcm3MessageKind {
    /// GPS observation (1001-1004).
    GpsObservation,
    /// GLONASS observation (1009-1012).
    GlonassObservation,
    /// Reference station coordinates (1005-1008, 1013, 1032, 1033).
    StationInfo,
    /// GPS broadcast ephemeris (1019).
    GpsEphemeris,
    /// GLONASS broadcast ephemeris (1020).
    GlonassEphemeris,
    /// Galileo F/NAV ephemeris (1042, 1045).
    GalileoEphemeris,
    /// BeiDou ephemeris (1042, 1046 — legacy overlap).
    BeidouEphemeris,
    /// QZSS ephemeris (1044).
    QzssEphemeris,
    /// SBAS ephemeris (1043).
    SbasEphemeris,
    /// MSM (Multiple Signal Messages, 1071-1077 GPS / 1081-1087 GLO / ...).
    MultipleSignalMessages,
    /// SSR (State Space Representation, 1057-1068 / 1240-1270).
    StateSpaceRepresentation,
    /// Proprietary / vendor-specific (4001-4095).
    Proprietary,
    /// Any other or unknown type.
    Other,
}

impl Rtcm3MessageKind {
    /// Classify by RTCM 3 message number.
    #[must_use]
    pub const fn from_msg_no(msg_no: u16) -> Self {
        match msg_no {
            1001..=1004 => Self::GpsObservation,
            1005..=1008 | 1013 | 1032 | 1033 => Self::StationInfo,
            1009..=1012 => Self::GlonassObservation,
            1019 => Self::GpsEphemeris,
            1020 => Self::GlonassEphemeris,
            1042 | 1045 => Self::GalileoEphemeris,
            1046 => Self::BeidouEphemeris,
            1044 => Self::QzssEphemeris,
            1043 => Self::SbasEphemeris,
            1071..=1137 => Self::MultipleSignalMessages,
            1057..=1068 | 1240..=1270 => Self::StateSpaceRepresentation,
            4001..=4095 => Self::Proprietary,
            _ => Self::Other,
        }
    }
}

// ---------------------------------------------------------------------------
// Rtcm3Frame
// ---------------------------------------------------------------------------

/// A parsed RTCM 3 frame.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct Rtcm3Frame {
    /// Payload length declared by the frame header.
    pub payload_len: u16,
    /// 12-bit message number (from bits 0-11 of the payload).
    pub message_no: u16,
    /// Classified message kind.
    pub kind: Rtcm3MessageKind,
    /// Payload bytes as broadcast (message-number stays in bits 0-11).
    pub payload: Vec<u8>,
    /// CRC-24Q value read from the trailer.
    pub crc24_declared: u32,
    /// CRC-24Q value recomputed over preamble + length + payload.
    pub crc24_computed: u32,
}

impl Rtcm3Frame {
    /// Whether the declared CRC-24Q matches the recomputed one.
    #[must_use]
    pub const fn crc_valid(&self) -> bool {
        self.crc24_declared == self.crc24_computed
    }
}

// ---------------------------------------------------------------------------
// Rtcm3ParseError
// ---------------------------------------------------------------------------

/// Errors returned by [`parse_frame`].
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Rtcm3ParseError {
    /// The buffer was too short to contain a header.
    TooShort { got: usize },
    /// The buffer did not start with the RTCM 3 preamble.
    BadPreamble { got: u8 },
    /// The high 6 bits of byte 1 were not zero (RTCM 3 reserved).
    NonZeroReserved,
    /// The declared payload length exceeded the buffer.
    LengthExceedsBuffer { declared: u16, available: usize },
}

impl core::fmt::Display for Rtcm3ParseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::TooShort { got } => write!(f, "RTCM3 buffer too short ({got} bytes)"),
            Self::BadPreamble { got } => write!(f, "RTCM3 bad preamble 0x{got:02X}"),
            Self::NonZeroReserved => write!(f, "RTCM3 reserved bits non-zero"),
            Self::LengthExceedsBuffer {
                declared,
                available,
            } => write!(
                f,
                "RTCM3 declared payload length {declared} exceeds available {available}"
            ),
        }
    }
}

impl std::error::Error for Rtcm3ParseError {}

// ---------------------------------------------------------------------------
// Parser
// ---------------------------------------------------------------------------

/// Parse a single RTCM 3 frame from the front of `buf`.
///
/// # Errors
///
/// Returns [`Rtcm3ParseError`] variants when the buffer is malformed.
pub fn parse_frame(buf: &[u8]) -> Result<Rtcm3Frame, Rtcm3ParseError> {
    if buf.len() < RTCM3_HEADER_BYTES + RTCM3_CRC_BYTES {
        return Err(Rtcm3ParseError::TooShort { got: buf.len() });
    }
    if buf[0] != RTCM3_PREAMBLE {
        return Err(Rtcm3ParseError::BadPreamble { got: buf[0] });
    }
    if (buf[1] & 0xFC) != 0 {
        return Err(Rtcm3ParseError::NonZeroReserved);
    }
    let len_high = u16::from(buf[1] & 0x03);
    let len_low = u16::from(buf[2]);
    let payload_len = (len_high << 8) | len_low;
    let payload_end = RTCM3_HEADER_BYTES + payload_len as usize;
    let frame_end = payload_end + RTCM3_CRC_BYTES;
    if frame_end > buf.len() {
        return Err(Rtcm3ParseError::LengthExceedsBuffer {
            declared: payload_len,
            available: buf.len(),
        });
    }
    let payload = buf[RTCM3_HEADER_BYTES..payload_end].to_vec();
    let message_no = if payload.len() >= 2 {
        let mn_high = u16::from(payload[0]);
        let mn_low = u16::from(payload[1] & 0xF0);
        (mn_high << 4) | (mn_low >> 4)
    } else {
        0
    };
    let kind = Rtcm3MessageKind::from_msg_no(message_no);

    let crc_start = payload_end;
    let crc24_declared = (u32::from(buf[crc_start]) << 16)
        | (u32::from(buf[crc_start + 1]) << 8)
        | u32::from(buf[crc_start + 2]);
    let crc24_computed = crc24_q(&buf[..payload_end]);

    Ok(Rtcm3Frame {
        payload_len,
        message_no,
        kind,
        payload,
        crc24_declared,
        crc24_computed,
    })
}

// ---------------------------------------------------------------------------
// CRC-24Q
// ---------------------------------------------------------------------------

/// CRC-24Q polynomial used by RTCM 3 (Qualcomm 24-bit CRC).
///
/// `Q(x) = x^24 + x^23 + x^18 + x^17 + x^14 + x^11 + x^10 + x^7 +
///         x^6 + x^5 + x^4 + x^3 + x + 1`.
///
/// Initial value 0, no output XOR, byte-oriented left-shift processing.
#[must_use]
pub fn crc24_q(data: &[u8]) -> u32 {
    const POLY: u32 = 0x0186_4CFB;
    let mut crc: u32 = 0;
    for &b in data {
        crc ^= u32::from(b) << 16;
        for _ in 0..8 {
            crc <<= 1;
            if crc & 0x0100_0000 != 0 {
                crc ^= POLY;
            }
        }
        crc &= 0x00FF_FFFF;
    }
    crc
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn make_frame(msg_no: u16, payload_body: &[u8]) -> Vec<u8> {
        let mut payload = Vec::new();
        payload.push((msg_no >> 4) as u8);
        payload.push(((msg_no & 0x000F) << 4) as u8);
        payload.extend_from_slice(payload_body);
        let mut frame = Vec::with_capacity(RTCM3_HEADER_BYTES + payload.len() + RTCM3_CRC_BYTES);
        let len = payload.len() as u16;
        frame.push(RTCM3_PREAMBLE);
        frame.push(((len >> 8) & 0x03) as u8);
        frame.push((len & 0xFF) as u8);
        frame.extend_from_slice(&payload);
        let crc = crc24_q(&frame);
        frame.push(((crc >> 16) & 0xFF) as u8);
        frame.push(((crc >> 8) & 0xFF) as u8);
        frame.push((crc & 0xFF) as u8);
        frame
    }

    #[test]
    fn parse_msm_type_1077() {
        let f = parse_frame(&make_frame(1077, &[0u8; 20])).unwrap();
        assert_eq!(f.message_no, 1077);
        assert_eq!(f.kind, Rtcm3MessageKind::MultipleSignalMessages);
    }

    #[test]
    fn parse_gps_ephemeris_1019() {
        let f = parse_frame(&make_frame(1019, &[0u8; 60])).unwrap();
        assert_eq!(f.kind, Rtcm3MessageKind::GpsEphemeris);
    }

    #[test]
    fn parse_qzss_ephemeris_1044() {
        let f = parse_frame(&make_frame(1044, &[0u8; 60])).unwrap();
        assert_eq!(f.kind, Rtcm3MessageKind::QzssEphemeris);
    }

    #[test]
    fn parse_station_info_1005() {
        let f = parse_frame(&make_frame(1005, &[0u8; 18])).unwrap();
        assert_eq!(f.kind, Rtcm3MessageKind::StationInfo);
    }

    #[test]
    fn parse_ssr_type_1057() {
        let f = parse_frame(&make_frame(1057, &[0u8; 12])).unwrap();
        assert_eq!(f.kind, Rtcm3MessageKind::StateSpaceRepresentation);
    }

    #[test]
    fn parse_proprietary_type_4062() {
        let f = parse_frame(&make_frame(4062, &[0u8; 4])).unwrap();
        assert_eq!(f.kind, Rtcm3MessageKind::Proprietary);
    }

    #[test]
    fn crc_matches_when_frame_is_intact() {
        let f = parse_frame(&make_frame(1077, &[0u8; 20])).unwrap();
        assert!(f.crc_valid());
    }

    #[test]
    fn crc_fails_when_frame_is_corrupted() {
        let mut frame = make_frame(1077, &[0u8; 20]);
        let idx = frame.len() - 5;
        frame[idx] ^= 0xFF;
        let f = parse_frame(&frame).unwrap();
        assert!(!f.crc_valid());
    }

    #[test]
    fn bad_preamble_is_rejected() {
        let mut frame = make_frame(1077, &[0u8; 4]);
        frame[0] = 0xFF;
        let err = parse_frame(&frame).unwrap_err();
        assert_eq!(err, Rtcm3ParseError::BadPreamble { got: 0xFF });
    }

    #[test]
    fn short_buffer_is_rejected() {
        let err = parse_frame(&[0xD3, 0x00]).unwrap_err();
        assert_eq!(err, Rtcm3ParseError::TooShort { got: 2 });
    }

    #[test]
    fn non_zero_reserved_bits_are_rejected() {
        let mut frame = make_frame(1077, &[0u8; 4]);
        frame[1] |= 0x80;
        let err = parse_frame(&frame).unwrap_err();
        assert_eq!(err, Rtcm3ParseError::NonZeroReserved);
    }

    #[test]
    fn length_beyond_buffer_is_rejected() {
        let buf = vec![0xD3, 0x03, 0xFF, 0x00, 0x00, 0x00];
        let err = parse_frame(&buf).unwrap_err();
        assert!(matches!(err, Rtcm3ParseError::LengthExceedsBuffer { .. }));
    }

    #[test]
    fn crc24q_of_empty_input_is_zero() {
        assert_eq!(crc24_q(&[]), 0);
    }

    #[test]
    fn crc24q_of_single_zero_byte_is_zero() {
        assert_eq!(crc24_q(&[0x00]), 0);
    }

    #[test]
    fn error_display_readable() {
        assert!(Rtcm3ParseError::BadPreamble { got: 0xAB }
            .to_string()
            .contains("AB"));
    }
}
