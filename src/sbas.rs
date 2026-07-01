//! `SBAS` message-type skeleton (MOPS `DO-229`, `RTCA` / `EUROCAE ED-259A`).
//!
//! `SBAS` (Satellite-Based Augmentation System) transmits `250-bit` messages
//! at `250 bps` from geostationary satellites. Japan's `MSAS`, the U.S. `WAAS`,
//! Europe's `EGNOS` and India's `GAGAN` all share the same message envelope:
//!
//! - Preamble byte + 6-bit message type ID + 218-bit payload +
//!   24-bit `CRC-24Q`.
//!
//! This module carries the message-type enumeration, a `SbasMessage` struct
//! that keeps the raw payload alongside decoded metadata, and a decoder for
//! the two most operationally important types:
//!
//! - Type 1 — PRN mask assignments (which satellites the corrections apply to).
//! - Type 2 – 5 — fast corrections for pseudorange bias.
//!
//! Ionospheric, long-term and integrity messages are represented but their
//! payloads are surfaced verbatim so downstream code can extend the decoder.

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------

/// Bit length of one `SBAS` message.
pub const SBAS_MESSAGE_BITS: usize = 250;
/// Length in bytes when a message is packed into a byte buffer.
pub const SBAS_MESSAGE_BYTES: usize = SBAS_MESSAGE_BITS.div_ceil(8);
/// Sync-preamble bytes used by MOPS `DO-229`.
pub const SBAS_PREAMBLES: [u8; 3] = [0x53, 0x9A, 0xC6];

// ---------------------------------------------------------------------------
// Message type
// ---------------------------------------------------------------------------

/// `SBAS` message-type identifier (6-bit field, values 0-63).
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum SbasMessageType {
    /// Type 0 — reserved for testing / do-not-use.
    DoNotUse,
    /// Type 1 — PRN mask assignments.
    PrnMask,
    /// Type 2-5 — fast corrections.
    FastCorrection(u8),
    /// Type 6 — integrity information.
    Integrity,
    /// Type 7 — fast correction degradation factor.
    FastCorrectionDegradation,
    /// Type 9 — GEO navigation message.
    GeoNavigation,
    /// Type 10 — degradation parameters.
    Degradation,
    /// Type 17 — GEO satellite almanacs.
    GeoAlmanac,
    /// Type 18 — ionospheric grid points mask.
    IonosphericMask,
    /// Type 24 — mixed fast / long-term correction.
    MixedCorrection,
    /// Type 25 — long-term satellite error corrections.
    LongTermCorrection,
    /// Type 26 — ionospheric delay corrections.
    IonosphericDelay,
    /// Type 27 — SBAS service message.
    ServiceMessage,
    /// Type 28 — clock ephemeris covariance matrix.
    ClockEphemerisCovariance,
    /// Type 63 — null message.
    Null,
    /// Anything else, kept as the raw 6-bit value.
    Unknown(u8),
}

impl SbasMessageType {
    /// Map a raw 6-bit identifier to the enum.
    #[must_use]
    pub const fn from_raw(id: u8) -> Self {
        match id {
            0 => Self::DoNotUse,
            1 => Self::PrnMask,
            2..=5 => Self::FastCorrection(id),
            6 => Self::Integrity,
            7 => Self::FastCorrectionDegradation,
            9 => Self::GeoNavigation,
            10 => Self::Degradation,
            17 => Self::GeoAlmanac,
            18 => Self::IonosphericMask,
            24 => Self::MixedCorrection,
            25 => Self::LongTermCorrection,
            26 => Self::IonosphericDelay,
            27 => Self::ServiceMessage,
            28 => Self::ClockEphemerisCovariance,
            63 => Self::Null,
            other => Self::Unknown(other),
        }
    }
}

// ---------------------------------------------------------------------------
// SbasMessage
// ---------------------------------------------------------------------------

/// One decoded (or partially decoded) `SBAS` message.
#[derive(Debug, Clone)]
pub struct SbasMessage {
    pub preamble: u8,
    pub message_type: SbasMessageType,
    /// Raw 218-bit payload padded to whole bytes (27 bytes actually used;
    /// the trailing 2 bits of the last byte are set to zero).
    pub payload: [u8; 27],
    /// Advertised CRC (24 bits, stored right-aligned in the lower 24 bits).
    pub crc24: u32,
}

impl SbasMessage {
    /// Extract a message from a packed 250-bit buffer.
    ///
    /// The first byte is the preamble; the next 6 bits carry the message
    /// type; the remaining 218 bits are the payload followed by the CRC.
    /// Returns `None` when the buffer is shorter than [`SBAS_MESSAGE_BYTES`].
    #[must_use]
    pub fn from_bytes(bytes: &[u8]) -> Option<Self> {
        if bytes.len() < SBAS_MESSAGE_BYTES {
            return None;
        }
        let preamble = bytes[0];
        // Message type: top 6 bits of byte 1.
        let raw_type = bytes[1] >> 2;
        let message_type = SbasMessageType::from_raw(raw_type);

        let mut payload = [0u8; 27];
        // Copy bits 6..224 of the message (26 full bytes + 6 bits of byte 28).
        payload[..26].copy_from_slice(&bytes[2..28]);
        payload[26] = bytes[28] & 0b1111_1100;

        // CRC24 = last 24 bits.
        let crc24 = (u32::from(bytes[28] & 0b0000_0011) << 22)
            | (u32::from(bytes[29]) << 14)
            | (u32::from(bytes[30]) << 6)
            | (u32::from(bytes[31]) >> 2);

        Some(Self {
            preamble,
            message_type,
            payload,
            crc24,
        })
    }
}

// ---------------------------------------------------------------------------
// Type 1 — PRN mask
// ---------------------------------------------------------------------------

/// Decoded Type 1 PRN mask.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct PrnMask {
    /// The bit vector: `mask[i]` is `true` when slot `i + 1` is included.
    /// The MOPS mask is 210 bits long; excess bytes are ignored.
    pub bits: [bool; 210],
    /// Issue-of-data (IODP) counter, 2 bits.
    pub iodp: u8,
}

impl PrnMask {
    /// Decode a Type 1 payload.
    ///
    /// # Panics
    ///
    /// Panics if `payload` is shorter than 27 bytes.
    #[must_use]
    pub fn decode(payload: &[u8; 27]) -> Self {
        let mut bits = [false; 210];
        // The first 210 bits are the mask.
        for i in 0..210 {
            let byte = payload[i / 8];
            let bit = 7 - (i % 8);
            bits[i] = (byte >> bit) & 1 == 1;
        }
        // IODP is at bit positions 210 and 211.
        let iodp_byte = payload[26];
        let iodp = (iodp_byte >> 4) & 0b11;
        Self { bits, iodp }
    }

    /// Count of active satellites in the mask.
    #[must_use]
    pub fn active_count(&self) -> usize {
        self.bits.iter().filter(|b| **b).count()
    }
}

// ---------------------------------------------------------------------------
// Type 2-5 — fast corrections
// ---------------------------------------------------------------------------

/// Decoded fast correction message. Only the first `N` slots (matched to a
/// `PrnMask`) carry meaningful values; the remainder are zero-filled.
#[derive(Debug, Clone, PartialEq)]
pub struct FastCorrection {
    /// Slot index of the first correction inside this message (0, 6, 12 …).
    pub slot_offset: u8,
    /// Pseudorange correction (metres) per active slot, 13 entries.
    pub prc_m: [f64; 13],
    /// User Differential Range Error index per active slot, 13 entries.
    pub udrei: [u8; 13],
    /// Issue-of-data-P (IODP), 2 bits.
    pub iodp: u8,
}

impl FastCorrection {
    /// Decode a Type 2-5 payload.
    ///
    /// # Panics
    ///
    /// Panics if `payload` is shorter than 27 bytes.
    #[must_use]
    pub fn decode(msg_type: u8, payload: &[u8; 27]) -> Self {
        // Type 2 -> slot_offset 0; Type 3 -> 13; Type 4 -> 26; Type 5 -> 39.
        let slot_offset = (msg_type - 2) * 13;
        let mut prc_m = [0.0_f64; 13];
        let mut udrei = [0_u8; 13];
        // Each PRC is a 12-bit two's-complement value; 13 PRCs = 156 bits.
        // Each UDREI is 4 bits; 13 UDREIs = 52 bits.
        // Followed by 2-bit IODP.
        for i in 0..13 {
            let bit_start = i * 12;
            let raw = extract_bits(payload, bit_start, 12) as i16;
            // Sign-extend from 12 bits.
            let extended = if raw & 0x800 != 0 {
                (raw as i32) | !0xFFF
            } else {
                raw as i32
            };
            // Resolution: 0.125 m.
            prc_m[i] = extended as f64 * 0.125;
        }
        for i in 0..13 {
            let bit_start = 156 + i * 4;
            udrei[i] = extract_bits(payload, bit_start, 4) as u8;
        }
        let iodp = extract_bits(payload, 208, 2) as u8;
        Self {
            slot_offset,
            prc_m,
            udrei,
            iodp,
        }
    }
}

// ---------------------------------------------------------------------------
// Bit helpers
// ---------------------------------------------------------------------------

fn extract_bits(payload: &[u8; 27], start_bit: usize, width: usize) -> u32 {
    let mut out: u32 = 0;
    for i in 0..width {
        let bit_index = start_bit + i;
        let byte = payload[bit_index / 8];
        let bit_val = (byte >> (7 - (bit_index % 8))) & 1;
        out = (out << 1) | u32::from(bit_val);
    }
    out
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn message_type_maps_common_values() {
        assert_eq!(SbasMessageType::from_raw(0), SbasMessageType::DoNotUse);
        assert_eq!(SbasMessageType::from_raw(1), SbasMessageType::PrnMask);
        assert_eq!(
            SbasMessageType::from_raw(3),
            SbasMessageType::FastCorrection(3)
        );
        assert_eq!(SbasMessageType::from_raw(6), SbasMessageType::Integrity);
        assert_eq!(SbasMessageType::from_raw(63), SbasMessageType::Null);
    }

    #[test]
    fn message_type_unknown_preserves_id() {
        assert_eq!(SbasMessageType::from_raw(45), SbasMessageType::Unknown(45));
    }

    #[test]
    fn message_from_bytes_needs_full_length() {
        let short = vec![0u8; 5];
        assert!(SbasMessage::from_bytes(&short).is_none());
        let full = vec![0u8; SBAS_MESSAGE_BYTES];
        assert!(SbasMessage::from_bytes(&full).is_some());
    }

    #[test]
    fn message_from_bytes_captures_preamble_and_type() {
        let mut buf = [0u8; SBAS_MESSAGE_BYTES];
        buf[0] = SBAS_PREAMBLES[0];
        // Type 1 << 2 = 4.
        buf[1] = 1 << 2;
        let msg = SbasMessage::from_bytes(&buf).unwrap();
        assert_eq!(msg.preamble, SBAS_PREAMBLES[0]);
        assert_eq!(msg.message_type, SbasMessageType::PrnMask);
    }

    #[test]
    fn prn_mask_decodes_active_bits() {
        // Craft a payload with bits 0, 5, 10 set.
        let mut payload = [0u8; 27];
        payload[0] = 0b1000_0100;
        payload[1] = 0b0010_0000;
        let mask = PrnMask::decode(&payload);
        assert!(mask.bits[0]);
        assert!(mask.bits[5]);
        assert!(mask.bits[10]);
        assert_eq!(mask.active_count(), 3);
    }

    #[test]
    fn prn_mask_iodp_decoded() {
        let mut payload = [0u8; 27];
        payload[26] = 0b0010_0000;
        let mask = PrnMask::decode(&payload);
        assert_eq!(mask.iodp, 0b10);
    }

    #[test]
    fn fast_correction_slot_offset_scales_with_type() {
        let payload = [0u8; 27];
        assert_eq!(FastCorrection::decode(2, &payload).slot_offset, 0);
        assert_eq!(FastCorrection::decode(3, &payload).slot_offset, 13);
        assert_eq!(FastCorrection::decode(5, &payload).slot_offset, 39);
    }

    #[test]
    fn fast_correction_decodes_zero_payload() {
        let payload = [0u8; 27];
        let corr = FastCorrection::decode(2, &payload);
        assert!(corr.prc_m.iter().all(|v| v.abs() < 1e-12));
        assert!(corr.udrei.iter().all(|v| *v == 0));
    }

    #[test]
    fn extract_bits_reads_high_order_first() {
        let mut payload = [0u8; 27];
        payload[0] = 0b1010_1100;
        assert_eq!(extract_bits(&payload, 0, 4), 0b1010);
        assert_eq!(extract_bits(&payload, 4, 4), 0b1100);
    }
}
