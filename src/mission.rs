//! Mission planning and event logging.

#[inline(always)]
fn fnv1a(data: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf29ce484222325;
    for &b in data { h ^= b as u64; h = h.wrapping_mul(0x100000001b3); }
    h
}

/// Mission phase classification.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MissionPhase {
    Launch,
    TransferOrbit,
    Insertion,
    Orbiting,
    Landing,
    Surface,
    Ascent,
    Return,
}

/// A recorded mission event.
#[derive(Debug, Clone)]
pub struct MissionEvent {
    pub sequence: u32,
    pub phase: MissionPhase,
    pub timestamp_ns: u64,
    pub delta_v_used: f64,
    pub fuel_remaining_kg: f64,
    pub content_hash: u64,
}

/// Append-only mission event log.
#[derive(Debug, Clone)]
pub struct MissionLog {
    events: Vec<MissionEvent>,
    next_sequence: u32,
}

impl MissionLog {
    pub fn new() -> Self {
        Self { events: Vec::new(), next_sequence: 0 }
    }

    /// Log a mission event, returns sequence number.
    pub fn log_event(&mut self, phase: MissionPhase, timestamp_ns: u64, delta_v: f64, fuel_kg: f64) -> u32 {
        let seq = self.next_sequence;
        let mut hash_data = [0u8; 25];
        hash_data[0..4].copy_from_slice(&seq.to_le_bytes());
        hash_data[4] = phase as u8;
        hash_data[5..13].copy_from_slice(&timestamp_ns.to_le_bytes());
        hash_data[13..21].copy_from_slice(&delta_v.to_bits().to_le_bytes());
        hash_data[21..25].copy_from_slice(&(fuel_kg as u32).to_le_bytes());

        self.events.push(MissionEvent {
            sequence: seq,
            phase,
            timestamp_ns,
            delta_v_used: delta_v,
            fuel_remaining_kg: fuel_kg,
            content_hash: fnv1a(&hash_data),
        });
        self.next_sequence += 1;
        seq
    }

    /// Total delta-v expended across all events.
    pub fn total_delta_v(&self) -> f64 {
        self.events.iter().map(|e| e.delta_v_used).sum()
    }

    /// Current (most recent) mission phase.
    pub fn current_phase(&self) -> Option<&MissionPhase> {
        self.events.last().map(|e| &e.phase)
    }

    /// All events for a specific phase.
    pub fn events_for_phase(&self, phase: MissionPhase) -> Vec<&MissionEvent> {
        self.events.iter().filter(|e| e.phase == phase).collect()
    }

    /// Average delta-v per event (fuel efficiency metric).
    pub fn fuel_efficiency(&self) -> f64 {
        if self.events.is_empty() { return 0.0; }
        self.total_delta_v() / self.events.len() as f64
    }

    /// Number of logged events.
    pub fn len(&self) -> usize {
        self.events.len()
    }

    /// Whether the log is empty.
    pub fn is_empty(&self) -> bool {
        self.events.is_empty()
    }
}

impl Default for MissionLog {
    fn default() -> Self { Self::new() }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn log_events() {
        let mut log = MissionLog::new();
        let seq = log.log_event(MissionPhase::Launch, 1000, 2.5, 500.0);
        assert_eq!(seq, 0);
        assert_eq!(log.len(), 1);
    }

    #[test]
    fn total_delta_v() {
        let mut log = MissionLog::new();
        log.log_event(MissionPhase::Launch, 1000, 2.5, 500.0);
        log.log_event(MissionPhase::TransferOrbit, 2000, 1.5, 450.0);
        assert!((log.total_delta_v() - 4.0).abs() < 1e-10);
    }

    #[test]
    fn current_phase() {
        let mut log = MissionLog::new();
        log.log_event(MissionPhase::Launch, 1000, 2.5, 500.0);
        log.log_event(MissionPhase::Orbiting, 2000, 0.0, 480.0);
        assert_eq!(log.current_phase(), Some(&MissionPhase::Orbiting));
    }

    #[test]
    fn events_for_phase_filter() {
        let mut log = MissionLog::new();
        log.log_event(MissionPhase::Launch, 1000, 2.5, 500.0);
        log.log_event(MissionPhase::Orbiting, 2000, 0.1, 490.0);
        log.log_event(MissionPhase::Orbiting, 3000, 0.2, 480.0);
        let orbiting = log.events_for_phase(MissionPhase::Orbiting);
        assert_eq!(orbiting.len(), 2);
    }

    #[test]
    fn fuel_efficiency() {
        let mut log = MissionLog::new();
        log.log_event(MissionPhase::Launch, 1000, 3.0, 500.0);
        log.log_event(MissionPhase::TransferOrbit, 2000, 1.0, 450.0);
        assert!((log.fuel_efficiency() - 2.0).abs() < 1e-10);
    }

    #[test]
    fn empty_log() {
        let log = MissionLog::new();
        assert!(log.is_empty());
        assert_eq!(log.total_delta_v(), 0.0);
        assert_eq!(log.fuel_efficiency(), 0.0);
        assert!(log.current_phase().is_none());
    }
}
