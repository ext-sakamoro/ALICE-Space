//! Deep-space communication protocol — model differentials over minimal bandwidth.

use crate::fnv1a;
use crate::orbit::BodyId;

/// Communication link between two bodies.
#[derive(Debug, Clone)]
pub struct CommLink {
    pub source_id: BodyId,
    pub target_id: BodyId,
    pub distance_km: f64,
    pub bandwidth_bps: f64,
}

impl CommLink {
    pub fn new(source: u64, target: u64, distance_km: f64, bandwidth_bps: f64) -> Self {
        Self {
            source_id: BodyId(source),
            target_id: BodyId(target),
            distance_km,
            bandwidth_bps,
        }
    }

    /// One-way latency in seconds.
    #[inline]
    pub fn latency_s(&self) -> f64 {
        self.distance_km / 299792.458
    }

    /// Total bits available in a transmission window.
    #[inline]
    pub fn bits_per_window(&self, window_s: f64) -> f64 {
        self.bandwidth_bps * window_s
    }
}

/// A compact model update sent over deep-space link.
#[derive(Debug, Clone)]
pub struct ModelDifferential {
    pub sequence: u64,
    pub timestamp_ns: u64,
    /// Parameter updates: (parameter_name_hash, new_value).
    pub param_updates: Vec<(u64, f64)>,
    pub content_hash: u64,
}

impl ModelDifferential {
    pub fn new(sequence: u64, timestamp_ns: u64) -> Self {
        Self { sequence, timestamp_ns, param_updates: Vec::new(), content_hash: 0 }
    }

    /// Add a parameter update (name is hashed via FNV-1a).
    pub fn add_param(&mut self, name: &str, value: f64) {
        self.param_updates.push((fnv1a(name.as_bytes()), value));
    }

    /// Estimated wire size in bytes.
    pub fn byte_size(&self) -> usize {
        8 + 8 + 8 + self.param_updates.len() * 16 // seq + ts + hash + N*(hash+f64)
    }

    /// Finalize content_hash over all parameter data.
    pub fn finalize(&mut self) {
        let mut data = Vec::with_capacity(16 + self.param_updates.len() * 16);
        data.extend_from_slice(&self.sequence.to_le_bytes());
        data.extend_from_slice(&self.timestamp_ns.to_le_bytes());
        for &(h, v) in &self.param_updates {
            data.extend_from_slice(&h.to_le_bytes());
            data.extend_from_slice(&v.to_bits().to_le_bytes());
        }
        self.content_hash = fnv1a(&data);
    }
}

/// Check if a differential can be transmitted within a given window.
pub fn can_transmit(diff: &ModelDifferential, link: &CommLink, window_s: f64) -> bool {
    (diff.byte_size() * 8) as f64 <= link.bits_per_window(window_s)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn link_latency() {
        let link = CommLink::new(1, 2, 299792.458, 1000.0);
        assert!((link.latency_s() - 1.0).abs() < 1e-6);
    }

    #[test]
    fn bits_per_window() {
        let link = CommLink::new(1, 2, 1000.0, 9600.0);
        assert!((link.bits_per_window(10.0) - 96000.0).abs() < 1e-6);
    }

    #[test]
    fn differential_byte_size() {
        let mut diff = ModelDifferential::new(1, 1000);
        assert_eq!(diff.byte_size(), 24); // 8+8+8+0
        diff.add_param("x", 1.0);
        assert_eq!(diff.byte_size(), 40); // 24+16
    }

    #[test]
    fn can_transmit_check() {
        let mut diff = ModelDifferential::new(1, 1000);
        diff.add_param("x", 1.0);
        diff.add_param("y", 2.0);
        // 56 bytes = 448 bits, need >= 448 bps in 1s window
        let link = CommLink::new(1, 2, 1000.0, 500.0);
        assert!(can_transmit(&diff, &link, 1.0));
        let slow_link = CommLink::new(1, 2, 1000.0, 100.0);
        assert!(!can_transmit(&diff, &slow_link, 1.0));
    }

    #[test]
    fn finalize_deterministic() {
        let mut d1 = ModelDifferential::new(1, 1000);
        d1.add_param("thrust_x", 3.14);
        d1.finalize();
        let mut d2 = ModelDifferential::new(1, 1000);
        d2.add_param("thrust_x", 3.14);
        d2.finalize();
        assert_eq!(d1.content_hash, d2.content_hash);
        assert_ne!(d1.content_hash, 0);
    }

    #[test]
    fn finalize_different_params_different_hash() {
        let mut d1 = ModelDifferential::new(1, 1000);
        d1.add_param("thrust_x", 3.14);
        d1.finalize();
        let mut d2 = ModelDifferential::new(1, 1000);
        d2.add_param("thrust_y", 3.14);
        d2.finalize();
        assert_ne!(d1.content_hash, d2.content_hash);
    }

    #[test]
    fn finalize_different_values_different_hash() {
        let mut d1 = ModelDifferential::new(1, 1000);
        d1.add_param("thrust_x", 1.0);
        d1.finalize();
        let mut d2 = ModelDifferential::new(1, 1000);
        d2.add_param("thrust_x", 2.0);
        d2.finalize();
        assert_ne!(d1.content_hash, d2.content_hash);
    }

    #[test]
    fn empty_differential_byte_size() {
        let diff = ModelDifferential::new(0, 0);
        assert_eq!(diff.byte_size(), 24);
    }

    #[test]
    fn differential_multiple_params_byte_size() {
        let mut diff = ModelDifferential::new(1, 100);
        for i in 0usize..10 {
            diff.add_param(&format!("p{i}"), i as f64);
        }
        // 24 base + 10 * 16 = 184
        assert_eq!(diff.byte_size(), 184);
    }

    #[test]
    fn link_latency_deep_space() {
        // Mars average distance: ~225 million km → ~750 s
        let link = CommLink::new(1, 2, 225_000_000.0, 1000.0);
        let latency = link.latency_s();
        assert!((latency - 750.5).abs() < 2.0, "Mars latency: {} s", latency);
    }

    #[test]
    fn can_transmit_exact_boundary() {
        // Exact boundary: diff size = link capacity
        let mut diff = ModelDifferential::new(1, 100);
        diff.add_param("x", 1.0);
        // 40 bytes = 320 bits; need exactly 320 bps in 1s window
        let link = CommLink::new(1, 2, 1000.0, 320.0);
        assert!(can_transmit(&diff, &link, 1.0));
    }

    #[test]
    fn can_transmit_just_below_boundary() {
        let mut diff = ModelDifferential::new(1, 100);
        diff.add_param("x", 1.0);
        // 40 bytes = 320 bits; 319 bps not enough in 1s
        let link = CommLink::new(1, 2, 1000.0, 319.0);
        assert!(!can_transmit(&diff, &link, 1.0));
    }

    #[test]
    fn comm_link_source_target_ids() {
        let link = CommLink::new(10, 20, 500.0, 9600.0);
        assert_eq!(link.source_id, BodyId(10));
        assert_eq!(link.target_id, BodyId(20));
    }
}
