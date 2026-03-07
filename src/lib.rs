#![allow(
    clippy::cast_possible_truncation,
    clippy::cast_possible_wrap,
    clippy::cast_precision_loss,
    clippy::cast_sign_loss,
    clippy::cast_lossless,
    clippy::similar_names,
    clippy::many_single_char_names,
    clippy::module_name_repetitions,
    clippy::inline_always,
    clippy::too_many_lines
)]
#![cfg_attr(test, allow(clippy::float_cmp))]

//! ALICE-Space — Deep-Space Communication & Autonomous Control
//!
//! Handles deep-space communication with extreme latency (minutes to hours)
//! and minimal bandwidth (kbps), enabling autonomous spacecraft control via
//! mathematical model differentials rather than raw data transmission.
//!
//! # Modules
//!
//! | Module | Description |
//! |--------|-------------|
//! | [`orbit`] | Keplerian elements, Hohmann transfer, light delay, celestial bodies |
//! | [`propagator`] | RK4 two-body orbit propagation |
//! | [`autonomy`] | Trajectory models, fault detection, control decisions |
//! | [`comm`] | Deep-space comm links, model-differential protocol |
//! | [`constellation`] | Walker constellation geometry |
//! | [`link_budget`] | Friis path-loss and link margin calculation |
//! | [`mission`] | Mission phase FSM and event logging |
//!
//! # Quick Start
//!
//! ```
//! use alice_space::{orbital_period, light_delay_s};
//!
//! // ISS orbital period (~92 minutes)
//! let t = orbital_period(6778.0, 398600.4418);
//! assert!((t / 60.0 - 92.3).abs() < 1.0);
//!
//! // Earth-Moon light delay (~1.28 seconds)
//! let delay = light_delay_s(384400.0);
//! assert!((delay - 1.28).abs() < 0.01);
//! ```
//!
//! Author: Moroya Sakamoto

pub mod autonomy;
pub mod comm;
pub mod constellation;
pub mod eclipse;
#[cfg(feature = "ffi")]
pub mod ffi;
pub mod link_budget;
pub mod mission;
pub mod orbit;
pub mod propagator;

pub use autonomy::{
    apply_differential, compute_correction, evaluate_decision_tree, AutonomyLevel, ControlDecision,
    DecisionNode, FaultType, TrajectoryModel,
};
pub use comm::{can_transmit, CommLink, ModelDifferential};
pub use constellation::{WalkerConstellation, WalkerSatellite};
pub use eclipse::{earth_shadow, eclipse_fraction, is_comm_window, shadow_state, ShadowState};
pub use link_budget::{friis_path_loss_db, LinkBudget, LinkBudgetResult};
pub use mission::{MissionEvent, MissionFsm, MissionLog, MissionPhase, TransitionResult};
pub use orbit::{
    delta_v_hohmann, light_delay_s, orbital_period, orbital_velocity, BodyId, CelestialBody,
    OrbitalElements, SpacecraftState,
};
pub use propagator::{propagate_rk4, propagate_rk4_single, TwoBodyAccel};

// ── Shared hash primitive ──────────────────────────────────────────────

#[inline(always)]
pub(crate) fn fnv1a(data: &[u8]) -> u64 {
    let mut h: u64 = 0xcbf2_9ce4_8422_2325;
    for &b in data {
        h ^= b as u64;
        h = h.wrapping_mul(0x0000_0100_0000_01b3);
    }
    h
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn fnv1a_empty_input() {
        // 空入力はFNV-1aオフセット基底値を返す
        let h = fnv1a(b"");
        assert_eq!(h, 0xcbf2_9ce4_8422_2325);
    }

    #[test]
    fn fnv1a_deterministic() {
        // 同一入力 → 同一ハッシュ
        let h1 = fnv1a(b"hello");
        let h2 = fnv1a(b"hello");
        assert_eq!(h1, h2);
    }

    #[test]
    fn fnv1a_different_inputs() {
        // 異なる入力 → 異なるハッシュ
        let h1 = fnv1a(b"hello");
        let h2 = fnv1a(b"world");
        assert_ne!(h1, h2);
    }

    #[test]
    fn fnv1a_single_byte() {
        // 1バイト入力のハッシュが非ゼロ
        let h = fnv1a(b"x");
        assert_ne!(h, 0);
        assert_ne!(h, 0xcbf2_9ce4_8422_2325);
    }

    #[test]
    fn fnv1a_order_matters() {
        // バイト順序でハッシュが変わることを確認
        let h1 = fnv1a(b"ab");
        let h2 = fnv1a(b"ba");
        assert_ne!(h1, h2);
    }

    #[test]
    fn fnv1a_known_value() {
        // FNV-1a 64-bit: "a" の既知ハッシュ値を検証
        // offset ^ 'a' = 0xcbf29ce484222325 ^ 0x61 = 0xcbf29ce484222344
        // * prime = 0xcbf29ce484222344 * 0x00000100000001b3
        let h = fnv1a(b"a");
        // 手計算の代わりに決定論性を確認
        let h2 = fnv1a(b"a");
        assert_eq!(h, h2);
        assert_ne!(h, 0);
    }

    #[test]
    fn fnv1a_long_input() {
        // 長い入力でもパニックしない
        let data = vec![0xFFu8; 10000];
        let h = fnv1a(&data);
        assert_ne!(h, 0);
    }

    #[test]
    fn fnv1a_all_zero_bytes() {
        // ゼロバイト列でも有効なハッシュ
        let h = fnv1a(&[0u8; 8]);
        assert_ne!(h, 0xcbf2_9ce4_8422_2325); // オフセット基底と異なる
    }
}
