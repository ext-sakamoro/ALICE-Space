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
pub use link_budget::{friis_path_loss_db, LinkBudget, LinkBudgetResult};
pub use mission::{MissionEvent, MissionLog, MissionPhase};
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
