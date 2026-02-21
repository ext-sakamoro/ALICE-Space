//! ALICE-Space — Deep-Space Communication & Autonomous Control
//!
//! Handles deep-space communication with extreme latency (minutes to hours)
//! and minimal bandwidth (kbps), enabling autonomous spacecraft control via
//! mathematical model differentials rather than raw data transmission.
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

pub mod orbit;
pub mod comm;
pub mod autonomy;
pub mod mission;

pub use orbit::{BodyId, OrbitalElements, CelestialBody, SpacecraftState, orbital_period, orbital_velocity, delta_v_hohmann, light_delay_s};
pub use comm::{CommLink, ModelDifferential, can_transmit};
pub use autonomy::{AutonomyLevel, TrajectoryModel, ControlDecision, apply_differential, compute_correction};
pub use mission::{MissionPhase, MissionEvent, MissionLog};
