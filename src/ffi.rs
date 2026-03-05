//! C FFI bindings for ALICE-Space
//!
//! 12 extern "C" functions for Unity / UE5 / any C-compatible host.
//!
//! Author: Moroya Sakamoto

use crate::comm::CommLink;
use crate::constellation::WalkerConstellation;
use crate::link_budget::LinkBudget;
use crate::orbit::SpacecraftState;
use crate::propagator::TwoBodyAccel;

// ---------------------------------------------------------------------------
// Orbit (4)
// ---------------------------------------------------------------------------

/// Orbital period T = 2π√(a³/μ) in seconds.
#[no_mangle]
pub extern "C" fn alice_space_orbital_period(semi_major_axis_km: f64, mu: f64) -> f64 {
    crate::orbit::orbital_period(semi_major_axis_km, mu)
}

/// Vis-viva orbital velocity in km/s.
#[no_mangle]
pub extern "C" fn alice_space_orbital_velocity(r_km: f64, a_km: f64, mu: f64) -> f64 {
    crate::orbit::orbital_velocity(r_km, a_km, mu)
}

/// One-way light delay in seconds for a given distance (km).
#[no_mangle]
pub extern "C" fn alice_space_light_delay_s(distance_km: f64) -> f64 {
    crate::orbit::light_delay_s(distance_km)
}

/// Hohmann transfer delta-v. Writes dv1 and dv2 (km/s) to `out`.
///
/// # Safety
/// `out` must point to a buffer of at least 2 f64.
#[no_mangle]
pub unsafe extern "C" fn alice_space_delta_v_hohmann(
    r1_km: f64,
    r2_km: f64,
    mu: f64,
    out: *mut f64,
) {
    let (dv1, dv2) = crate::orbit::delta_v_hohmann(r1_km, r2_km, mu);
    *out = dv1;
    *out.add(1) = dv2;
}

// ---------------------------------------------------------------------------
// LinkBudget (2)
// ---------------------------------------------------------------------------

/// Friis free-space path loss in dB.
#[no_mangle]
pub extern "C" fn alice_space_friis_path_loss_db(distance_km: f64, frequency_ghz: f64) -> f64 {
    crate::link_budget::friis_path_loss_db(distance_km, frequency_ghz)
}

/// Compute link budget with default parameters, returning `margin_db`.
/// Use this for quick go/no-go link closure checks.
#[no_mangle]
pub extern "C" fn alice_space_link_budget_margin(
    distance_km: f64,
    frequency_ghz: f64,
    data_rate_bps: f64,
) -> f64 {
    let lb = LinkBudget {
        frequency_ghz,
        distance_km,
        data_rate_bps,
        ..LinkBudget::default()
    };
    lb.compute().margin_db
}

// ---------------------------------------------------------------------------
// Propagator (2)
// ---------------------------------------------------------------------------

/// Single RK4 propagation step. Reads pos[3]+vel[3] from `state_in`,
/// writes new pos[3]+vel[3] to `state_out`.
///
/// # Safety
/// `state_in` must point to 6 contiguous f64 [px,py,pz,vx,vy,vz].
/// `state_out` must point to a buffer of at least 6 f64.
#[no_mangle]
pub unsafe extern "C" fn alice_space_rk4_single(
    state_in: *const f64,
    mu: f64,
    dt_s: f64,
    state_out: *mut f64,
) {
    let pos = [*state_in, *state_in.add(1), *state_in.add(2)];
    let vel = [*state_in.add(3), *state_in.add(4), *state_in.add(5)];
    let sc = SpacecraftState {
        position_km: pos,
        velocity_km_s: vel,
        timestamp_ns: 0,
        fuel_kg: 0.0,
    };
    let result = crate::propagator::propagate_rk4_single(&sc, mu, dt_s);
    *state_out = result.position_km[0];
    *state_out.add(1) = result.position_km[1];
    *state_out.add(2) = result.position_km[2];
    *state_out.add(3) = result.velocity_km_s[0];
    *state_out.add(4) = result.velocity_km_s[1];
    *state_out.add(5) = result.velocity_km_s[2];
}

/// Two-body gravitational acceleration. Writes [ax,ay,az] to `out`.
///
/// # Safety
/// `r` must point to 3 contiguous f64 [x,y,z].
/// `out` must point to a buffer of at least 3 f64.
#[no_mangle]
pub unsafe extern "C" fn alice_space_two_body_accel(r: *const f64, mu: f64, out: *mut f64) {
    let pos = [*r, *r.add(1), *r.add(2)];
    let accel = TwoBodyAccel::new(mu);
    let a = accel.acceleration(&pos);
    *out = a[0];
    *out.add(1) = a[1];
    *out.add(2) = a[2];
}

// ---------------------------------------------------------------------------
// Comm (2)
// ---------------------------------------------------------------------------

/// Create a `CommLink` and return its one-way latency in seconds.
#[no_mangle]
pub extern "C" fn alice_space_comm_latency_s(distance_km: f64) -> f64 {
    let link = CommLink::new(0, 1, distance_km, 1000.0);
    link.latency_s()
}

/// Bits available in a transmission window.
#[no_mangle]
pub extern "C" fn alice_space_comm_bits_per_window(bandwidth_bps: f64, window_s: f64) -> f64 {
    let link = CommLink::new(0, 1, 0.0, bandwidth_bps);
    link.bits_per_window(window_s)
}

// ---------------------------------------------------------------------------
// Constellation (2)
// ---------------------------------------------------------------------------

/// Walker constellation plane spacing in radians.
#[no_mangle]
pub extern "C" fn alice_space_walker_plane_spacing_rad(num_planes: u32) -> f64 {
    let w = WalkerConstellation::new(num_planes, num_planes, 0, 7000.0, 0.0);
    w.plane_spacing_rad()
}

/// Walker constellation ground track period in seconds.
#[no_mangle]
pub extern "C" fn alice_space_walker_ground_track_period_s(
    semi_major_axis_km: f64,
    mu: f64,
) -> f64 {
    let w = WalkerConstellation::new(1, 1, 0, semi_major_axis_km, 0.0);
    w.ground_track_period_s(mu)
}

#[cfg(test)]
mod tests {
    use super::*;

    const MU_EARTH: f64 = 398600.4418;

    #[test]
    fn test_orbital_period_ffi() {
        let t = alice_space_orbital_period(6778.0, MU_EARTH);
        assert!((t / 60.0 - 92.3).abs() < 1.0);
    }

    #[test]
    fn test_orbital_velocity_ffi() {
        let v = alice_space_orbital_velocity(6778.0, 6778.0, MU_EARTH);
        assert!(v > 7.0 && v < 8.0);
    }

    #[test]
    fn test_light_delay_ffi() {
        let d = alice_space_light_delay_s(384400.0);
        assert!((d - 1.28).abs() < 0.01);
    }

    #[test]
    fn test_delta_v_hohmann_ffi() {
        let mut out = [0.0f64; 2];
        unsafe { alice_space_delta_v_hohmann(6578.0, 42164.0, MU_EARTH, out.as_mut_ptr()) };
        assert!(out[0] > 0.0);
        assert!(out[1] > 0.0);
    }

    #[test]
    fn test_friis_path_loss_ffi() {
        let loss = alice_space_friis_path_loss_db(384400.0, 8.4);
        assert!(loss > 200.0);
    }

    #[test]
    fn test_link_budget_margin_ffi() {
        let margin = alice_space_link_budget_margin(384400.0, 8.4, 1000.0);
        assert!(margin.is_finite());
    }

    #[test]
    fn test_rk4_single_ffi() {
        let state_in: [f64; 6] = [6778.0, 0.0, 0.0, 0.0, 7.67, 0.0];
        let mut state_out = [0.0f64; 6];
        unsafe {
            alice_space_rk4_single(state_in.as_ptr(), MU_EARTH, 60.0, state_out.as_mut_ptr())
        };
        let r = (state_out[0].powi(2) + state_out[1].powi(2) + state_out[2].powi(2)).sqrt();
        assert!((r - 6778.0).abs() < 50.0);
    }

    #[test]
    fn test_two_body_accel_ffi() {
        let r: [f64; 3] = [6778.0, 0.0, 0.0];
        let mut out = [0.0f64; 3];
        unsafe { alice_space_two_body_accel(r.as_ptr(), MU_EARTH, out.as_mut_ptr()) };
        assert!(out[0] < 0.0); // acceleration toward center
    }

    #[test]
    fn test_comm_latency_ffi() {
        let lat = alice_space_comm_latency_s(384400.0);
        assert!((lat - 1.28).abs() < 0.01);
    }

    #[test]
    fn test_comm_bits_per_window_ffi() {
        let bits = alice_space_comm_bits_per_window(1000.0, 10.0);
        assert!((bits - 10000.0).abs() < 1e-6);
    }

    #[test]
    fn test_walker_plane_spacing_ffi() {
        let sp = alice_space_walker_plane_spacing_rad(6);
        assert!((sp - std::f64::consts::PI / 3.0).abs() < 1e-10);
    }

    #[test]
    fn test_walker_ground_track_ffi() {
        let t = alice_space_walker_ground_track_period_s(6778.0, MU_EARTH);
        assert!((t / 60.0 - 92.3).abs() < 1.0);
    }
}
