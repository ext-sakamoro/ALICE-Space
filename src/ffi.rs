//! C FFI bindings for ALICE-Space
//!
//! 12 extern "C" functions for Unity / UE5 / any C-compatible host.
//!
//! The functions that take buffers pass each one with its length (in `f64`
//! elements) and return a status: [`ALICE_SPACE_OK`], or a nonzero error. They
//! check every pointer and length, read every input and compute every output
//! before writing; on a nonzero status nothing is written. A non-finite input
//! (NaN or ±∞) or a non-finite result is an error, not a NaN in the output.
//! The functions that return one `f64` take no buffer and return NaN for an
//! input outside their domain, as the Rust functions do.
//!
//! Author: Moroya Sakamoto

use crate::comm::CommLink;
use crate::constellation::WalkerConstellation;
use crate::link_budget::LinkBudget;
use crate::orbit::SpacecraftState;
use crate::propagator::TwoBodyAccel;

/// The call succeeded and the output buffer holds the result.
pub const ALICE_SPACE_OK: i32 = 0;
/// A buffer pointer is null; nothing was written.
pub const ALICE_SPACE_ERR_NULL: i32 = 1;
/// A buffer has another length than the function reads or writes; nothing was written.
pub const ALICE_SPACE_ERR_LENGTH: i32 = 2;
/// An input is NaN or infinite; nothing was written.
pub const ALICE_SPACE_ERR_NON_FINITE_INPUT: i32 = 3;
/// The inputs are finite but the result is not (for example a zero radius); nothing was written.
pub const ALICE_SPACE_ERR_NON_FINITE_RESULT: i32 = 4;

/// `N` values read from `ptr`, which the caller says holds `len`.
///
/// # Safety
/// When `ptr` is not null and `len == N`, `ptr` must point to `N` readable `f64`.
unsafe fn read<const N: usize>(ptr: *const f64, len: usize) -> Result<[f64; N], i32> {
    if ptr.is_null() {
        return Err(ALICE_SPACE_ERR_NULL);
    }
    if len != N {
        return Err(ALICE_SPACE_ERR_LENGTH);
    }
    let mut v = [0.0; N];
    for (i, x) in v.iter_mut().enumerate() {
        *x = *ptr.add(i);
    }
    if v.iter().all(|x| x.is_finite()) {
        Ok(v)
    } else {
        Err(ALICE_SPACE_ERR_NON_FINITE_INPUT)
    }
}

/// Whether `ptr` can take `N` values, for a buffer the caller says holds `len`.
fn writable<const N: usize>(ptr: *mut f64, len: usize) -> Result<(), i32> {
    if ptr.is_null() {
        return Err(ALICE_SPACE_ERR_NULL);
    }
    if len != N {
        return Err(ALICE_SPACE_ERR_LENGTH);
    }
    Ok(())
}

fn finite(scalars: &[f64]) -> Result<(), i32> {
    if scalars.iter().all(|x| x.is_finite()) {
        Ok(())
    } else {
        Err(ALICE_SPACE_ERR_NON_FINITE_INPUT)
    }
}

/// Writes `v` to `ptr` when every value is finite, else nothing.
///
/// # Safety
/// `ptr` must have passed [`writable`] for `N`.
unsafe fn write<const N: usize>(ptr: *mut f64, v: [f64; N]) -> i32 {
    if !v.iter().all(|x| x.is_finite()) {
        return ALICE_SPACE_ERR_NON_FINITE_RESULT;
    }
    for (i, x) in v.into_iter().enumerate() {
        *ptr.add(i) = x;
    }
    ALICE_SPACE_OK
}

fn status(r: Result<i32, i32>) -> i32 {
    r.unwrap_or_else(|e| e)
}

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

/// Hohmann transfer delta-v. Writes dv1 and dv2 (km/s) to `out` (`out_len` = 2)
/// and returns a status ([`ALICE_SPACE_OK`] or an error, with nothing written).
///
/// # Safety
/// When `out` is not null and `out_len` is 2, `out` must point to 2 writable f64.
#[no_mangle]
pub unsafe extern "C" fn alice_space_delta_v_hohmann(
    r1_km: f64,
    r2_km: f64,
    mu: f64,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    status((|| {
        writable::<2>(out, out_len)?;
        finite(&[r1_km, r2_km, mu])?;
        let (dv1, dv2) = crate::orbit::delta_v_hohmann(r1_km, r2_km, mu);
        Ok(write(out, [dv1, dv2]))
    })())
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

/// Single RK4 propagation step. Reads `[px, py, pz, vx, vy, vz]` from
/// `state_in` (`state_in_len` = 6), writes the new state in the same order to
/// `state_out` (`state_out_len` = 6) and returns a status ([`ALICE_SPACE_OK`] or
/// an error, with nothing written). `state_in` and `state_out` may be the same
/// buffer: every input is read before the first write.
///
/// # Safety
/// When the pointers are not null and the lengths are 6, `state_in` must point
/// to 6 readable f64 and `state_out` to 6 writable f64.
#[no_mangle]
pub unsafe extern "C" fn alice_space_rk4_single(
    state_in: *const f64,
    state_in_len: usize,
    mu: f64,
    dt_s: f64,
    state_out: *mut f64,
    state_out_len: usize,
) -> i32 {
    status((|| {
        writable::<6>(state_out, state_out_len)?;
        let s = read::<6>(state_in, state_in_len)?;
        finite(&[mu, dt_s])?;
        let sc = SpacecraftState {
            position_km: [s[0], s[1], s[2]],
            velocity_km_s: [s[3], s[4], s[5]],
            timestamp_ns: 0,
            fuel_kg: 0.0,
        };
        let r = crate::propagator::propagate_rk4_single(&sc, mu, dt_s);
        let (p, v) = (r.position_km, r.velocity_km_s);
        Ok(write(state_out, [p[0], p[1], p[2], v[0], v[1], v[2]]))
    })())
}

/// Two-body gravitational acceleration. Reads `[x, y, z]` from `r` (`r_len` =
/// 3), writes `[ax, ay, az]` to `out` (`out_len` = 3) and returns a status
/// ([`ALICE_SPACE_OK`] or an error, with nothing written).
///
/// # Safety
/// When the pointers are not null and the lengths are 3, `r` must point to 3
/// readable f64 and `out` to 3 writable f64.
#[no_mangle]
pub unsafe extern "C" fn alice_space_two_body_accel(
    r: *const f64,
    r_len: usize,
    mu: f64,
    out: *mut f64,
    out_len: usize,
) -> i32 {
    status((|| {
        writable::<3>(out, out_len)?;
        let pos = read::<3>(r, r_len)?;
        finite(&[mu])?;
        Ok(write(out, TwoBodyAccel::new(mu).acceleration(&pos)))
    })())
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
    use alice_det_math::sqrt64;

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
        let s = unsafe { alice_space_delta_v_hohmann(6578.0, 42164.0, MU_EARTH, out.as_mut_ptr(), 2) };
        assert_eq!(s, ALICE_SPACE_OK);
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
        let s = unsafe {
            alice_space_rk4_single(state_in.as_ptr(), 6, MU_EARTH, 60.0, state_out.as_mut_ptr(), 6)
        };
        assert_eq!(s, ALICE_SPACE_OK);
        let r = sqrt64(
            (state_out[0] * state_out[0])
                + state_out[1] * state_out[1]
                + state_out[2] * state_out[2],
        );
        assert!((r - 6778.0).abs() < 50.0);
    }

    #[test]
    fn test_two_body_accel_ffi() {
        let r: [f64; 3] = [6778.0, 0.0, 0.0];
        let mut out = [0.0f64; 3];
        let s = unsafe { alice_space_two_body_accel(r.as_ptr(), 3, MU_EARTH, out.as_mut_ptr(), 3) };
        assert_eq!(s, ALICE_SPACE_OK);
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
