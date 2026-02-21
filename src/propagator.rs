//! RK4 numerical orbit propagator
//!
//! 4th-order Runge-Kutta integration for two-body orbital mechanics.
//! Propagates spacecraft state vectors forward in time using the
//! gravitational acceleration a = -μ r / |r|³.
//!
//! Author: Moroya Sakamoto

use crate::orbit::SpacecraftState;

/// Two-body gravitational acceleration calculator.
///
/// Computes a = -μ r / |r|³ where μ is the gravitational parameter
/// and r is the position vector from the central body.
#[derive(Debug, Clone, Copy)]
pub struct TwoBodyAccel {
    /// Gravitational parameter μ = G*M (km³/s²).
    pub mu: f64,
}

impl TwoBodyAccel {
    /// Create a new two-body acceleration model.
    pub fn new(mu: f64) -> Self {
        Self { mu }
    }

    /// Compute gravitational acceleration at position r.
    /// Returns [ax, ay, az] in km/s².
    #[inline]
    pub fn acceleration(&self, r: &[f64; 3]) -> [f64; 3] {
        let r2 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];
        let r_mag = r2.sqrt();
        if r_mag < 1e-15 {
            return [0.0, 0.0, 0.0];
        }
        // -mu / |r|^3, using reciprocal to avoid 3 divisions
        let factor = -self.mu / (r_mag * r2);
        [r[0] * factor, r[1] * factor, r[2] * factor]
    }
}

/// 6-element state: [x, y, z, vx, vy, vz]
type State6 = [f64; 6];

/// Pack position and velocity into a 6-element state.
#[inline]
fn pack(pos: &[f64; 3], vel: &[f64; 3]) -> State6 {
    [pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]]
}

/// Compute derivative: d/dt [r, v] = [v, a(r)]
#[inline]
fn deriv(state: &State6, accel: &TwoBodyAccel) -> State6 {
    let r = [state[0], state[1], state[2]];
    let a = accel.acceleration(&r);
    [state[3], state[4], state[5], a[0], a[1], a[2]]
}

/// Scale a 6-element array by a scalar.
#[inline]
fn scale(s: &State6, h: f64) -> State6 {
    [s[0] * h, s[1] * h, s[2] * h, s[3] * h, s[4] * h, s[5] * h]
}

/// Add two 6-element arrays.
#[inline]
fn add(a: &State6, b: &State6) -> State6 {
    [
        a[0] + b[0], a[1] + b[1], a[2] + b[2],
        a[3] + b[3], a[4] + b[4], a[5] + b[5],
    ]
}

/// Perform a single RK4 step.
///
/// Returns the new state after advancing by `dt` seconds.
fn rk4_step(state: &State6, dt: f64, accel: &TwoBodyAccel) -> State6 {
    let k1 = deriv(state, accel);
    let k2 = deriv(&add(state, &scale(&k1, dt * 0.5)), accel);
    let k3 = deriv(&add(state, &scale(&k2, dt * 0.5)), accel);
    let k4 = deriv(&add(state, &scale(&k3, dt)), accel);

    // y_{n+1} = y_n + (dt/6)(k1 + 2*k2 + 2*k3 + k4)
    let h6 = dt / 6.0;
    [
        state[0] + h6 * (k1[0] + 2.0 * k2[0] + 2.0 * k3[0] + k4[0]),
        state[1] + h6 * (k1[1] + 2.0 * k2[1] + 2.0 * k3[1] + k4[1]),
        state[2] + h6 * (k1[2] + 2.0 * k2[2] + 2.0 * k3[2] + k4[2]),
        state[3] + h6 * (k1[3] + 2.0 * k2[3] + 2.0 * k3[3] + k4[3]),
        state[4] + h6 * (k1[4] + 2.0 * k2[4] + 2.0 * k3[4] + k4[4]),
        state[5] + h6 * (k1[5] + 2.0 * k2[5] + 2.0 * k3[5] + k4[5]),
    ]
}

/// Propagate a single RK4 step from a `SpacecraftState`.
///
/// Returns the state after `dt_s` seconds. Fuel is unchanged.
pub fn propagate_rk4_single(
    state: &SpacecraftState,
    mu: f64,
    dt_s: f64,
) -> SpacecraftState {
    let accel = TwoBodyAccel::new(mu);
    let s = pack(&state.position_km, &state.velocity_km_s);
    let s_new = rk4_step(&s, dt_s, &accel);
    let dt_ns = (dt_s * 1e9) as u64;

    SpacecraftState {
        position_km: [s_new[0], s_new[1], s_new[2]],
        velocity_km_s: [s_new[3], s_new[4], s_new[5]],
        timestamp_ns: state.timestamp_ns + dt_ns,
        fuel_kg: state.fuel_kg,
    }
}

/// Propagate a spacecraft state forward using RK4 integration.
///
/// Advances the state by `steps` steps of `dt_s` seconds each.
/// Returns a Vec of states at each step (length = steps + 1, including initial).
pub fn propagate_rk4(
    initial: &SpacecraftState,
    mu: f64,
    dt_s: f64,
    steps: usize,
) -> Vec<SpacecraftState> {
    let accel = TwoBodyAccel::new(mu);
    let mut results = Vec::with_capacity(steps + 1);
    let mut s = pack(&initial.position_km, &initial.velocity_km_s);
    let mut t_ns = initial.timestamp_ns;
    let dt_ns = (dt_s * 1e9) as u64;

    results.push(initial.clone());

    for _ in 0..steps {
        s = rk4_step(&s, dt_s, &accel);
        t_ns += dt_ns;
        results.push(SpacecraftState {
            position_km: [s[0], s[1], s[2]],
            velocity_km_s: [s[3], s[4], s[5]],
            timestamp_ns: t_ns,
            fuel_kg: initial.fuel_kg,
        });
    }

    results
}

/// Compute specific orbital energy: E = v²/2 - μ/r
#[inline]
pub fn specific_energy(pos: &[f64; 3], vel: &[f64; 3], mu: f64) -> f64 {
    let r = (pos[0] * pos[0] + pos[1] * pos[1] + pos[2] * pos[2]).sqrt();
    let v2 = vel[0] * vel[0] + vel[1] * vel[1] + vel[2] * vel[2];
    v2 * 0.5 - mu / r
}

// ── Tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const MU_EARTH: f64 = 398600.4418;

    fn circular_leo_state() -> SpacecraftState {
        // ISS-like orbit: r = 6778 km, v = sqrt(mu/r) ≈ 7.67 km/s
        let r = 6778.0;
        let v = (MU_EARTH / r).sqrt();
        SpacecraftState {
            position_km: [r, 0.0, 0.0],
            velocity_km_s: [0.0, v, 0.0],
            timestamp_ns: 0,
            fuel_kg: 100.0,
        }
    }

    #[test]
    fn acceleration_at_earth_surface() {
        let accel = TwoBodyAccel::new(MU_EARTH);
        let r = [6371.0, 0.0, 0.0];
        let a = accel.acceleration(&r);
        // g ≈ 9.82e-3 km/s² at Earth surface
        assert!((a[0] + 9.82e-3).abs() < 0.1e-3);
        assert!(a[1].abs() < 1e-15);
        assert!(a[2].abs() < 1e-15);
    }

    #[test]
    fn acceleration_zero_position() {
        let accel = TwoBodyAccel::new(MU_EARTH);
        let a = accel.acceleration(&[0.0, 0.0, 0.0]);
        assert_eq!(a, [0.0, 0.0, 0.0]);
    }

    #[test]
    fn single_step_preserves_energy() {
        let state = circular_leo_state();
        let e_before = specific_energy(&state.position_km, &state.velocity_km_s, MU_EARTH);
        let new_state = propagate_rk4_single(&state, MU_EARTH, 10.0);
        let e_after = specific_energy(&new_state.position_km, &new_state.velocity_km_s, MU_EARTH);
        // RK4 should conserve energy well over short steps
        assert!((e_before - e_after).abs() / e_before.abs() < 1e-10);
    }

    #[test]
    fn propagate_one_orbit_returns_to_start() {
        let state = circular_leo_state();
        let r = 6778.0;
        let period = 2.0 * std::f64::consts::PI * (r * r * r / MU_EARTH).sqrt();
        let dt = 1.0; // 1-second steps for higher accuracy
        let steps = (period / dt).round() as usize;
        let trajectory = propagate_rk4(&state, MU_EARTH, dt, steps);

        // After one full orbit, should return near initial position
        let final_state = trajectory.last().unwrap();
        let dx = final_state.position_km[0] - state.position_km[0];
        let dy = final_state.position_km[1] - state.position_km[1];
        let dz = final_state.position_km[2] - state.position_km[2];
        let error = (dx * dx + dy * dy + dz * dz).sqrt();
        // Error < 5 km after full orbit with 1s steps (RK4 truncation + rounding)
        assert!(error < 5.0, "Position error after one orbit: {} km", error);
    }

    #[test]
    fn propagate_energy_conservation() {
        let state = circular_leo_state();
        let e_initial = specific_energy(&state.position_km, &state.velocity_km_s, MU_EARTH);
        let trajectory = propagate_rk4(&state, MU_EARTH, 10.0, 100);

        // Check energy at every step
        for s in &trajectory {
            let e = specific_energy(&s.position_km, &s.velocity_km_s, MU_EARTH);
            let rel_error = (e - e_initial).abs() / e_initial.abs();
            assert!(rel_error < 1e-7, "Energy drift: {}", rel_error);
        }
    }

    #[test]
    fn propagate_radius_stays_constant_circular() {
        let state = circular_leo_state();
        let r_initial = 6778.0;
        let trajectory = propagate_rk4(&state, MU_EARTH, 30.0, 50);

        for s in &trajectory {
            let r = (s.position_km[0].powi(2) + s.position_km[1].powi(2) + s.position_km[2].powi(2)).sqrt();
            assert!((r - r_initial).abs() < 0.1, "Radius drift: {} km", (r - r_initial).abs());
        }
    }

    #[test]
    fn propagate_timestamps_advance() {
        let state = circular_leo_state();
        let trajectory = propagate_rk4(&state, MU_EARTH, 60.0, 10);
        assert_eq!(trajectory.len(), 11);
        assert_eq!(trajectory[0].timestamp_ns, 0);
        assert_eq!(trajectory[1].timestamp_ns, 60_000_000_000);
        assert_eq!(trajectory[10].timestamp_ns, 600_000_000_000);
    }

    #[test]
    fn propagate_fuel_unchanged() {
        let state = circular_leo_state();
        let trajectory = propagate_rk4(&state, MU_EARTH, 60.0, 5);
        for s in &trajectory {
            assert!((s.fuel_kg - 100.0).abs() < 1e-15);
        }
    }

    #[test]
    fn specific_energy_circular() {
        // Circular orbit: E = -μ/(2a)
        let r = 6778.0;
        let v = (MU_EARTH / r).sqrt();
        let e = specific_energy(&[r, 0.0, 0.0], &[0.0, v, 0.0], MU_EARTH);
        let e_expected = -MU_EARTH / (2.0 * r);
        assert!((e - e_expected).abs() < 1e-6);
    }

    #[test]
    fn propagate_empty_steps() {
        let state = circular_leo_state();
        let trajectory = propagate_rk4(&state, MU_EARTH, 60.0, 0);
        assert_eq!(trajectory.len(), 1);
    }

    #[test]
    fn propagate_single_matches_batch() {
        let state = circular_leo_state();
        let single = propagate_rk4_single(&state, MU_EARTH, 60.0);
        let batch = propagate_rk4(&state, MU_EARTH, 60.0, 1);
        let b = &batch[1];
        for (s_p, b_p) in single.position_km.iter().zip(b.position_km.iter()) {
            assert!((s_p - b_p).abs() < 1e-12);
        }
        for (s_v, b_v) in single.velocity_km_s.iter().zip(b.velocity_km_s.iter()) {
            assert!((s_v - b_v).abs() < 1e-12);
        }
    }

    #[test]
    fn acceleration_direction_opposes_position() {
        let accel = TwoBodyAccel::new(MU_EARTH);
        // Position along +Y axis
        let r = [0.0, 7000.0, 0.0];
        let a = accel.acceleration(&r);
        // Acceleration should point toward origin → negative Y
        assert!(a[1] < 0.0, "Accel should be -Y, got {}", a[1]);
        assert!(a[0].abs() < 1e-15);
        assert!(a[2].abs() < 1e-15);
    }

    #[test]
    fn acceleration_scales_with_mu() {
        let accel1 = TwoBodyAccel::new(MU_EARTH);
        let accel2 = TwoBodyAccel::new(MU_EARTH * 2.0);
        let r = [7000.0, 0.0, 0.0];
        let a1 = accel1.acceleration(&r);
        let a2 = accel2.acceleration(&r);
        // Double mu → double acceleration
        assert!((a2[0] / a1[0] - 2.0).abs() < 1e-10);
    }

    #[test]
    fn specific_energy_negative_for_bound_orbit() {
        let r = 6778.0;
        let v = (MU_EARTH / r).sqrt(); // circular orbit
        let e = specific_energy(&[r, 0.0, 0.0], &[0.0, v, 0.0], MU_EARTH);
        assert!(e < 0.0, "Bound orbit should have negative energy: {}", e);
    }

    #[test]
    fn specific_energy_positive_for_escape() {
        let r = 6778.0;
        let v_escape = (2.0 * MU_EARTH / r).sqrt();
        let v = v_escape * 1.1; // above escape velocity
        let e = specific_energy(&[r, 0.0, 0.0], &[0.0, v, 0.0], MU_EARTH);
        assert!(e > 0.0, "Escape trajectory should have positive energy: {}", e);
    }
}
