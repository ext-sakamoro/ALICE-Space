//! Atmospheric drag on a satellite, with the density supplied by the caller.
//!
//! a = −½ ρ Cd (A/m) |v_r| v_r, with v_r = v − ω × r the velocity relative to
//! an atmosphere that rotates with the body about z at rate ω. Positions and
//! velocities are in km and km/s, the density in kg/m³, the area-to-mass ratio
//! in m²/kg, the result in km/s². The density model (exponential, Harris-
//! Priester, MSIS, ...) is the caller's: this function is the drag law alone.
//!
//! Author: Moroya Sakamoto

use alice_det_math::sqrt64;

/// The Earth's rotation rate of GRS80 / WGS84, rad/s.
pub const EARTH_ROTATION_RATE: f64 = 7.292_115e-5;

/// Why no acceleration was computed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum DragError {
    /// A non-finite state, or a negative or non-finite density, drag
    /// coefficient or area-to-mass ratio.
    InvalidInput,
}

impl core::fmt::Display for DragError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str("invalid input to the drag")
    }
}

impl std::error::Error for DragError {}

/// The drag acceleration (km/s²) at position `r` (km) and velocity `v`
/// (km/s) in density `rho` (kg/m³), with drag coefficient `cd`, area-to-mass
/// ratio `area_to_mass` (m²/kg) and atmospheric rotation rate `omega` (rad/s,
/// about z; zero for a non-rotating atmosphere).
pub fn atmospheric_drag(
    r: [f64; 3],
    v: [f64; 3],
    rho: f64,
    cd: f64,
    area_to_mass: f64,
    omega: f64,
) -> Result<[f64; 3], DragError> {
    // a non-negative finite coefficient: NaN fails `>= 0.0`, infinity fails
    // `is_finite`
    let usable = |x: f64| x >= 0.0 && x.is_finite();
    let finite = r.iter().chain(&v).all(|x| x.is_finite()) && omega.is_finite();
    if !(finite && usable(rho) && usable(cd) && usable(area_to_mass)) {
        return Err(DragError::InvalidInput);
    }
    let vr = [v[0] + omega * r[1], v[1] - omega * r[0], v[2]];
    let speed = sqrt64(vr[2] * vr[2] + (vr[0] * vr[0] + vr[1] * vr[1]));
    // ½ ρ Cd (A/m) |v_r| v_r in m/s² with v in m/s is 10⁶ times the km/s
    // product, and 10³ m/s² per km/s²
    let k = -0.5 * rho * cd * area_to_mass * speed * 1_000.0;
    Ok([k * vr[0], k * vr[1], k * vr[2]])
}
