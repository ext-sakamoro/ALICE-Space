//! Solar radiation pressure on a satellite: the cannonball model scaled by the
//! visible fraction of the solar disc.
//!
//! a = ν · (S/c) · Cr · (A/m) · (AU/d)² · d̂, with d the vector from the Sun to
//! the satellite, S the nominal total solar irradiance at 1 au (1361 W/m²,
//! IAU 2015 Resolution B3), c the speed of light, Cr the reflectivity
//! coefficient (1 absorbs, 2 reflects back), A/m the area-to-mass ratio, and ν
//! the visible fraction of the Sun's disc from [`crate::eclipse`]'s conical
//! shadow of each occulting body. With several occulting bodies ν is the
//! product of their fractions, exact when at most one of them covers part of
//! the Sun at a time (two discs overlapping the Sun at once are not combined
//! geometrically).
//!
//! Author: Moroya Sakamoto

use alice_det_math::sqrt64;

use crate::eclipse::{conical_illumination, EclipseError};

/// Nominal total solar irradiance at 1 au (IAU 2015 Resolution B3), W/m².
pub const SOLAR_IRRADIANCE_W_M2: f64 = 1361.0;
/// Speed of light, m/s (exact).
pub const SPEED_OF_LIGHT_M_S: f64 = 299_792_458.0;
/// The astronomical unit (IAU 2012), km.
pub const AU_KM: f64 = 149_597_870.7;

/// A spherical body that can cast a shadow.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Occulter {
    /// Position (km), in the same frame as the satellite and the Sun.
    pub position: [f64; 3],
    /// Radius (km).
    pub radius: f64,
}

/// Why no acceleration was computed.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum SrpError {
    /// A non-finite coordinate, a negative or non-finite coefficient or
    /// area-to-mass ratio, or the satellite at the Sun's centre.
    InvalidInput,
    /// The shadow of an occulting body could not be computed.
    Eclipse(EclipseError),
}

impl core::fmt::Display for SrpError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::InvalidInput => f.write_str("invalid input to the radiation pressure"),
            Self::Eclipse(e) => write!(f, "shadow: {e}"),
        }
    }
}

impl std::error::Error for SrpError {}

/// The acceleration (km/s²) of solar radiation pressure on a satellite at
/// `sat` with the Sun at `sun` (km, same frame), shadowed by `occulters`, with
/// reflectivity coefficient `cr` and area-to-mass ratio `area_to_mass`
/// (m²/kg).
pub fn solar_radiation_pressure(
    sat: [f64; 3],
    sun: [f64; 3],
    occulters: &[Occulter],
    cr: f64,
    area_to_mass: f64,
) -> Result<[f64; 3], SrpError> {
    // a non-negative finite coefficient: NaN fails `>= 0.0`, infinity fails
    // `is_finite`
    let usable = |x: f64| x >= 0.0 && x.is_finite();
    let finite = sat.iter().chain(&sun).all(|v| v.is_finite());
    if !(finite && usable(cr) && usable(area_to_mass)) {
        return Err(SrpError::InvalidInput);
    }
    let d = [sat[0] - sun[0], sat[1] - sun[1], sat[2] - sun[2]];
    let dn = sqrt64(d[2] * d[2] + (d[0] * d[0] + d[1] * d[1]));
    if dn == 0.0 {
        return Err(SrpError::InvalidInput);
    }
    let mut visible = 1.0;
    for o in occulters {
        visible *= conical_illumination(sat, sun, o.position, o.radius)
            .map_err(SrpError::Eclipse)?
            .fraction;
    }
    let scale = AU_KM / dn;
    // (S/c) in N/m² = kg/(m·s²), times m²/kg gives m/s², /1000 for km/s²
    let a = visible
        * (SOLAR_IRRADIANCE_W_M2 / SPEED_OF_LIGHT_M_S)
        * cr
        * area_to_mass
        * (scale * scale)
        / 1000.0;
    let k = a / dn;
    Ok([k * d[0], k * d[1], k * d[2]])
}
