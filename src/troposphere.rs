//! Tropospheric delay correction — Saastamoinen model.
//!
//! The Saastamoinen (1972) model splits the tropospheric delay into a dry
//! (hydrostatic) and a wet component, each mapped from the zenith direction
//! to the receiver-to-satellite line-of-sight by a mapping function.
//!
//! For SPACID's near-Earth GNSS use cases the standard-atmosphere / Niell
//! wet-dry decomposition is sufficient. The interface accepts either
//! measured meteorology or an ISO 2533 default derived from the user
//! altitude alone.

use crate::geodetic::Geodetic;

// ---------------------------------------------------------------------------
// Meteorological state
// ---------------------------------------------------------------------------

/// Surface meteorology at the receiver.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct MetState {
    /// Pressure in hPa (mbar).
    pub pressure_hpa: f64,
    /// Temperature in Kelvin.
    pub temperature_k: f64,
    /// Partial pressure of water vapour in hPa.
    pub water_vapour_hpa: f64,
}

impl MetState {
    /// Constructor.
    #[must_use]
    pub const fn new(pressure_hpa: f64, temperature_k: f64, water_vapour_hpa: f64) -> Self {
        Self {
            pressure_hpa,
            temperature_k,
            water_vapour_hpa,
        }
    }

    /// Standard atmosphere at sea level: 1013.25 hPa, 288.15 K, 50% RH.
    #[must_use]
    pub fn standard_at_altitude(altitude_m: f64) -> Self {
        // ISO 2533 dry-adiabatic lapse below the tropopause (11 km).
        let sea_level_pressure = 1013.25;
        let sea_level_temperature = 288.15;
        let lapse_rate = 6.5e-3; // K / m
        let alt_clamped = altitude_m.clamp(-500.0, 11_000.0);
        let t = sea_level_temperature - lapse_rate * alt_clamped;
        let p = sea_level_pressure * (t / sea_level_temperature).powf(5.2559);
        // Assume 50% relative humidity, using saturation vapour pressure from
        // Magnus formula.
        let t_c = t - 273.15;
        let es = 6.112 * ((17.62 * t_c) / (243.12 + t_c)).exp();
        let e = 0.5 * es;
        Self::new(p, t, e)
    }
}

// ---------------------------------------------------------------------------
// Saastamoinen zenith delays
// ---------------------------------------------------------------------------

/// Zenith hydrostatic (dry) delay in metres.
///
/// The Saastamoinen expression scales with surface pressure and a
/// latitude/altitude term encoding the local gravity anomaly.
#[must_use]
pub fn zenith_hydrostatic_delay_m(met: MetState, lat_rad: f64, altitude_m: f64) -> f64 {
    let cos_2phi = (2.0 * lat_rad).cos();
    let f = 1.0 - 0.002_66 * cos_2phi - 0.000_28 * altitude_m / 1000.0;
    // 0.0022768 m/hPa is the Saastamoinen constant for dry delay.
    0.002_276_8 * met.pressure_hpa / f
}

/// Zenith wet delay in metres.
///
/// Saastamoinen's original wet model uses partial pressure of water vapour
/// and receiver temperature.
#[must_use]
pub fn zenith_wet_delay_m(met: MetState) -> f64 {
    let t = met.temperature_k;
    0.002_277 * (1255.0 / t + 0.05) * met.water_vapour_hpa
}

// ---------------------------------------------------------------------------
// Mapping functions
// ---------------------------------------------------------------------------

/// Niell / secant mapping function `m(elev) = 1/sin(elev)`.
///
/// Suitable above about 5° elevation. For SPACID's cm-scale use case a
/// full Niell decomposition can replace this; the simple form is exposed
/// here as a public API for downstream refinement.
#[must_use]
pub fn simple_mapping_function(elevation_rad: f64) -> f64 {
    let sin_e = elevation_rad.sin();
    if sin_e < 1e-6 {
        // Guard against horizon rays.
        1e6
    } else {
        1.0 / sin_e
    }
}

// ---------------------------------------------------------------------------
// Full slant delay
// ---------------------------------------------------------------------------

/// Total tropospheric delay along the line-of-sight in metres.
///
/// Combines the Saastamoinen zenith hydrostatic + wet delays with the
/// [`simple_mapping_function`].
#[must_use]
pub fn saastamoinen_delay_m(met: MetState, user: Geodetic, elevation_rad: f64) -> f64 {
    let zhd = zenith_hydrostatic_delay_m(met, user.lat_rad, user.alt_m);
    let zwd = zenith_wet_delay_m(met);
    let m = simple_mapping_function(elevation_rad);
    (zhd + zwd) * m
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use core::f64::consts::TAU;

    fn tokyo() -> Geodetic {
        Geodetic::from_degrees(35.6895, 139.6917, 40.0)
    }

    #[test]
    fn standard_atmosphere_at_sea_level_matches_iso_2533() {
        let met = MetState::standard_at_altitude(0.0);
        assert!((met.pressure_hpa - 1013.25).abs() < 1e-6);
        assert!((met.temperature_k - 288.15).abs() < 1e-6);
    }

    #[test]
    fn zenith_hydrostatic_delay_is_around_2_3_m_at_sea_level() {
        let met = MetState::standard_at_altitude(0.0);
        let zhd = zenith_hydrostatic_delay_m(met, tokyo().lat_rad, 0.0);
        // Textbook value is 2.30 m ± 0.02 m.
        assert!((2.25..=2.35).contains(&zhd), "zhd = {zhd}");
    }

    #[test]
    fn zenith_wet_delay_is_smaller_than_dry() {
        let met = MetState::standard_at_altitude(0.0);
        let zwd = zenith_wet_delay_m(met);
        let zhd = zenith_hydrostatic_delay_m(met, tokyo().lat_rad, 0.0);
        assert!(zwd < zhd);
        // Wet delay is typically 5-40 cm at sea level.
        assert!(zwd > 0.0 && zwd < 0.5, "zwd = {zwd}");
    }

    #[test]
    fn simple_mapping_grows_toward_horizon() {
        assert!(simple_mapping_function(TAU / 4.0) < simple_mapping_function(TAU / 12.0));
    }

    #[test]
    fn saastamoinen_delay_increases_at_low_elevation() {
        let met = MetState::standard_at_altitude(40.0);
        let high = saastamoinen_delay_m(met, tokyo(), 80.0_f64.to_radians());
        let low = saastamoinen_delay_m(met, tokyo(), 10.0_f64.to_radians());
        assert!(low > high);
    }

    #[test]
    fn high_altitude_reduces_zhd() {
        let sea = MetState::standard_at_altitude(0.0);
        let mtn = MetState::standard_at_altitude(3000.0);
        let sea_zhd = zenith_hydrostatic_delay_m(sea, tokyo().lat_rad, 0.0);
        let mtn_zhd = zenith_hydrostatic_delay_m(mtn, tokyo().lat_rad, 3000.0);
        assert!(mtn_zhd < sea_zhd);
    }
}
