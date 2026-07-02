//! Hopfield tropospheric delay model.
//!
//! The Hopfield model is a two-layer (dry + wet) empirical tropospheric
//! delay estimator introduced by Helen S. Hopfield in 1969. It uses
//! surface meteorological measurements (pressure, temperature, humidity)
//! to compute the zenith delay and a mapping function to slant it.
//!
//! Compared to Saastamoinen (implemented in [`crate::troposphere`]),
//! Hopfield gives similar accuracy near sea level but retains its
//! validity at higher altitudes because it explicitly integrates through
//! the troposphere top layer. Combined use in a PPP engine hedges
//! against modelling error at the level of ~5 cm zenith / ~1 m at
//! 5° elevation.
//!
//! # References
//!
//! - Hopfield, H. S. (1969), "Two-quartic tropospheric refractivity
//!   profile for correcting satellite data", J. Geophys. Res., 74(18),
//!   4487-4499.
//! - Kaplan, E. D. & Hegarty, C. J. (2017), §7.2.4.2.
//! - Leick, A., Rapoport, L., & Tatarnikov, D. (2015), "GPS Satellite
//!   Surveying", 4th ed., Wiley.

/// Surface meteorological conditions used by the Hopfield model.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SurfaceConditions {
    /// Atmospheric pressure in hectopascals (hPa = millibar).
    pub pressure_hpa: f64,
    /// Temperature in kelvin.
    pub temperature_k: f64,
    /// Relative humidity in percent (0-100).
    pub relative_humidity_percent: f64,
    /// Receiver antenna height in metres above the reference ellipsoid.
    pub height_m: f64,
}

impl SurfaceConditions {
    /// Standard atmosphere at sea level (`ICAO`).
    #[must_use]
    pub const fn standard_atmosphere() -> Self {
        Self {
            pressure_hpa: 1013.25,
            temperature_k: 288.15,
            relative_humidity_percent: 50.0,
            height_m: 0.0,
        }
    }

    /// Saturation water vapour pressure in hPa (Magnus formula).
    #[must_use]
    pub fn saturation_pressure_hpa(&self) -> f64 {
        let t_c = self.temperature_k - 273.15;
        6.1078 * (17.27 * t_c / (t_c + 237.3)).exp()
    }

    /// Partial water vapour pressure in hPa.
    #[must_use]
    pub fn water_vapour_pressure_hpa(&self) -> f64 {
        self.saturation_pressure_hpa() * self.relative_humidity_percent / 100.0
    }
}

// ---------------------------------------------------------------------------
// Zenith delay
// ---------------------------------------------------------------------------

/// Hopfield zenith delay components.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ZenithDelay {
    /// Dry / hydrostatic component in metres.
    pub dry_m: f64,
    /// Wet component in metres.
    pub wet_m: f64,
}

impl ZenithDelay {
    /// Total zenith delay in metres.
    #[must_use]
    pub const fn total_m(&self) -> f64 {
        self.dry_m + self.wet_m
    }
}

/// Compute the Hopfield zenith delay given surface conditions.
///
/// `Δd = 155.2·10⁻⁷ · P / T · (h_d − h)` where `h_d = 40136 + 148.72·T_c`
/// is the top of the dry troposphere (`T_c` in °C).
///
/// `Δw = 155.2·10⁻⁷ · 4810 · e / T² · (h_w − h)` with `h_w ≈ 11000 m`.
#[must_use]
pub fn zenith_delay(cond: SurfaceConditions) -> ZenithDelay {
    let t_c = cond.temperature_k - 273.15;
    let hd = 148.72_f64.mul_add(t_c, 40_136.0);
    let hw = 11_000.0_f64;
    let e = cond.water_vapour_pressure_hpa();

    let dry_m = 155.2e-7 * cond.pressure_hpa / cond.temperature_k * (hd - cond.height_m).max(0.0);
    let wet_m = 155.2e-7 * 4810.0 * e / (cond.temperature_k * cond.temperature_k)
        * (hw - cond.height_m).max(0.0);
    ZenithDelay { dry_m, wet_m }
}

// ---------------------------------------------------------------------------
// Slant mapping
// ---------------------------------------------------------------------------

/// Simple 1/sin(e) mapping function with lower elevation cap at 3°.
///
/// Production PPP engines use the Vienna Mapping Function 1 (`VMF1`) or
/// Global Mapping Function (`GMF`) instead.
#[must_use]
pub fn simple_mapping(elev_rad: f64) -> f64 {
    let min = 3.0_f64.to_radians();
    let e = elev_rad.max(min);
    1.0 / e.sin()
}

/// Slant delay in metres for the given zenith delay and elevation angle
/// (radians).
#[must_use]
pub fn slant_delay_m(zd: ZenithDelay, elev_rad: f64) -> f64 {
    zd.total_m() * simple_mapping(elev_rad)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn standard_atmosphere_zenith_delay_matches_2_3_m() {
        let cond = SurfaceConditions::standard_atmosphere();
        let zd = zenith_delay(cond);
        // Textbook Hopfield gives ~2.3 m at sea-level standard atmosphere.
        assert!((zd.total_m() - 2.3).abs() < 0.4, "got {}", zd.total_m());
    }

    #[test]
    fn zenith_delay_dry_dominates_wet_by_order_of_magnitude() {
        let cond = SurfaceConditions::standard_atmosphere();
        let zd = zenith_delay(cond);
        assert!(zd.dry_m > 10.0 * zd.wet_m);
    }

    #[test]
    fn saturation_pressure_positive_at_room_temperature() {
        let cond = SurfaceConditions {
            temperature_k: 293.15,
            ..SurfaceConditions::standard_atmosphere()
        };
        assert!(cond.saturation_pressure_hpa() > 20.0);
    }

    #[test]
    fn zero_humidity_gives_zero_wet_delay() {
        let cond = SurfaceConditions {
            relative_humidity_percent: 0.0,
            ..SurfaceConditions::standard_atmosphere()
        };
        let zd = zenith_delay(cond);
        assert!((zd.wet_m - 0.0).abs() < 1e-9);
    }

    #[test]
    fn increasing_pressure_increases_dry_delay() {
        let a = SurfaceConditions::standard_atmosphere();
        let mut b = a;
        b.pressure_hpa += 50.0;
        assert!(zenith_delay(b).dry_m > zenith_delay(a).dry_m);
    }

    #[test]
    fn slant_delay_grows_with_low_elevation() {
        let zd = zenith_delay(SurfaceConditions::standard_atmosphere());
        let low = slant_delay_m(zd, 0.1);
        let high = slant_delay_m(zd, 1.5);
        assert!(low > high);
    }

    #[test]
    fn mapping_function_at_zenith_is_one() {
        let m = simple_mapping(std::f64::consts::FRAC_PI_2);
        assert!((m - 1.0).abs() < 1e-9);
    }

    #[test]
    fn mapping_function_clamps_low_elevation() {
        // Below 3° we cap at 3° so mapping cannot explode.
        let below_cap = simple_mapping(0.001);
        let at_cap = simple_mapping(3.0_f64.to_radians());
        assert!((below_cap - at_cap).abs() < 1e-9);
    }

    #[test]
    fn zenith_delay_zero_at_high_altitude() {
        let cond = SurfaceConditions {
            height_m: 50_000.0,
            ..SurfaceConditions::standard_atmosphere()
        };
        let zd = zenith_delay(cond);
        assert!((zd.dry_m - 0.0).abs() < 1e-9);
        assert!((zd.wet_m - 0.0).abs() < 1e-9);
    }

    #[test]
    fn total_m_equals_sum_of_components() {
        let zd = ZenithDelay {
            dry_m: 2.0,
            wet_m: 0.3,
        };
        assert!((zd.total_m() - 2.3).abs() < 1e-12);
    }
}
