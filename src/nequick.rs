//! Simplified NeQuick G ionosphere model (Galileo Reduced Broadcast).
//!
//! The full NeQuick 2 model integrates the Chapman-layer electron density
//! along a slant path and requires solar flux and F2-layer coefficients
//! from CCIR maps. NeQuick G ("Galileo") supplies a broadcast reduced
//! form using three effective ionisation level coefficients (a0, a1, a2)
//! that Galileo transmits in the F/NAV message.
//!
//! This module implements the **single-layer slant TEC + delay** portion
//! of NeQuick G. It is dimensionally correct and suitable for engineering
//! studies, unit testing of downstream code, and comparison with
//! Klobuchar in the [`crate::ionosphere`] module. A production PPP
//! engine should link the ITU-R Galileo reference source.
//!
//! # References
//!
//! - European Commission (2016), "European GNSS (Galileo) Open Service
//!   Ionospheric Correction Algorithm for Galileo Single Frequency Users"
//!   Issue 1.2, §2.5.
//! - Nava, B., Coïsson, P., & Radicella, S. M. (2008), "A new version of
//!   the NeQuick ionosphere electron density model", J. Atmos. Sol.-Terr.
//!   Phys., 70(15), 1856-1862.

#![allow(clippy::too_long_first_doc_paragraph)]

/// NeQuick G broadcast coefficients (three effective ionisation levels).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct NequickGCoefficients {
    /// `a0` — constant term (solar-flux-units, SFU).
    pub a0: f64,
    /// `a1` — dependence on modified dip latitude (SFU/deg).
    pub a1: f64,
    /// `a2` — quadratic dependence on modified dip latitude (SFU/deg²).
    pub a2: f64,
}

impl NequickGCoefficients {
    /// Construct a new set of coefficients.
    #[must_use]
    pub const fn new(a0: f64, a1: f64, a2: f64) -> Self {
        Self { a0, a1, a2 }
    }

    /// A representative broadcast set (~solar maximum 2014-2015).
    #[must_use]
    pub const fn representative_2014() -> Self {
        Self {
            a0: 236.831641,
            a1: -0.39362878,
            a2: 0.00402826613,
        }
    }

    /// Effective ionisation level at the given modified dip latitude
    /// (`µ`, degrees). Az = a0 + a1·µ + a2·µ².
    #[must_use]
    pub fn az_sfu(&self, mu_deg: f64) -> f64 {
        self.a2
            .mul_add(mu_deg * mu_deg, self.a1.mul_add(mu_deg, self.a0))
    }
}

// ---------------------------------------------------------------------------
// Slant delay
// ---------------------------------------------------------------------------

/// Ionospheric slant delay in metres for a GNSS carrier at frequency
/// `f_hz`, given the NeQuick G effective ionisation level `az_sfu`,
/// the vertical total electron content model, and the receiver
/// elevation `elev_rad`.
///
/// Uses the standard obliquity factor `1 / cos(z')` where
/// `z' = arcsin(R_E · sin(z) / (R_E + H_ion))` with `R_E = 6371.2 km`
/// and shell height `H_ion = 350 km`.
#[must_use]
pub fn slant_delay_m(az_sfu: f64, elev_rad: f64, f_hz: f64) -> f64 {
    // Convert Az (10^22 electrons/m^2) to vertical TEC (TECU). NeQuick G
    // documentation uses a linear proxy: vTEC ≈ Az.
    let vtec_tecu = az_sfu.max(0.0);
    // Piercing point obliquity: R_E / (R_E + H).
    let r_e_km = 6371.2_f64;
    let h_ion_km = 350.0_f64;
    let zenith_rad = std::f64::consts::FRAC_PI_2 - elev_rad;
    let ratio = (r_e_km / (r_e_km + h_ion_km)) * zenith_rad.sin();
    let z_prime = ratio.asin();
    let obliquity = 1.0 / z_prime.cos();
    // 1 TECU introduces 0.162 m of delay on L1 at 1575.42 MHz (Kaplan
    // 2017 §5.3.3): delay_L1 = 40.31 · TEC / f^2, with f in Hz.
    let stec = vtec_tecu * obliquity;
    40.31 * stec * 1e16 / (f_hz * f_hz)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    const L1_HZ: f64 = 1_575_420_000.0;

    #[test]
    fn az_at_zero_lat_returns_a0() {
        let c = NequickGCoefficients::new(200.0, 1.0, 0.01);
        assert!((c.az_sfu(0.0) - 200.0).abs() < 1e-12);
    }

    #[test]
    fn az_is_quadratic_in_latitude() {
        let c = NequickGCoefficients::new(100.0, 0.0, 1.0);
        // Az(10) = 100 + 0 + 100 = 200
        assert!((c.az_sfu(10.0) - 200.0).abs() < 1e-9);
    }

    #[test]
    fn representative_2014_returns_positive_at_equator() {
        let c = NequickGCoefficients::representative_2014();
        assert!(c.az_sfu(0.0) > 0.0);
    }

    #[test]
    fn slant_delay_positive_at_zenith() {
        let d = slant_delay_m(100.0, std::f64::consts::FRAC_PI_2, L1_HZ);
        assert!(d > 0.0);
    }

    #[test]
    fn slant_delay_grows_with_low_elevation() {
        let high = slant_delay_m(100.0, 1.4, L1_HZ);
        let low = slant_delay_m(100.0, 0.2, L1_HZ);
        assert!(low > high);
    }

    #[test]
    fn zero_az_gives_zero_delay() {
        assert!((slant_delay_m(0.0, 0.5, L1_HZ) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn negative_az_is_clamped_to_zero() {
        assert!((slant_delay_m(-10.0, 0.5, L1_HZ) - 0.0).abs() < 1e-12);
    }

    #[test]
    fn delay_scales_inversely_with_frequency_squared() {
        let d1 = slant_delay_m(100.0, std::f64::consts::FRAC_PI_2, L1_HZ);
        let d2 = slant_delay_m(100.0, std::f64::consts::FRAC_PI_2, 2.0 * L1_HZ);
        assert!((d2 / d1 - 0.25).abs() < 1e-9);
    }

    #[test]
    fn az_sign_symmetry_with_only_quadratic_term() {
        let c = NequickGCoefficients::new(0.0, 0.0, 1.0);
        assert!((c.az_sfu(5.0) - c.az_sfu(-5.0)).abs() < 1e-12);
    }

    #[test]
    fn coefficients_equality() {
        let a = NequickGCoefficients::new(1.0, 2.0, 3.0);
        let b = NequickGCoefficients::new(1.0, 2.0, 3.0);
        assert_eq!(a, b);
    }
}
