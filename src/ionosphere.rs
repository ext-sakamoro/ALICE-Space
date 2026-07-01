//! Ionospheric delay correction — Klobuchar single-frequency model.
//!
//! Implements the Klobuchar model as specified in `IS-GPS-200` Section 20.3.3.5.2.5:
//! given the 8 broadcast coefficients (`α0..α3, β0..β3`), the user's geodetic
//! position, the line-of-sight to a satellite, and the observation time,
//! compute the extra pseudorange delay caused by the ionosphere.
//!
//! The model is used by single-frequency GPS/QZSS receivers and by every
//! signal-quality tool that needs to isolate ionospheric anomalies (a common
//! spoofing indicator).

use core::f64::consts::PI;

// ---------------------------------------------------------------------------
// Coefficients broadcast in GPS/QZSS navigation message
// ---------------------------------------------------------------------------

/// Klobuchar coefficients broadcast in the GPS/QZSS navigation message.
///
/// `alpha` are the four cosine-amplitude terms (seconds); `beta` are the
/// four cosine-period terms (seconds). Values are typically decoded from
/// LNAV subframe 4 page 18 (GPS) or the corresponding QZSS L1 C/A page.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct KlobucharCoefficients {
    pub alpha: [f64; 4],
    pub beta: [f64; 4],
}

impl KlobucharCoefficients {
    /// Constructor.
    #[must_use]
    pub const fn new(alpha: [f64; 4], beta: [f64; 4]) -> Self {
        Self { alpha, beta }
    }

    /// Return the mid-2020 broadcast set for demo/testing purposes. These
    /// are representative real values but should not be relied on for
    /// operational corrections — decode the live navigation message instead.
    #[must_use]
    pub const fn representative_2020() -> Self {
        Self {
            alpha: [1.164e-8, 2.235e-8, -5.960e-8, -1.192e-7],
            beta: [1.229e5, 1.310e5, -6.554e4, -5.243e5],
        }
    }
}

// ---------------------------------------------------------------------------
// Klobuchar model
// ---------------------------------------------------------------------------

/// Speed of light in vacuum (m/s), used to convert seconds to meters.
pub const SPEED_OF_LIGHT_MS: f64 = 299_792_458.0;

/// Compute the ionospheric delay in seconds for the given geometry.
///
/// Arguments:
///
/// - `coeffs`: current Klobuchar coefficients.
/// - `user_lat_rad`, `user_lon_rad`: user geodetic latitude and longitude (radians).
/// - `azimuth_rad`, `elevation_rad`: line-of-sight to the satellite in the
///   user's local ENU frame (radians). Elevation must be positive.
/// - `gps_time_sec`: GPS time of week in seconds.
///
/// Returns the additional propagation delay in seconds. Convert to metres
/// by multiplying by [`SPEED_OF_LIGHT_MS`] or use [`klobuchar_delay_m`].
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn klobuchar_delay_s(
    coeffs: KlobucharCoefficients,
    user_lat_rad: f64,
    user_lon_rad: f64,
    azimuth_rad: f64,
    elevation_rad: f64,
    gps_time_sec: f64,
) -> f64 {
    // 1. Earth-centred angle (semicircles).
    let user_lat_sc = user_lat_rad / PI;
    let user_lon_sc = user_lon_rad / PI;
    let el_sc = elevation_rad / PI;
    let az_sc = azimuth_rad / PI;

    let psi = 0.0137 / (el_sc + 0.11) - 0.022;

    // 2. Sub-ionospheric latitude (semicircles), clamped to ±0.416.
    let lat_i = (user_lat_sc + psi * azimuth_rad.cos()).clamp(-0.416, 0.416);

    // 3. Sub-ionospheric longitude (semicircles).
    let denom = (lat_i * PI).cos();
    let lon_i = user_lon_sc + psi * azimuth_rad.sin() / denom;
    let _ = az_sc; // Not used directly (kept for spec alignment).

    // 4. Geomagnetic latitude of the sub-ionospheric point (semicircles).
    let lat_m = lat_i + 0.064 * ((lon_i - 1.617) * PI).cos();

    // 5. Local time at the sub-ionospheric point (seconds), 0 <= t < 86400.
    let mut t = 43_200.0 * lon_i + gps_time_sec;
    t = t.rem_euclid(86_400.0);

    // 6. Slant factor (dimensionless).
    let f = 1.0 + 16.0 * (0.53 - el_sc).powi(3);

    // 7. Amplitude of the cosine term (seconds).
    let mut amp = coeffs.alpha[0]
        + coeffs.alpha[1] * lat_m
        + coeffs.alpha[2] * lat_m * lat_m
        + coeffs.alpha[3] * lat_m * lat_m * lat_m;
    if amp < 0.0 {
        amp = 0.0;
    }

    // 8. Period of the cosine term (seconds), minimum 72 000.
    let mut per = coeffs.beta[0]
        + coeffs.beta[1] * lat_m
        + coeffs.beta[2] * lat_m * lat_m
        + coeffs.beta[3] * lat_m * lat_m * lat_m;
    if per < 72_000.0 {
        per = 72_000.0;
    }

    // 9. Phase (radians).
    let x = 2.0 * PI * (t - 50_400.0) / per;

    // 10. Ionospheric delay (seconds).
    if x.abs() < 1.57 {
        f * (5e-9 + amp * (1.0 - x * x / 2.0 + x.powi(4) / 24.0))
    } else {
        f * 5e-9
    }
}

/// Convenience: [`klobuchar_delay_s`] converted to metres by multiplying by
/// the speed of light.
#[must_use]
#[allow(clippy::too_many_arguments)]
pub fn klobuchar_delay_m(
    coeffs: KlobucharCoefficients,
    user_lat_rad: f64,
    user_lon_rad: f64,
    azimuth_rad: f64,
    elevation_rad: f64,
    gps_time_sec: f64,
) -> f64 {
    klobuchar_delay_s(
        coeffs,
        user_lat_rad,
        user_lon_rad,
        azimuth_rad,
        elevation_rad,
        gps_time_sec,
    ) * SPEED_OF_LIGHT_MS
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn tokyo_user() -> (f64, f64) {
        (35.6895_f64.to_radians(), 139.6917_f64.to_radians())
    }

    #[test]
    fn delay_is_positive_for_typical_geometry() {
        let (lat, lon) = tokyo_user();
        let d = klobuchar_delay_s(
            KlobucharCoefficients::representative_2020(),
            lat,
            lon,
            45.0_f64.to_radians(),
            30.0_f64.to_radians(),
            30_000.0,
        );
        assert!(d > 0.0);
    }

    #[test]
    fn delay_is_larger_at_low_elevation() {
        let (lat, lon) = tokyo_user();
        let coeffs = KlobucharCoefficients::representative_2020();
        let high = klobuchar_delay_s(coeffs, lat, lon, 0.0, 80.0_f64.to_radians(), 43_200.0);
        let low = klobuchar_delay_s(coeffs, lat, lon, 0.0, 10.0_f64.to_radians(), 43_200.0);
        assert!(
            low > high,
            "low-elevation delay should exceed high-elevation delay: high={high} low={low}"
        );
    }

    #[test]
    fn delay_matches_meter_conversion() {
        let (lat, lon) = tokyo_user();
        let coeffs = KlobucharCoefficients::representative_2020();
        let s = klobuchar_delay_s(coeffs, lat, lon, 0.5, 0.6, 40_000.0);
        let m = klobuchar_delay_m(coeffs, lat, lon, 0.5, 0.6, 40_000.0);
        assert!((m - s * SPEED_OF_LIGHT_MS).abs() < 1e-9);
    }

    #[test]
    fn delay_uses_night_floor() {
        // Nighttime (all-zero coefficients) saturates the delay at the
        // 5 ns floor scaled by the slant factor. At near-zenith the slant
        // factor is very close to 1, so the delay is ~5 ns.
        let (lat, lon) = tokyo_user();
        let coeffs = KlobucharCoefficients {
            alpha: [0.0; 4],
            beta: [0.0; 4],
        };
        let d = klobuchar_delay_s(coeffs, lat, lon, 0.0, 90.0_f64.to_radians(), 0.0);
        // 5 ns floor with slant factor ≈ 1.0004 gives d ≈ 5.002 ns.
        assert!((d - 5e-9).abs() < 1e-11, "d = {d}");
    }
}
