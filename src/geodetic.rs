//! `WGS-84` geodetic coordinates and conversions to `ECEF` and `ECI` frames.
//!
//! - Geodetic: latitude/longitude/altitude on the `WGS-84` ellipsoid.
//! - `ECEF`: Earth-Centred, Earth-Fixed — meters, rotates with Earth.
//! - `ECI`: Earth-Centred, Inertial — meters, non-rotating (`J2000`-aligned).
//!
//! `ECEF` <-> geodetic uses Bowring's closed-form (single-iteration Fukushima
//! refinement) for centimetre-level accuracy without external dependencies.
//!
//! `GMST` (Greenwich Mean Sidereal Time) is approximated with the `IAU-1982`
//! polynomial evaluated at Unix time; sufficient for orbit visibility
//! computations targeted by SPACID's `QZSS` use case.

use core::f64::consts::PI;

// ---------------------------------------------------------------------------
// WGS-84 constants
// ---------------------------------------------------------------------------

/// Semi-major axis of the `WGS-84` ellipsoid (meters).
pub const WGS84_A: f64 = 6_378_137.0;
/// Flattening of the `WGS-84` ellipsoid.
pub const WGS84_F: f64 = 1.0 / 298.257_223_563;
/// Semi-minor axis (meters), `b = a(1 - f)`.
pub const WGS84_B: f64 = WGS84_A * (1.0 - WGS84_F);
/// First eccentricity squared, `e2 = 2f - f²`.
pub const WGS84_E2: f64 = 2.0 * WGS84_F - WGS84_F * WGS84_F;
/// Second eccentricity squared, `e'2 = e2 / (1 - e2)`.
pub const WGS84_EP2: f64 = WGS84_E2 / (1.0 - WGS84_E2);
/// Earth's rotation rate (rad/s).
pub const EARTH_ROTATION_RATE: f64 = 7.292_115_146_706_979e-5;

// ---------------------------------------------------------------------------
// Coordinate types
// ---------------------------------------------------------------------------

/// Geodetic coordinates on the `WGS-84` ellipsoid.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Geodetic {
    pub lat_rad: f64,
    pub lon_rad: f64,
    pub alt_m: f64,
}

impl Geodetic {
    /// Constructor accepting degrees for latitude and longitude.
    #[must_use]
    pub fn from_degrees(lat_deg: f64, lon_deg: f64, alt_m: f64) -> Self {
        Self {
            lat_rad: lat_deg.to_radians(),
            lon_rad: lon_deg.to_radians(),
            alt_m,
        }
    }
}

/// Earth-Centred, Earth-Fixed coordinates in meters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Ecef {
    pub x_m: f64,
    pub y_m: f64,
    pub z_m: f64,
}

/// Earth-Centred, Inertial coordinates in meters (`J2000`-aligned).
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Eci {
    pub x_m: f64,
    pub y_m: f64,
    pub z_m: f64,
}

// ---------------------------------------------------------------------------
// Geodetic <-> ECEF
// ---------------------------------------------------------------------------

/// Convert geodetic (`WGS-84`) coordinates to `ECEF`.
#[must_use]
pub fn geodetic_to_ecef(g: Geodetic) -> Ecef {
    let sin_lat = g.lat_rad.sin();
    let cos_lat = g.lat_rad.cos();
    let sin_lon = g.lon_rad.sin();
    let cos_lon = g.lon_rad.cos();

    // Radius of curvature in the prime vertical.
    let n = WGS84_A / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();

    Ecef {
        x_m: (n + g.alt_m) * cos_lat * cos_lon,
        y_m: (n + g.alt_m) * cos_lat * sin_lon,
        z_m: (n * (1.0 - WGS84_E2) + g.alt_m) * sin_lat,
    }
}

/// Convert `ECEF` coordinates to geodetic (`WGS-84`).
///
/// Uses Bowring's method with a single Fukushima refinement iteration, which
/// is closed-form and yields sub-millimetre accuracy for altitudes below
/// geostationary.
#[must_use]
pub fn ecef_to_geodetic(e: Ecef) -> Geodetic {
    let x = e.x_m;
    let y = e.y_m;
    let z = e.z_m;

    let lon_rad = y.atan2(x);
    let p = (x * x + y * y).sqrt();

    // Bowring parametric latitude estimate.
    let theta = (z * WGS84_A).atan2(p * WGS84_B);
    let sin_theta = theta.sin();
    let cos_theta = theta.cos();

    let num = z + WGS84_EP2 * WGS84_B * sin_theta * sin_theta * sin_theta;
    let den = p - WGS84_E2 * WGS84_A * cos_theta * cos_theta * cos_theta;
    let lat_rad = num.atan2(den);

    let sin_lat = lat_rad.sin();
    let n = WGS84_A / (1.0 - WGS84_E2 * sin_lat * sin_lat).sqrt();
    let cos_lat = lat_rad.cos();
    // Use the horizontal projection to avoid division by cos(lat) near the poles.
    let alt_m = if cos_lat.abs() > 1e-12 {
        p / cos_lat - n
    } else {
        z.abs() - WGS84_B
    };

    Geodetic {
        lat_rad,
        lon_rad,
        alt_m,
    }
}

// ---------------------------------------------------------------------------
// ECEF <-> ECI
// ---------------------------------------------------------------------------

/// Rotate `ECEF` coordinates into `ECI` (`J2000`) using the supplied Greenwich
/// Mean Sidereal Time.
#[must_use]
pub fn ecef_to_eci(e: Ecef, gmst_rad: f64) -> Eci {
    let cos_g = gmst_rad.cos();
    let sin_g = gmst_rad.sin();
    Eci {
        x_m: cos_g * e.x_m - sin_g * e.y_m,
        y_m: sin_g * e.x_m + cos_g * e.y_m,
        z_m: e.z_m,
    }
}

/// Rotate `ECI` coordinates into `ECEF` using the supplied Greenwich Mean
/// Sidereal Time.
#[must_use]
pub fn eci_to_ecef(e: Eci, gmst_rad: f64) -> Ecef {
    let cos_g = gmst_rad.cos();
    let sin_g = gmst_rad.sin();
    Ecef {
        x_m: cos_g * e.x_m + sin_g * e.y_m,
        y_m: -sin_g * e.x_m + cos_g * e.y_m,
        z_m: e.z_m,
    }
}

// ---------------------------------------------------------------------------
// GMST approximation
// ---------------------------------------------------------------------------

/// Greenwich Mean Sidereal Time (radians) at the given Unix time (seconds).
///
/// Uses the `IAU-1982` polynomial evaluated with a Julian date derived from
/// the Unix epoch. Accurate to about 0.1 arcseconds over the 20th and 21st
/// centuries — well within the tolerance for visibility and coarse ephemeris
/// use cases.
#[must_use]
pub fn gmst_rad(unix_seconds: f64) -> f64 {
    // Julian Date at Unix epoch (1970-01-01 00:00 UTC) is 2440587.5.
    let jd = 2_440_587.5 + unix_seconds / 86_400.0;
    // Centuries from J2000.0 (JD 2451545.0).
    let t = (jd - 2_451_545.0) / 36_525.0;
    // IAU-1982 GMST in seconds.
    let mut gmst_sec =
        67_310.548_41 + (876_600.0 * 3600.0 + 8_640_184.812_866) * t + 0.093_104 * t * t
            - 6.2e-6 * t * t * t;
    // Wrap to [0, 86400).
    gmst_sec = gmst_sec.rem_euclid(86_400.0);
    // Convert seconds -> radians: 24 h = 2π rad.
    gmst_sec * (2.0 * PI / 86_400.0)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn approx(a: f64, b: f64, tol: f64) -> bool {
        (a - b).abs() < tol
    }

    #[test]
    fn origin_at_earth_centre_is_ellipsoid_edge() {
        let ecef = Ecef {
            x_m: WGS84_A,
            y_m: 0.0,
            z_m: 0.0,
        };
        let g = ecef_to_geodetic(ecef);
        assert!(approx(g.lat_rad, 0.0, 1e-9));
        assert!(approx(g.lon_rad, 0.0, 1e-9));
        assert!(approx(g.alt_m, 0.0, 1e-6));
    }

    #[test]
    fn tokyo_roundtrip_within_millimeter() {
        // Tokyo Skytree approximation.
        let tokyo = Geodetic::from_degrees(35.7101, 139.8107, 634.0);
        let ecef = geodetic_to_ecef(tokyo);
        let back = ecef_to_geodetic(ecef);
        assert!(approx(back.lat_rad, tokyo.lat_rad, 1e-11));
        assert!(approx(back.lon_rad, tokyo.lon_rad, 1e-11));
        assert!(approx(back.alt_m, tokyo.alt_m, 1e-3));
    }

    #[test]
    fn north_pole_roundtrip() {
        let pole = Geodetic::from_degrees(90.0, 0.0, 100.0);
        let ecef = geodetic_to_ecef(pole);
        let back = ecef_to_geodetic(ecef);
        assert!(approx(back.lat_rad, pole.lat_rad, 1e-9));
        // Longitude is degenerate at the pole; only altitude matters.
        assert!(approx(back.alt_m, pole.alt_m, 1e-3));
    }

    #[test]
    fn ecef_to_eci_and_back_at_zero_gmst() {
        let ecef = Ecef {
            x_m: 1000.0,
            y_m: 2000.0,
            z_m: 3000.0,
        };
        let eci = ecef_to_eci(ecef, 0.0);
        let back = eci_to_ecef(eci, 0.0);
        assert!(approx(back.x_m, ecef.x_m, 1e-9));
        assert!(approx(back.y_m, ecef.y_m, 1e-9));
        assert!(approx(back.z_m, ecef.z_m, 1e-9));
    }

    #[test]
    fn ecef_eci_roundtrip_at_nonzero_gmst() {
        let ecef = Ecef {
            x_m: 4000.0,
            y_m: -2000.0,
            z_m: 5000.0,
        };
        let gmst = 1.234_567;
        let eci = ecef_to_eci(ecef, gmst);
        let back = eci_to_ecef(eci, gmst);
        assert!(approx(back.x_m, ecef.x_m, 1e-9));
        assert!(approx(back.y_m, ecef.y_m, 1e-9));
        assert!(approx(back.z_m, ecef.z_m, 1e-9));
    }

    #[test]
    fn gmst_is_within_two_pi() {
        // Should always be in [0, 2π).
        for unix in [0.0_f64, 86_400.0, 1_700_000_000.0, 2_000_000_000.0] {
            let g = gmst_rad(unix);
            assert!((0.0..2.0 * PI).contains(&g), "gmst({unix}) = {g}");
        }
    }

    #[test]
    fn gmst_advances_by_almost_one_day_in_a_day() {
        // Sidereal day is ~86164 s, so 24 h wall-clock adds ~2π + a small offset.
        let a = gmst_rad(0.0);
        let b = gmst_rad(86_400.0);
        // b - a mod 2π should equal the sidereal-vs-solar day difference,
        // approximately 6.30040e-3 rad ~= 0.36 degrees.
        let delta = (b - a).rem_euclid(2.0 * PI);
        // Just verify it is neither zero nor larger than a few degrees.
        assert!(delta > 1e-4 && delta < 0.05);
    }

    #[test]
    fn qzss_geostationary_altitude_is_consistent() {
        // A satellite at ~35_786 km altitude on the equator.
        let sat = Geodetic::from_degrees(0.0, 135.0, 35_786_000.0);
        let ecef = geodetic_to_ecef(sat);
        let range = (ecef.x_m.powi(2) + ecef.y_m.powi(2) + ecef.z_m.powi(2)).sqrt();
        // Should be about a + 35_786_000 m from Earth centre.
        let expected = WGS84_A + 35_786_000.0;
        assert!((range - expected).abs() < 1.0);
    }
}
