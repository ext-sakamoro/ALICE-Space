//! Orbital mechanics primitives.

use std::f64::consts::PI;

/// Unique body identifier.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct BodyId(pub u64);

/// Classical Keplerian orbital elements.
#[derive(Debug, Clone)]
pub struct OrbitalElements {
    pub semi_major_axis_km: f64,
    pub eccentricity: f64,
    pub inclination_rad: f64,
    pub raan_rad: f64,
    pub arg_periapsis_rad: f64,
    pub true_anomaly_rad: f64,
}

/// A celestial body with gravitational parameter.
#[derive(Debug, Clone)]
pub struct CelestialBody {
    pub id: BodyId,
    pub name_hash: u64,
    pub mass_kg: f64,
    pub radius_km: f64,
    /// Gravitational parameter μ = G*M (km³/s²).
    pub mu: f64,
}

/// Spacecraft position and velocity state vector.
#[derive(Debug, Clone)]
pub struct SpacecraftState {
    pub position_km: [f64; 3],
    pub velocity_km_s: [f64; 3],
    pub timestamp_ns: u64,
    pub fuel_kg: f64,
}

use crate::fnv1a;

impl CelestialBody {
    #[must_use]
    pub fn new(id: u64, name: &str, mass_kg: f64, radius_km: f64, mu: f64) -> Self {
        Self {
            id: BodyId(id),
            name_hash: fnv1a(name.as_bytes()),
            mass_kg,
            radius_km,
            mu,
        }
    }
}

/// Orbital period T = 2π√(a³/μ) in seconds.
#[inline]
#[must_use]
pub fn orbital_period(semi_major_axis_km: f64, mu: f64) -> f64 {
    2.0 * PI * (semi_major_axis_km.powi(3) / mu).sqrt()
}

/// Vis-viva velocity: v = √(μ*(2/r - 1/a)) in km/s.
#[inline]
#[must_use]
pub fn orbital_velocity(r_km: f64, a_km: f64, mu: f64) -> f64 {
    (mu * (2.0 / r_km - 1.0 / a_km)).sqrt()
}

/// Hohmann transfer delta-v: returns (dv1, dv2) in km/s.
#[must_use]
pub fn delta_v_hohmann(r1_km: f64, r2_km: f64, mu: f64) -> (f64, f64) {
    let a_transfer = f64::midpoint(r1_km, r2_km);
    let v1_circular = (mu / r1_km).sqrt();
    let v1_transfer = orbital_velocity(r1_km, a_transfer, mu);
    let dv1 = (v1_transfer - v1_circular).abs();

    let v2_circular = (mu / r2_km).sqrt();
    let v2_transfer = orbital_velocity(r2_km, a_transfer, mu);
    let dv2 = (v2_circular - v2_transfer).abs();

    (dv1, dv2)
}

/// One-way light delay in seconds. c = 299792.458 km/s.
#[inline]
#[must_use]
pub fn light_delay_s(distance_km: f64) -> f64 {
    distance_km / 299_792.458
}

#[cfg(test)]
mod tests {
    use super::*;

    const MU_EARTH: f64 = 398_600.441_8; // km³/s²

    #[test]
    fn iss_orbital_period() {
        // ISS: ~408 km altitude → a ≈ 6778 km, T ≈ 92.3 min
        let t = orbital_period(6778.0, MU_EARTH);
        assert!((t / 60.0 - 92.3).abs() < 1.0);
    }

    #[test]
    fn circular_velocity_leo() {
        // LEO at 400 km: v ≈ 7.67 km/s
        let v = orbital_velocity(6778.0, 6778.0, MU_EARTH);
        assert!((v - 7.67).abs() < 0.1);
    }

    #[test]
    fn hohmann_leo_to_geo() {
        // LEO (6578 km) to GEO (42164 km)
        let (dv1, dv2) = delta_v_hohmann(6578.0, 42164.0, MU_EARTH);
        // dv1 ≈ 2.46 km/s, dv2 ≈ 1.48 km/s, total ≈ 3.94 km/s
        let total = dv1 + dv2;
        assert!((total - 3.94).abs() < 0.1);
    }

    #[test]
    fn light_delay_earth_moon() {
        // Earth-Moon: 384400 km → ~1.28 s
        let d = light_delay_s(384_400.0);
        assert!((d - 1.28).abs() < 0.01);
    }

    #[test]
    fn light_delay_earth_mars() {
        // Closest approach: ~55.7 million km → ~186 s (~3.1 min)
        let d = light_delay_s(55_700_000.0);
        assert!((d / 60.0 - 3.1).abs() < 0.1);
    }

    #[test]
    fn celestial_body_hash() {
        let earth = CelestialBody::new(3, "Earth", 5.972e24, 6371.0, MU_EARTH);
        let mars = CelestialBody::new(4, "Mars", 6.39e23, 3389.5, 42828.37);
        assert_ne!(earth.name_hash, mars.name_hash);
    }

    #[test]
    fn celestial_body_hash_deterministic() {
        let e1 = CelestialBody::new(3, "Earth", 5.972e24, 6371.0, MU_EARTH);
        let e2 = CelestialBody::new(3, "Earth", 5.972e24, 6371.0, MU_EARTH);
        assert_eq!(e1.name_hash, e2.name_hash);
        assert_ne!(e1.name_hash, 0);
    }

    #[test]
    fn orbital_period_geo() {
        // GEO: a = 42164 km, T ≈ 86164 s (sidereal day)
        let t = orbital_period(42164.0, MU_EARTH);
        assert!((t - 86164.0).abs() < 100.0, "GEO period: {t} s");
    }

    #[test]
    fn orbital_velocity_increases_closer_to_body() {
        let v_near = orbital_velocity(6578.0, 6578.0, MU_EARTH);
        let v_far = orbital_velocity(42164.0, 42164.0, MU_EARTH);
        assert!(v_near > v_far, "LEO v={v_near} > GEO v={v_far}");
    }

    #[test]
    fn hohmann_symmetric_transfer() {
        // Hohmann from r1 to r2 should give same total dv regardless of direction
        let (dv1_up, dv2_up) = delta_v_hohmann(6578.0, 42164.0, MU_EARTH);
        let (dv1_down, dv2_down) = delta_v_hohmann(42164.0, 6578.0, MU_EARTH);
        let total_up = dv1_up + dv2_up;
        let total_down = dv1_down + dv2_down;
        assert!(
            (total_up - total_down).abs() < 0.01,
            "up={total_up} down={total_down}"
        );
    }

    #[test]
    fn hohmann_zero_transfer() {
        // Same orbit: dv should be ~0
        let (dv1, dv2) = delta_v_hohmann(6778.0, 6778.0, MU_EARTH);
        assert!(dv1 < 1e-10, "dv1={dv1}");
        assert!(dv2 < 1e-10, "dv2={dv2}");
    }

    #[test]
    fn light_delay_zero_distance() {
        let d = light_delay_s(0.0);
        assert!((d).abs() < 1e-15);
    }

    #[test]
    fn light_delay_earth_saturn() {
        // Earth-Saturn closest ~1.2 billion km → ~4000 s (~67 min)
        let d = light_delay_s(1_200_000_000.0);
        assert!(
            (d / 60.0 - 66.7).abs() < 1.0,
            "Saturn delay: {} min",
            d / 60.0
        );
    }

    #[test]
    fn spacecraft_state_clone() {
        let s = SpacecraftState {
            position_km: [1.0, 2.0, 3.0],
            velocity_km_s: [4.0, 5.0, 6.0],
            timestamp_ns: 42,
            fuel_kg: 100.0,
        };
        let s2 = s.clone();
        assert_eq!(s.position_km, s2.position_km);
        assert_eq!(s.velocity_km_s, s2.velocity_km_s);
        assert_eq!(s.timestamp_ns, s2.timestamp_ns);
        assert_eq!(s.fuel_kg, s2.fuel_kg);
    }

    #[test]
    fn orbital_elements_eccentricity_bounds() {
        // Parabolic escape: e = 1.0 means v = escape velocity at periapsis
        let elements = OrbitalElements {
            semi_major_axis_km: 6778.0,
            eccentricity: 0.0,
            inclination_rad: 0.0,
            raan_rad: 0.0,
            arg_periapsis_rad: 0.0,
            true_anomaly_rad: 0.0,
        };
        assert!(elements.eccentricity >= 0.0);
    }

    #[test]
    fn body_id_equality() {
        let a = BodyId(42);
        let b = BodyId(42);
        let c = BodyId(43);
        assert_eq!(a, b);
        assert_ne!(a, c);
    }

    #[test]
    fn orbital_velocity_vis_viva_elliptical() {
        // 楕円軌道: 近地点でr<aなので速度は円軌道より大きい
        let r = 6578.0; // 近地点
        let a = 24371.0; // 楕円軌道の半長軸
        let v_ellip = orbital_velocity(r, a, MU_EARTH);
        let v_circ = orbital_velocity(r, r, MU_EARTH);
        assert!(v_ellip > v_circ, "楕円近地点速度 > 円軌道速度");
    }

    #[test]
    fn orbital_period_kepler_third_law() {
        // ケプラー第三法則: T² ∝ a³
        let a1 = 7000.0;
        let a2 = 14000.0; // 2倍の半長軸
        let t1 = orbital_period(a1, MU_EARTH);
        let t2 = orbital_period(a2, MU_EARTH);
        // T2/T1 = (a2/a1)^(3/2) = 2^1.5 ≈ 2.828
        let ratio = t2 / t1;
        assert!(
            (ratio - 2.828).abs() < 0.01,
            "Kepler 3rd law ratio: {ratio}"
        );
    }

    #[test]
    fn hohmann_dv1_positive_outward() {
        // 外側への遷移: dv1は正（加速）
        let (dv1, dv2) = delta_v_hohmann(6578.0, 42164.0, MU_EARTH);
        assert!(dv1 > 0.0);
        assert!(dv2 > 0.0);
    }

    #[test]
    fn light_delay_voyager_distance() {
        // ボイジャー1号: ~24 billion km → ~22 hours
        let d = light_delay_s(24_000_000_000.0);
        let hours = d / 3600.0;
        assert!((hours - 22.2).abs() < 1.0, "Voyager delay: {hours} hours");
    }

    #[test]
    fn orbital_elements_clone() {
        let e = OrbitalElements {
            semi_major_axis_km: 7000.0,
            eccentricity: 0.1,
            inclination_rad: 0.5,
            raan_rad: 1.0,
            arg_periapsis_rad: 2.0,
            true_anomaly_rad: 3.0,
        };
        let e2 = e.clone();
        assert_eq!(e.semi_major_axis_km, e2.semi_major_axis_km);
        assert_eq!(e.eccentricity, e2.eccentricity);
    }

    mod prop {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            #[test]
            fn orbital_period_positive(
                a in 6400.0f64..100_000.0,
                mu in 1000.0f64..500_000.0,
            ) {
                let t = orbital_period(a, mu);
                prop_assert!(t > 0.0, "period must be positive: {t}");
                prop_assert!(t.is_finite());
            }

            #[test]
            fn light_delay_non_negative(d in 0.0f64..1e12) {
                let delay = light_delay_s(d);
                prop_assert!(delay >= 0.0);
                prop_assert!(delay.is_finite());
            }

            #[test]
            fn hohmann_total_dv_symmetric(
                r1 in 6400.0f64..50000.0,
                r2 in 6400.0f64..50000.0,
            ) {
                let (dv1a, dv2a) = delta_v_hohmann(r1, r2, MU_EARTH);
                let (dv1b, dv2b) = delta_v_hohmann(r2, r1, MU_EARTH);
                let total_a = dv1a + dv2a;
                let total_b = dv1b + dv2b;
                prop_assert!((total_a - total_b).abs() < 0.01,
                    "asymmetric: {total_a} vs {total_b}");
            }

            #[test]
            fn orbital_velocity_positive(
                r in 6400.0f64..50000.0,
            ) {
                // Circular orbit: r == a → always valid
                let v = orbital_velocity(r, r, MU_EARTH);
                prop_assert!(v > 0.0);
                prop_assert!(v.is_finite());
            }
        }
    }
}
