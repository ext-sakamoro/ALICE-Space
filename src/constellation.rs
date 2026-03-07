//! Walker Delta Pattern constellation generator
//!
//! Generates uniform satellite constellation layouts using the
//! Walker Delta Pattern (T/P/F notation). Used for coverage analysis
//! of communication and Earth observation constellations.
//!
//! Author: Moroya Sakamoto

use crate::fnv1a;
use crate::orbit::OrbitalElements;
use std::f64::consts::PI;

/// Walker Delta constellation configuration.
///
/// Walker notation: i:T/P/F where
/// - i = inclination (rad)
/// - T = total number of satellites
/// - P = number of equally-spaced orbital planes
/// - F = phasing factor (0 <= F < P)
#[derive(Debug, Clone, Copy)]
pub struct WalkerConstellation {
    /// Total number of satellites.
    pub total_satellites: u32,
    /// Number of orbital planes.
    pub num_planes: u32,
    /// Relative phasing factor (0 <= F < P).
    pub phasing_factor: u32,
    /// Semi-major axis in km (all satellites share the same altitude).
    pub semi_major_axis_km: f64,
    /// Inclination in radians (all planes share the same inclination).
    pub inclination_rad: f64,
    /// Eccentricity (typically 0 for Walker constellations).
    pub eccentricity: f64,
}

/// A generated satellite with its orbital elements and identifiers.
#[derive(Debug, Clone)]
pub struct WalkerSatellite {
    /// Plane index (0-based).
    pub plane_index: u32,
    /// Satellite index within plane (0-based).
    pub sat_index: u32,
    /// Global satellite index (0-based).
    pub global_index: u32,
    /// Keplerian orbital elements.
    pub elements: OrbitalElements,
    /// Deterministic content hash.
    pub content_hash: u64,
}

impl WalkerConstellation {
    /// Create a new Walker Delta constellation.
    #[must_use]
    pub const fn new(total: u32, planes: u32, phasing: u32, sma_km: f64, inc_rad: f64) -> Self {
        Self {
            total_satellites: total,
            num_planes: planes,
            phasing_factor: if planes > 0 { phasing % planes } else { 0 },
            semi_major_axis_km: sma_km,
            inclination_rad: inc_rad,
            eccentricity: 0.0,
        }
    }

    /// Satellites per plane (T / P).
    #[must_use]
    pub const fn sats_per_plane(&self) -> u32 {
        if self.num_planes == 0 {
            return 0;
        }
        self.total_satellites / self.num_planes
    }

    /// Generate all satellite orbital elements.
    #[must_use]
    pub fn generate(&self) -> Vec<WalkerSatellite> {
        if self.num_planes == 0 || self.total_satellites == 0 {
            return Vec::new();
        }

        let spp = self.sats_per_plane();
        let two_pi = 2.0 * PI;
        let delta_raan = two_pi / self.num_planes as f64;
        let delta_anomaly = two_pi / spp as f64;
        let phase_offset = if self.total_satellites > 0 {
            two_pi * self.phasing_factor as f64 / self.total_satellites as f64
        } else {
            0.0
        };

        let mut satellites = Vec::with_capacity(self.total_satellites as usize);
        let mut global_idx = 0u32;

        for p in 0..self.num_planes {
            let raan = delta_raan * p as f64;

            for s in 0..spp {
                let base_anomaly = delta_anomaly * s as f64;
                let phased_anomaly = base_anomaly + phase_offset * p as f64;
                // Normalize to [0, 2π)
                let true_anomaly = if two_pi > 0.0 {
                    phased_anomaly % two_pi
                } else {
                    0.0
                };

                let elements = OrbitalElements {
                    semi_major_axis_km: self.semi_major_axis_km,
                    eccentricity: self.eccentricity,
                    inclination_rad: self.inclination_rad,
                    raan_rad: raan,
                    arg_periapsis_rad: 0.0,
                    true_anomaly_rad: true_anomaly,
                };

                // Content hash
                let mut buf = [0u8; 20];
                buf[..4].copy_from_slice(&p.to_le_bytes());
                buf[4..8].copy_from_slice(&s.to_le_bytes());
                buf[8..16].copy_from_slice(&raan.to_bits().to_le_bytes());
                buf[16..20].copy_from_slice(&global_idx.to_le_bytes());
                let content_hash = fnv1a(&buf);

                satellites.push(WalkerSatellite {
                    plane_index: p,
                    sat_index: s,
                    global_index: global_idx,
                    elements,
                    content_hash,
                });

                global_idx += 1;
            }
        }

        satellites
    }

    /// Compute the angular spacing between adjacent orbital planes (rad).
    #[must_use]
    pub fn plane_spacing_rad(&self) -> f64 {
        if self.num_planes == 0 {
            return 0.0;
        }
        2.0 * PI / self.num_planes as f64
    }

    /// Compute the in-plane angular spacing between satellites (rad).
    #[must_use]
    pub fn in_plane_spacing_rad(&self) -> f64 {
        let spp = self.sats_per_plane();
        if spp == 0 {
            return 0.0;
        }
        2.0 * PI / spp as f64
    }

    /// Approximate ground track repeat period in seconds.
    /// Assumes circular orbit with `μ_Earth`.
    #[must_use]
    pub fn ground_track_period_s(&self, mu: f64) -> f64 {
        let a = self.semi_major_axis_km;
        2.0 * PI * (a * a * a / mu).sqrt()
    }

    /// Minimum elevation angle (deg) for single-coverage at equator
    /// (approximate). Based on half-cone angle of the satellite footprint.
    ///
    /// Returns 0.0 if geometry doesn't allow coverage computation.
    #[must_use]
    pub fn min_elevation_equator_deg(&self) -> f64 {
        if self.num_planes == 0 || self.total_satellites == 0 {
            return 0.0;
        }
        // Earth radius
        let r_earth = 6371.0;
        let altitude = self.semi_major_axis_km - r_earth;
        if altitude <= 0.0 {
            return 0.0;
        }

        // Half-cone angle for single coverage: approximate as
        // θ = asin(R_earth / (R_earth + h)) → minimum elevation = 90 - θ - half_spacing
        let sin_theta = r_earth / (r_earth + altitude);
        let theta = sin_theta.asin();
        let half_spacing = PI / self.total_satellites as f64;
        let el_rad = (PI / 2.0 - theta - half_spacing).max(0.0);
        el_rad.to_degrees()
    }
}

// ── Tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    const MU_EARTH: f64 = 398_600.441_8;

    #[test]
    fn gps_constellation_24_6_1() {
        // GPS: 24/6/1, MEO, ~55° inclination
        let walker = WalkerConstellation::new(
            24,
            6,
            1,
            26559.7,               // GPS semi-major axis
            55.0_f64.to_radians(), // 55° inclination
        );
        let sats = walker.generate();
        assert_eq!(sats.len(), 24);
        assert_eq!(walker.sats_per_plane(), 4);
    }

    #[test]
    fn correct_raan_spacing() {
        let walker = WalkerConstellation::new(24, 6, 0, 26559.7, 0.0);
        let sats = walker.generate();
        // RAAN should be 0, 60, 120, 180, 240, 300 degrees
        let expected_raan = 2.0 * PI / 6.0;
        let plane0 = &sats[0];
        let plane1 = &sats[4]; // first sat of plane 1
        let raan_diff = plane1.elements.raan_rad - plane0.elements.raan_rad;
        assert!((raan_diff - expected_raan).abs() < 1e-10);
    }

    #[test]
    fn in_plane_spacing() {
        let walker = WalkerConstellation::new(24, 6, 0, 26559.7, 0.0);
        let sats = walker.generate();
        // 4 sats per plane → spacing = 90°
        let s0 = &sats[0];
        let s1 = &sats[1];
        let anomaly_diff = s1.elements.true_anomaly_rad - s0.elements.true_anomaly_rad;
        assert!((anomaly_diff - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn phasing_factor_offsets_anomaly() {
        let w0 = WalkerConstellation::new(24, 6, 0, 26559.7, 0.0);
        let w1 = WalkerConstellation::new(24, 6, 1, 26559.7, 0.0);
        let sats0 = w0.generate();
        let sats1 = w1.generate();
        // Plane 1, sat 0 should have different true anomaly
        assert!(
            (sats0[4].elements.true_anomaly_rad - sats1[4].elements.true_anomaly_rad).abs() > 1e-10
        );
    }

    #[test]
    fn iridium_constellation() {
        // Iridium: 66/6/2
        let walker = WalkerConstellation::new(66, 6, 2, 7159.0, 86.4_f64.to_radians());
        let sats = walker.generate();
        assert_eq!(sats.len(), 66);
        assert_eq!(walker.sats_per_plane(), 11);
    }

    #[test]
    fn all_elements_share_sma_and_inc() {
        let inc = 55.0_f64.to_radians();
        let sma = 26559.7;
        let walker = WalkerConstellation::new(24, 6, 1, sma, inc);
        let sats = walker.generate();
        for s in &sats {
            assert!((s.elements.semi_major_axis_km - sma).abs() < 1e-10);
            assert!((s.elements.inclination_rad - inc).abs() < 1e-10);
            assert!((s.elements.eccentricity).abs() < 1e-15);
        }
    }

    #[test]
    fn global_indices_unique() {
        let walker = WalkerConstellation::new(24, 6, 1, 26559.7, 0.0);
        let sats = walker.generate();
        let mut indices: Vec<u32> = sats.iter().map(|s| s.global_index).collect();
        indices.sort_unstable();
        indices.dedup();
        assert_eq!(indices.len(), 24);
    }

    #[test]
    fn content_hash_deterministic() {
        let walker = WalkerConstellation::new(24, 6, 1, 26559.7, 0.0);
        let sats1 = walker.generate();
        let sats2 = walker.generate();
        for (a, b) in sats1.iter().zip(sats2.iter()) {
            assert_eq!(a.content_hash, b.content_hash);
            assert_ne!(a.content_hash, 0);
        }
    }

    #[test]
    fn empty_constellation() {
        let walker = WalkerConstellation::new(0, 0, 0, 26559.7, 0.0);
        let sats = walker.generate();
        assert!(sats.is_empty());
    }

    #[test]
    fn plane_spacing() {
        let walker = WalkerConstellation::new(24, 6, 0, 26559.7, 0.0);
        assert!((walker.plane_spacing_rad() - PI / 3.0).abs() < 1e-10);
        assert!((walker.in_plane_spacing_rad() - PI / 2.0).abs() < 1e-10);
    }

    #[test]
    fn ground_track_period_geo() {
        // GEO: a = 42164 km → T ≈ 86164 s (sidereal day)
        let walker = WalkerConstellation::new(3, 1, 0, 42164.0, 0.0);
        let period = walker.ground_track_period_s(MU_EARTH);
        assert!((period - 86164.0).abs() < 100.0);
    }

    #[test]
    fn min_elevation_leo() {
        // LEO constellation at 780 km altitude
        let walker = WalkerConstellation::new(66, 6, 2, 7159.0, 86.4_f64.to_radians());
        let el = walker.min_elevation_equator_deg();
        // Should be some positive value
        assert!(el > 0.0);
        assert!(el < 90.0);
    }

    #[test]
    fn phasing_wraps_modulo_planes() {
        // Phasing factor >= num_planes should wrap
        let w = WalkerConstellation::new(24, 6, 8, 26559.7, 0.0);
        assert_eq!(w.phasing_factor, 8 % 6); // should be 2
    }

    #[test]
    fn single_plane_constellation() {
        let walker = WalkerConstellation::new(4, 1, 0, 7000.0, 0.5);
        let sats = walker.generate();
        assert_eq!(sats.len(), 4);
        // All in plane 0
        for s in &sats {
            assert_eq!(s.plane_index, 0);
            assert!(s.elements.raan_rad.abs() < 1e-10);
        }
    }

    #[test]
    fn min_elevation_below_surface_returns_zero() {
        // SMA below Earth radius → altitude < 0 → should return 0.0
        let walker = WalkerConstellation::new(10, 2, 0, 5000.0, 0.5);
        let el = walker.min_elevation_equator_deg();
        assert!((el).abs() < 1e-10, "Below surface should give 0, got {el}");
    }

    #[test]
    fn zero_planes_sats_per_plane() {
        let walker = WalkerConstellation::new(24, 0, 0, 7000.0, 0.5);
        assert_eq!(walker.sats_per_plane(), 0);
        assert!((walker.plane_spacing_rad()).abs() < 1e-10);
        assert!((walker.in_plane_spacing_rad()).abs() < 1e-10);
    }

    #[test]
    fn starlink_like_constellation() {
        // Starlink第1シェル的な構成: 72面×22衛星 = 1584
        let walker = WalkerConstellation::new(1584, 72, 1, 6921.0, 53.0_f64.to_radians());
        let sats = walker.generate();
        assert_eq!(sats.len(), 1584);
        assert_eq!(walker.sats_per_plane(), 22);
    }

    #[test]
    fn true_anomaly_in_range() {
        // 全衛星のtrue_anomalyが[0, 2π)の範囲内
        let walker = WalkerConstellation::new(24, 6, 3, 26559.7, 0.5);
        let sats = walker.generate();
        for s in &sats {
            assert!(
                s.elements.true_anomaly_rad >= 0.0,
                "anomaly={} < 0",
                s.elements.true_anomaly_rad
            );
            assert!(
                s.elements.true_anomaly_rad < 2.0f64.mul_add(PI, 1e-10),
                "anomaly={} >= 2π",
                s.elements.true_anomaly_rad
            );
        }
    }

    #[test]
    fn ground_track_period_leo() {
        // LEO: a=7000km → T ≈ 5800s
        let walker = WalkerConstellation::new(1, 1, 0, 7000.0, 0.0);
        let period = walker.ground_track_period_s(MU_EARTH);
        assert!((period - 5834.0).abs() < 50.0, "LEO period: {period}");
    }

    #[test]
    fn content_hash_unique_per_satellite() {
        // 各衛星のcontent_hashがユニーク
        let walker = WalkerConstellation::new(24, 6, 1, 26559.7, 0.5);
        let sats = walker.generate();
        let mut hashes: Vec<u64> = sats.iter().map(|s| s.content_hash).collect();
        hashes.sort_unstable();
        hashes.dedup();
        assert_eq!(hashes.len(), 24, "ハッシュの重複あり");
    }

    #[test]
    fn plane_index_matches_raan() {
        // plane_indexが異なる衛星はRAANも異なる
        let walker = WalkerConstellation::new(12, 3, 0, 7000.0, 0.5);
        let sats = walker.generate();
        let plane0_raan = sats[0].elements.raan_rad;
        let plane1_raan = sats[4].elements.raan_rad;
        assert!(
            (plane1_raan - plane0_raan - 2.0 * PI / 3.0).abs() < 1e-10,
            "RAAN差: {}",
            plane1_raan - plane0_raan
        );
    }

    mod prop {
        use super::*;
        use proptest::prelude::*;

        proptest! {
            #[test]
            fn generate_count_equals_total(
                planes in 1u32..10,
                spp in 1u32..20,
            ) {
                let total = planes * spp;
                let w = WalkerConstellation::new(total, planes, 0, 7000.0, 0.5);
                let sats = w.generate();
                prop_assert_eq!(sats.len(), total as usize);
            }

            #[test]
            fn global_indices_sequential(
                planes in 1u32..8,
                spp in 1u32..12,
            ) {
                let total = planes * spp;
                let w = WalkerConstellation::new(total, planes, 0, 7000.0, 0.5);
                let sats = w.generate();
                for (i, s) in sats.iter().enumerate() {
                    prop_assert_eq!(s.global_index, i as u32);
                }
            }

            #[test]
            fn all_sats_share_sma_and_ecc(
                planes in 1u32..6,
                spp in 1u32..10,
                sma in 6500.0f64..50000.0,
            ) {
                let total = planes * spp;
                let w = WalkerConstellation::new(total, planes, 0, sma, 0.5);
                let sats = w.generate();
                for s in &sats {
                    prop_assert!((s.elements.semi_major_axis_km - sma).abs() < 1e-10);
                    prop_assert!(s.elements.eccentricity.abs() < 1e-15);
                }
            }
        }
    }
}
