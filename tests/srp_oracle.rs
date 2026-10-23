//! Oracle: solar radiation pressure (`srp::solar_radiation_pressure`).
//!
//! - In full sunlight the acceleration is (S/c)·Cr·(A/m)·(AU/d)², directed
//!   away from the Sun, with the nominal total solar irradiance S = 1361 W/m²
//!   (IAU 2015 Resolution B3) and c = 299,792,458 m/s: 4.5398 × 10⁻⁶ m/s² for
//!   Cr = 1, A/m = 1 m²/kg at 1 au. Computed here from those numbers.
//! - Inverse square in the distance from the Sun.
//! - In a shadow the pressure is scaled by the visible fraction of the solar
//!   disc from `eclipse::conical_illumination` (zero in the umbra, the same
//!   fraction in the penumbra), and by the product of the fractions for two
//!   occulting bodies.
//!
//! Author: Moroya Sakamoto

use alice_det_math::sqrt64;
use alice_space::eclipse::{
    conical_illumination, EARTH_EQUATORIAL_RADIUS_KM, MOON_RADIUS_KM, SUN_RADIUS_KM,
};
use alice_space::srp::{solar_radiation_pressure, Occulter, SrpError};

const AU_KM: f64 = 149_597_870.7;
const EPS: f64 = f64::EPSILON;

fn norm(v: [f64; 3]) -> f64 {
    sqrt64(v[0] * v[0] + v[1] * v[1] + v[2] * v[2])
}

fn earth() -> Occulter {
    Occulter {
        position: [0.0; 3],
        radius: EARTH_EQUATORIAL_RADIUS_KM,
    }
}

/// (S/c)·Cr·(A/m)·(AU/d)² in km/s².
fn full(cr: f64, am: f64, d_km: f64) -> f64 {
    1361.0 / 299_792_458.0 * cr * am * (AU_KM / d_km) * (AU_KM / d_km) / 1000.0
}

#[test]
fn full_sunlight_matches_the_irradiance_over_c() {
    let sun = [AU_KM, 0.0, 0.0];
    let mut compared = 0;
    for (sat, cr, am) in [
        ([7_000.0, 0.0, 0.0], 1.0, 1.0),
        ([0.0, 42_164.0, 0.0], 1.3, 0.02),
        ([-1.0e6, 3.0e5, 2.0e5], 1.8, 0.01),
        ([1.5e6, 0.0, 0.0], 1.0, 25.0),
    ] {
        let a = solar_radiation_pressure(sat, sun, &[earth()], cr, am).unwrap();
        let d = [sat[0] - sun[0], sat[1] - sun[1], sat[2] - sun[2]];
        let dn = norm(d);
        let want = full(cr, am, dn);
        assert!(
            (norm(a) - want).abs() <= 8.0 * EPS * want,
            "{sat:?}: {} vs {want}",
            norm(a)
        );
        // directed away from the Sun
        for c in 0..3 {
            assert!(
                (a[c] / norm(a) - d[c] / dn).abs() <= 8.0 * EPS,
                "{sat:?} {c}"
            );
        }
        compared += 1;
    }
    assert_eq!(compared, 4);
    // the reference number at 1 au for Cr = 1, A/m = 1 m²/kg
    let at_au = full(1.0, 1.0, AU_KM) * 1000.0;
    assert!((at_au - 4.5398e-6).abs() < 1e-10, "{at_au}");
}

#[test]
fn pressure_falls_with_the_square_of_the_distance() {
    let sat = [0.0, 0.0, 50_000.0];
    let near = solar_radiation_pressure(sat, [AU_KM, 0.0, 0.0], &[], 1.2, 0.05).unwrap();
    let far = solar_radiation_pressure(sat, [2.0 * AU_KM, 0.0, 0.0], &[], 1.2, 0.05).unwrap();
    let ratio = norm(near) / norm(far);
    // the satellite is off the Sun line, so the distances are not exactly 1 : 2
    let (dn, df) = (
        sqrt64(AU_KM * AU_KM + 50_000.0 * 50_000.0),
        sqrt64(4.0 * AU_KM * AU_KM + 50_000.0 * 50_000.0),
    );
    assert!(
        (ratio - (df / dn) * (df / dn)).abs() <= 16.0 * EPS * ratio,
        "{ratio}"
    );
}

#[test]
fn shadows_scale_the_pressure_by_the_visible_fraction() {
    let sun = [AU_KM, 0.0, 0.0];
    // centre of the umbra in LEO: no pressure at all
    let umbra = solar_radiation_pressure([-7_000.0, 0.0, 0.0], sun, &[earth()], 1.0, 1.0).unwrap();
    assert_eq!(umbra, [0.0; 3]);
    // across the penumbra: |a| / |a without the Earth| is the visible fraction
    let mut partial = 0;
    for k in 0..200 {
        let y = EARTH_EQUATORIAL_RADIUS_KM * (0.9 + 0.001 * f64::from(k));
        let sat = [-40_000.0, y, 0.0];
        let lit = solar_radiation_pressure(sat, sun, &[], 1.0, 1.0).unwrap();
        let shaded = solar_radiation_pressure(sat, sun, &[earth()], 1.0, 1.0).unwrap();
        let f = conical_illumination(sat, sun, [0.0; 3], EARTH_EQUATORIAL_RADIUS_KM)
            .unwrap()
            .fraction;
        for c in 0..3 {
            assert!(
                (shaded[c] - f * lit[c]).abs() <= 4.0 * EPS * norm(lit),
                "y {y} {c}: {} vs {}",
                shaded[c],
                f * lit[c]
            );
        }
        partial += usize::from(f > 0.0 && f < 1.0);
    }
    assert!(partial > 20, "only {partial} penumbral points");
    // two occulting bodies: the product of the fractions, in either order. The
    // Moon at quadrature, the satellite 5600 km behind it and 1740 km off its
    // shadow axis: inside the lunar penumbra (umbra ≈ 1711 km, penumbra ≈
    // 1763 km there) and far from the Earth's shadow
    let moon = Occulter {
        position: [0.0, 384_400.0, 0.0],
        radius: MOON_RADIUS_KM,
    };
    let sat = [-5_600.0, 384_400.0 + 1_740.0, 0.0];
    let lit = solar_radiation_pressure(sat, sun, &[], 1.0, 1.0).unwrap();
    let fe = conical_illumination(sat, sun, [0.0; 3], EARTH_EQUATORIAL_RADIUS_KM)
        .unwrap()
        .fraction;
    let fm = conical_illumination(sat, sun, moon.position, moon.radius)
        .unwrap()
        .fraction;
    assert!(
        fm > 0.0 && fm < 1.0,
        "the Moon only partly covers the Sun: {fm}"
    );
    for order in [[earth(), moon], [moon, earth()]] {
        let both = solar_radiation_pressure(sat, sun, &order, 1.0, 1.0).unwrap();
        for c in 0..3 {
            assert!((both[c] - fe * fm * lit[c]).abs() <= 8.0 * EPS * norm(lit));
        }
    }
    let _ = SUN_RADIUS_KM;
}

/// A satellite inside an occulting body, a non-finite input, or a negative
/// coefficient or area give an error, not a value.
#[test]
fn degenerate_inputs() {
    let sun = [AU_KM, 0.0, 0.0];
    assert_eq!(
        solar_radiation_pressure([100.0, 0.0, 0.0], sun, &[earth()], 1.0, 1.0),
        Err(SrpError::Eclipse(
            alice_space::eclipse::EclipseError::InsideBody
        ))
    );
    assert_eq!(
        solar_radiation_pressure([f64::NAN, 7_000.0, 0.0], sun, &[], 1.0, 1.0),
        Err(SrpError::InvalidInput)
    );
    assert_eq!(
        solar_radiation_pressure([7_000.0, 0.0, 0.0], sun, &[], -1.0, 1.0),
        Err(SrpError::InvalidInput)
    );
    assert_eq!(
        solar_radiation_pressure([7_000.0, 0.0, 0.0], sun, &[], 1.0, -0.1),
        Err(SrpError::InvalidInput)
    );
    assert_eq!(
        solar_radiation_pressure(sun, sun, &[], 1.0, 1.0),
        Err(SrpError::InvalidInput)
    );
}
