//! Eclipses of the Sun seen from a spacecraft.
//!
//! The conical model ([`conical_illumination`], [`earth_illumination`]) takes
//! the positions of the satellite, the Sun and the occulting body in one frame
//! (kilometres) and compares the apparent disks of the Sun (radius `a`) and
//! of the body (radius `b`) at their angular separation `c` (Montenbruck &
//! Gill, *Satellite Orbits*, §3.4.2):
//!
//! - `c ≥ a + b`: no eclipse, the whole disk is visible;
//! - `c ≤ b − a`: total (umbra), nothing is visible;
//! - `c ≤ a − b`: annular, the body is inside the disk, `1 − b²/a²` visible;
//! - otherwise partial (penumbra), `1 − A/(πa²)` visible with `A` the area of
//!   the overlap of the two disks.
//!
//! The disks are treated as flat circles of angular radius (the angles are
//! below 0.27° for the Sun seen from near the Earth) and the body as a sphere
//! without atmosphere. Every transcendental is alice-det-math's.
//!
//! The older functions ([`shadow_state`], [`earth_shadow`],
//! [`eclipse_fraction`]) are a cylinder approximation and are deprecated.

use alice_det_math::{acos64, asin64, atan2_64, sin64, sqrt64};

// ── Constants ──────────────────────────────────────────────────────────

/// Mean Earth radius used by the deprecated cylinder functions (km).
const EARTH_RADIUS_KM: f64 = 6371.0;

/// Nominal solar radius (IAU 2015 Resolution B3), km.
pub const SUN_RADIUS_KM: f64 = 695_700.0;

/// Earth equatorial radius (WGS-84), km.
pub const EARTH_EQUATORIAL_RADIUS_KM: f64 = 6378.137;

/// Mean lunar radius (IAU), km.
pub const MOON_RADIUS_KM: f64 = 1737.4;

/// Mean Earth–Sun distance used by the deprecated cylinder functions (km).
const EARTH_SUN_DISTANCE_KM: f64 = 149_597_870.7;

// ── Conical model ──────────────────────────────────────────────────────

/// What the occulting body hides of the Sun, seen from the satellite.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum Eclipse {
    /// The whole disk of the Sun is visible.
    None,
    /// The body covers part of the disk (penumbra).
    Partial,
    /// The body lies inside the disk, a ring of the Sun is visible.
    Annular,
    /// The body covers the whole disk (umbra).
    Total,
}

/// The illumination of a satellite by the Sun past one occulting body.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Illumination {
    /// The kind of eclipse.
    pub eclipse: Eclipse,
    /// Visible fraction of the Sun's disk, 0 (total) to 1 (none).
    pub fraction: f64,
    /// Apparent radius of the Sun, radians.
    pub sun_angular_radius: f64,
    /// Apparent radius of the body, radians.
    pub body_angular_radius: f64,
    /// Angle between the directions to the Sun and to the body, radians.
    pub separation: f64,
}

/// Why the illumination cannot be computed.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub enum EclipseError {
    /// A coordinate or the radius is not finite, or the radius is not positive.
    NonFinite,
    /// The satellite is inside (or on) the occulting body.
    InsideBody,
    /// The satellite is inside (or on) the Sun.
    InsideSun,
}

impl core::fmt::Display for EclipseError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.write_str(match self {
            Self::NonFinite => "non-finite position or non-positive radius",
            Self::InsideBody => "the satellite is inside the occulting body",
            Self::InsideSun => "the satellite is inside the Sun",
        })
    }
}

impl std::error::Error for EclipseError {}

fn sub(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] - b[0], a[1] - b[1], a[2] - b[2]]
}

fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0] * b[0] + a[1] * b[1] + a[2] * b[2]
}

fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}

fn norm(a: [f64; 3]) -> f64 {
    sqrt64(dot(a, a))
}

/// The illumination of a satellite at `sat` by the Sun at `sun` (radius
/// [`SUN_RADIUS_KM`]) past a spherical body of radius `body_radius` at
/// `body`, all positions in one frame, in kilometres.
pub fn conical_illumination(
    sat: [f64; 3],
    sun: [f64; 3],
    body: [f64; 3],
    body_radius: f64,
) -> Result<Illumination, EclipseError> {
    let finite = sat.iter().chain(&sun).chain(&body).all(|v| v.is_finite());
    if !finite || !body_radius.is_finite() || body_radius <= 0.0 {
        return Err(EclipseError::NonFinite);
    }
    let to_sun = sub(sun, sat);
    let to_body = sub(body, sat);
    let (d_sun, d_body) = (norm(to_sun), norm(to_body));
    if d_body <= body_radius {
        return Err(EclipseError::InsideBody);
    }
    if d_sun <= SUN_RADIUS_KM {
        return Err(EclipseError::InsideSun);
    }
    let a = asin64(SUN_RADIUS_KM / d_sun);
    let b = asin64(body_radius / d_body);
    // atan2 of |u × v| and u · v keeps the digits of small separations
    let c = atan2_64(norm(cross(to_sun, to_body)), dot(to_sun, to_body));
    let (eclipse, fraction) = if c >= a + b {
        (Eclipse::None, 1.0)
    } else if c <= b - a {
        (Eclipse::Total, 0.0)
    } else if c <= a - b {
        (Eclipse::Annular, 1.0 - (b * b) / (a * a))
    } else {
        // the overlap of disks a and b at distance c is two circular segments
        // cut by the common chord, at xa from the centre of the Sun's disk and
        // xb from the body's. The textbook form a²·θa + b²·θb − c·y subtracts
        // the triangle c·y from the sectors and, near tangency, leaves a tiny
        // difference of large terms (a visible fraction above 1 by 6e-6 was
        // measured); each segment a²·(θ − sin θ cos θ) is formed on its own
        // instead, and the visible part from the segments on its own side, so
        // neither end (fraction near 1 or near 0) cancels
        // c² − b² and c² − a² as products: on the edge of a body much larger
        // than the Sun's disk c ≈ b, and the squared form would lose the digits
        // of c − b
        let xa = (a * a + (c - b) * (c + b)) / (2.0 * c);
        let xb = (b * b + (c - a) * (c + a)) / (2.0 * c);
        let (ta, tb) = (acos64(xa / a), acos64(xb / b));
        let disk = core::f64::consts::PI * a * a;
        let fraction = if ta <= core::f64::consts::FRAC_PI_2 {
            1.0 - (a * a * segment(ta) + b * b * segment(tb)) / disk
        } else {
            (a * a * segment(core::f64::consts::PI - ta) - b * b * segment(tb)) / disk
        };
        (Eclipse::Partial, fraction)
    };
    Ok(Illumination {
        eclipse,
        fraction,
        sun_angular_radius: a,
        body_angular_radius: b,
        separation: c,
    })
}

/// Area of the circular segment of the unit circle with half-angle `t`
/// (`0 ≤ t ≤ π`): `t − sin t · cos t = (u − sin u) / 2` with `u = 2t`. Its
/// absolute error is a few ulp of `u`, which is what the visible fraction
/// sees; a series for small `u` would only sharpen the relative error of a
/// segment far below the conditioning of the fraction on the angles.
fn segment(t: f64) -> f64 {
    let u = 2.0 * t;
    (u - sin64(u)) / 2.0
}

/// [`conical_illumination`] with the Earth (equatorial radius
/// [`EARTH_EQUATORIAL_RADIUS_KM`]) at the origin: `sat` and `sun` are
/// geocentric, in kilometres.
pub fn earth_illumination(sat: [f64; 3], sun: [f64; 3]) -> Result<Illumination, EclipseError> {
    conical_illumination(sat, sun, [0.0; 3], EARTH_EQUATORIAL_RADIUS_KM)
}

// ── Cylinder approximation (deprecated) ────────────────────────────────

/// 日食/影の種別 (cylinder approximation, see [`shadow_state`]).
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShadowState {
    /// 完全に太陽光が当たっている。
    Sunlit,
    /// 本影 (cylinder approximation)。
    Umbra,
    /// 半影 (fixed-factor approximation)。
    Penumbra,
}

/// Shadow by a cylinder approximation, **not a cone**: the umbra is the
/// cylinder of radius `body_radius` behind the body along `sun_dir`, at any
/// distance, and the penumbra a cylinder wider by the fixed factor
/// `1 + R_sun / 1 au`. It ignores that the umbra narrows and the penumbra
/// widens with distance, and has no annular case.
///
/// - `sat_pos`: satellite position (km), centred on the body
/// - `sun_dir`: unit vector from the body towards the Sun
/// - `body_radius`: radius of the occulting body (km)
#[deprecated(note = "cylinder approximation; use conical_illumination (removed in the next minor)")]
#[must_use]
pub fn shadow_state(
    sat_pos: (f64, f64, f64),
    sun_dir: (f64, f64, f64),
    body_radius: f64,
) -> ShadowState {
    let dot = sat_pos.0 * sun_dir.0 + sat_pos.1 * sun_dir.1 + sat_pos.2 * sun_dir.2;
    // on the Sun's side of the body: sunlit
    if dot >= 0.0 {
        return ShadowState::Sunlit;
    }
    // distance from the axis through the body along the Sun direction
    let perp_x = sat_pos.0 - dot * sun_dir.0;
    let perp_y = sat_pos.1 - dot * sun_dir.1;
    let perp_z = sat_pos.2 - dot * sun_dir.2;
    let perp_dist_sq = perp_x * perp_x + perp_y * perp_y + perp_z * perp_z;
    let umbra_radius = body_radius;
    let penumbra_radius = body_radius * (1.0 + SUN_RADIUS_KM / EARTH_SUN_DISTANCE_KM);
    if perp_dist_sq < umbra_radius * umbra_radius {
        ShadowState::Umbra
    } else if perp_dist_sq < penumbra_radius * penumbra_radius {
        ShadowState::Penumbra
    } else {
        ShadowState::Sunlit
    }
}

/// [`shadow_state`] for the Earth (mean radius 6371 km), the same cylinder
/// approximation, **not a cone**.
#[deprecated(note = "cylinder approximation; use earth_illumination (removed in the next minor)")]
#[allow(deprecated)]
#[must_use]
pub fn earth_shadow(sat_pos: (f64, f64, f64), sun_dir: (f64, f64, f64)) -> ShadowState {
    shadow_state(sat_pos, sun_dir, EARTH_RADIUS_KM)
}

/// Fraction of a circular orbit inside a cylindrical Earth shadow (mean
/// radius 6371 km), worst case of the Sun in the orbit plane:
/// `asin(R/r)/π`. A cylinder approximation, **not a cone**, without
/// penumbra.
///
/// - `altitude_km`: orbit altitude (km)
///
/// Returns 0.0 (always sunlit) to 1.0.
#[deprecated(
    note = "cylinder approximation; sample conical_illumination along the orbit (removed in the next minor)"
)]
#[must_use]
pub fn eclipse_fraction(altitude_km: f64) -> f64 {
    let r = EARTH_RADIUS_KM + altitude_km;
    if r <= EARTH_RADIUS_KM {
        return 1.0;
    }
    let sin_half = EARTH_RADIUS_KM / r;
    if sin_half >= 1.0 {
        return 1.0;
    }
    asin64(sin_half) / core::f64::consts::PI
}

/// 通信ウィンドウの判定: 衛星が日照中かどうか。
///
/// 多くの深宇宙通信は太陽光パネルに依存するため、
/// 日食中はパワーが制限され通信不可となる場合がある。
#[must_use]
pub const fn is_comm_window(state: ShadowState) -> bool {
    matches!(state, ShadowState::Sunlit | ShadowState::Penumbra)
}

// ── Tests ──────────────────────────────────────────────────────────────

#[cfg(test)]
#[allow(deprecated)]
mod tests {
    use super::*;

    #[test]
    fn sunlit_when_facing_sun() {
        // 衛星が太陽方向にいる → 日照
        let sat = (7000.0, 0.0, 0.0);
        let sun_dir = (1.0, 0.0, 0.0);
        assert_eq!(earth_shadow(sat, sun_dir), ShadowState::Sunlit);
    }

    #[test]
    fn umbra_behind_earth() {
        // 衛星が地球の反対側（太陽の真裏）
        let sat = (-7000.0, 0.0, 0.0);
        let sun_dir = (1.0, 0.0, 0.0);
        assert_eq!(earth_shadow(sat, sun_dir), ShadowState::Umbra);
    }

    #[test]
    fn sunlit_when_offset_from_shadow() {
        // 衛星が影の外側にオフセット
        let sat = (-7000.0, 10000.0, 0.0);
        let sun_dir = (1.0, 0.0, 0.0);
        assert_eq!(earth_shadow(sat, sun_dir), ShadowState::Sunlit);
    }

    #[test]
    fn eclipse_fraction_leo() {
        // ISS高度 (~400km) の日食割合: 約35%程度
        let frac = eclipse_fraction(400.0);
        assert!(frac > 0.2);
        assert!(frac < 0.5);
    }

    #[test]
    fn eclipse_fraction_geo() {
        // GEO高度 (~35786km) の日食割合: 非常に小さい
        let frac = eclipse_fraction(35_786.0);
        assert!(frac < 0.05);
    }

    #[test]
    fn eclipse_fraction_zero_altitude() {
        // 地表: 割合 1.0
        let frac = eclipse_fraction(0.0);
        assert!((frac - 1.0).abs() < 1e-10);
    }

    #[test]
    fn eclipse_fraction_very_high() {
        // 非常に高い軌道: 割合はほぼ0
        let frac = eclipse_fraction(100_000.0);
        assert!(frac < 0.02);
    }

    #[test]
    fn comm_window_sunlit() {
        assert!(is_comm_window(ShadowState::Sunlit));
    }

    #[test]
    fn comm_window_penumbra() {
        assert!(is_comm_window(ShadowState::Penumbra));
    }

    #[test]
    fn comm_window_umbra_blocked() {
        assert!(!is_comm_window(ShadowState::Umbra));
    }

    #[test]
    fn shadow_state_equality() {
        assert_eq!(ShadowState::Sunlit, ShadowState::Sunlit);
        assert_ne!(ShadowState::Sunlit, ShadowState::Umbra);
        assert_ne!(ShadowState::Umbra, ShadowState::Penumbra);
    }

    #[test]
    fn custom_body_radius() {
        // 月の半径で影を計算
        let moon_radius = 1737.4;
        let sat = (-3000.0, 0.0, 0.0);
        let sun_dir = (1.0, 0.0, 0.0);
        assert_eq!(shadow_state(sat, sun_dir, moon_radius), ShadowState::Umbra);
    }

    #[test]
    fn custom_body_large_offset_sunlit() {
        let sat = (-3000.0, 5000.0, 0.0);
        let sun_dir = (1.0, 0.0, 0.0);
        assert_eq!(shadow_state(sat, sun_dir, 1737.4), ShadowState::Sunlit);
    }

    #[test]
    fn shadow_state_behind_along_y() {
        // Y軸方向の太陽で-Y方向の衛星 → 影
        let sat = (0.0, -7000.0, 0.0);
        let sun_dir = (0.0, 1.0, 0.0);
        assert_eq!(earth_shadow(sat, sun_dir), ShadowState::Umbra);
    }

    #[test]
    fn shadow_state_behind_along_z() {
        // Z軸方向の太陽で-Z方向の衛星 → 影
        let sat = (0.0, 0.0, -7000.0);
        let sun_dir = (0.0, 0.0, 1.0);
        assert_eq!(earth_shadow(sat, sun_dir), ShadowState::Umbra);
    }

    #[test]
    fn shadow_perpendicular_to_sun_sunlit() {
        // 太陽方向に対して垂直 → dot=0 → 日照
        let sat = (0.0, 7000.0, 0.0);
        let sun_dir = (1.0, 0.0, 0.0);
        assert_eq!(earth_shadow(sat, sun_dir), ShadowState::Sunlit);
    }

    #[test]
    fn eclipse_fraction_negative_altitude() {
        // 負の高度（地下） → 割合1.0
        let frac = eclipse_fraction(-100.0);
        assert!((frac - 1.0).abs() < 1e-10);
    }

    #[test]
    fn eclipse_fraction_between_0_and_1() {
        // 任意の正高度で0 < frac < 1
        for alt in [100.0, 500.0, 1000.0, 5000.0, 20000.0] {
            let frac = eclipse_fraction(alt);
            assert!(frac > 0.0 && frac < 1.0, "altitude={alt}, frac={frac}");
        }
    }

    #[test]
    fn penumbra_detection() {
        // 本影の境界外、半影の境界内の位置を作成
        // body_radius=6371, penumbra_radius ≈ 6371 * 1.00465
        let penumbra_factor = 1.0 + SUN_RADIUS_KM / EARTH_SUN_DISTANCE_KM;
        let penumbra_r = EARTH_RADIUS_KM * penumbra_factor;
        let mid = f64::midpoint(EARTH_RADIUS_KM, penumbra_r);
        let sat = (-10000.0, mid, 0.0);
        let sun_dir = (1.0, 0.0, 0.0);
        assert_eq!(
            shadow_state(sat, sun_dir, EARTH_RADIUS_KM),
            ShadowState::Penumbra
        );
    }

    #[test]
    fn eclipse_fraction_monotonically_decreasing() {
        let f400 = eclipse_fraction(400.0);
        let f1000 = eclipse_fraction(1000.0);
        let f5000 = eclipse_fraction(5000.0);
        assert!(f400 > f1000);
        assert!(f1000 > f5000);
    }
}
