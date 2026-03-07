//! Eclipse / shadow detection for spacecraft visibility windows.
//!
//! Determines whether a spacecraft at orbital altitude is in the
//! shadow of a celestial body (eclipse) or illuminated by the sun.
//! Uses cylindrical shadow model for computational efficiency.

// ── Constants ──────────────────────────────────────────────────────────

/// 地球の平均半径 (km)。
const EARTH_RADIUS_KM: f64 = 6371.0;

/// 太陽の平均半径 (km)。
const SUN_RADIUS_KM: f64 = 695_700.0;

/// 地球−太陽の平均距離 (km)。
const EARTH_SUN_DISTANCE_KM: f64 = 149_597_870.7;

// ── Shadow Model ───────────────────────────────────────────────────────

/// 日食/影の種別。
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ShadowState {
    /// 完全に太陽光が当たっている。
    Sunlit,
    /// 本影（完全な影）に入っている。
    Umbra,
    /// 半影（部分的な影）に入っている。
    Penumbra,
}

/// 衛星の位置と太陽方向から影の状態を判定する。
///
/// 円筒影モデル: 遮蔽天体（地球）の影を円錐として近似する。
///
/// # Arguments
///
/// - `sat_pos`: 衛星位置 (x, y, z) [km] — 地球中心座標系
/// - `sun_dir`: 太陽方向の単位ベクトル (x, y, z) — 地球から太陽への方向
/// - `body_radius`: 遮蔽天体の半径 [km]
#[must_use]
pub fn shadow_state(
    sat_pos: (f64, f64, f64),
    sun_dir: (f64, f64, f64),
    body_radius: f64,
) -> ShadowState {
    // 衛星位置の太陽方向成分
    let dot = (sat_pos.0).mul_add(
        sun_dir.0,
        (sat_pos.1).mul_add(sun_dir.1, sat_pos.2 * sun_dir.2),
    );

    // 太陽の反対側にいない場合は日照
    if dot >= 0.0 {
        return ShadowState::Sunlit;
    }

    // 遮蔽天体軸からの距離（太陽方向に対する垂直成分）
    let perp_x = (-dot).mul_add(sun_dir.0, sat_pos.0);
    let perp_y = (-dot).mul_add(sun_dir.1, sat_pos.1);
    let perp_z = (-dot).mul_add(sun_dir.2, sat_pos.2);
    let perp_dist_sq = perp_x.mul_add(perp_x, perp_y.mul_add(perp_y, perp_z * perp_z));

    // 本影半径（円錐モデル）
    let umbra_radius = body_radius;
    // 半影半径（太陽の見かけの大きさを考慮）
    let penumbra_factor = 1.0 + SUN_RADIUS_KM / EARTH_SUN_DISTANCE_KM;
    let penumbra_radius = body_radius * penumbra_factor;

    if perp_dist_sq < umbra_radius * umbra_radius {
        ShadowState::Umbra
    } else if perp_dist_sq < penumbra_radius * penumbra_radius {
        ShadowState::Penumbra
    } else {
        ShadowState::Sunlit
    }
}

/// 地球影モデルのショートカット（`body_radius` = 地球半径）。
#[must_use]
pub fn earth_shadow(sat_pos: (f64, f64, f64), sun_dir: (f64, f64, f64)) -> ShadowState {
    shadow_state(sat_pos, sun_dir, EARTH_RADIUS_KM)
}

/// LEO 衛星の1軌道あたりの日食割合を推定する。
///
/// 円軌道を仮定し、軌道面が太陽方向に対して垂直な最悪ケースを計算する。
///
/// # Arguments
///
/// - `altitude_km`: 軌道高度 [km]
///
/// # Returns
///
/// 日食割合（0.0〜1.0）。0.0 = 常に日照、1.0 = 常に影（実際には起こらない）。
#[must_use]
pub fn eclipse_fraction(altitude_km: f64) -> f64 {
    let r = EARTH_RADIUS_KM + altitude_km;
    if r <= EARTH_RADIUS_KM {
        return 1.0;
    }
    // 最悪ケース: 軌道面が太陽方向に垂直
    // 影角度 = 2 * arcsin(R_earth / r)
    // 日食割合 = 影角度 / (2π)
    let sin_half = EARTH_RADIUS_KM / r;
    if sin_half >= 1.0 {
        return 1.0;
    }
    let half_angle = sin_half.asin();
    half_angle / core::f64::consts::PI
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
        let mid = (EARTH_RADIUS_KM + penumbra_r) / 2.0;
        let sat = (-10000.0, mid, 0.0);
        let sun_dir = (1.0, 0.0, 0.0);
        assert_eq!(shadow_state(sat, sun_dir, EARTH_RADIUS_KM), ShadowState::Penumbra);
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
