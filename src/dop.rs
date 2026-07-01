//! Dilution of Precision (`DOP`) computation for `GNSS` receivers.
//!
//! Given the position of a receiver and the positions of a set of visible
//! satellites in the same Cartesian frame (`ECEF` in meters is the common
//! choice), this module computes:
//!
//! - `GDOP` — Geometric Dilution of Precision
//! - `PDOP` — Position Dilution of Precision
//! - `HDOP` — Horizontal Dilution of Precision (local ENU frame)
//! - `VDOP` — Vertical Dilution of Precision (local ENU frame)
//! - `TDOP` — Time Dilution of Precision
//!
//! The algorithm builds the geometry matrix `G` whose i-th row is
//! `[-ux_i, -uy_i, -uz_i, 1]`, where `u_i` is the unit line-of-sight vector
//! from the receiver to satellite i. `Q = (Gᵀ G)⁻¹` is a 4x4 covariance-shaped
//! matrix; the DOP metrics are square roots of trace subsets.
//!
//! For `HDOP` / `VDOP` the horizontal covariance block is first rotated into
//! the local East-North-Up frame around the receiver's geodetic latitude and
//! longitude.

use crate::geodetic::{ecef_to_geodetic, Ecef};

// ---------------------------------------------------------------------------
// Output type
// ---------------------------------------------------------------------------

/// A bundle of Dilution-of-Precision metrics.
///
/// Smaller numbers indicate better geometry.  A `PDOP` of 1.0 is ideal;
/// `PDOP > 6` is generally considered poor.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DopValues {
    pub gdop: f64,
    pub pdop: f64,
    pub hdop: f64,
    pub vdop: f64,
    pub tdop: f64,
}

// ---------------------------------------------------------------------------
// Public API
// ---------------------------------------------------------------------------

/// Compute all `DOP` metrics for the receiver observing `satellites`.
///
/// Returns `None` when fewer than four satellites are supplied or when the
/// geometry matrix is (numerically) singular.
#[must_use]
pub fn compute_dop(receiver: Ecef, satellites: &[Ecef]) -> Option<DopValues> {
    if satellites.len() < 4 {
        return None;
    }

    // Geometry matrix G (N x 4).
    let mut g_rows: Vec<[f64; 4]> = Vec::with_capacity(satellites.len());
    for sat in satellites {
        let dx = sat.x_m - receiver.x_m;
        let dy = sat.y_m - receiver.y_m;
        let dz = sat.z_m - receiver.z_m;
        let range = (dx * dx + dy * dy + dz * dz).sqrt();
        if range < 1.0 {
            // Coincident satellite/receiver.
            return None;
        }
        g_rows.push([-dx / range, -dy / range, -dz / range, 1.0]);
    }

    // Compute Gᵀ G (4 x 4 symmetric).
    let mut ata = [[0.0_f64; 4]; 4];
    for row in &g_rows {
        for i in 0..4 {
            for j in 0..4 {
                ata[i][j] += row[i] * row[j];
            }
        }
    }

    // Invert.
    let q = invert_4x4(ata)?;

    // Diagonal-based DOP metrics in ECEF.
    let sigma_x2 = q[0][0];
    let sigma_y2 = q[1][1];
    let sigma_z2 = q[2][2];
    let sigma_t2 = q[3][3];

    if sigma_x2 < 0.0 || sigma_y2 < 0.0 || sigma_z2 < 0.0 || sigma_t2 < 0.0 {
        return None;
    }

    let pdop = (sigma_x2 + sigma_y2 + sigma_z2).sqrt();
    let tdop = sigma_t2.sqrt();
    let gdop = (sigma_x2 + sigma_y2 + sigma_z2 + sigma_t2).sqrt();

    // Rotate the 3x3 position block of Q into East-North-Up around the
    // receiver's geodetic latitude and longitude.
    let g_pos = ecef_to_geodetic(receiver);
    let (sin_lat, cos_lat) = g_pos.lat_rad.sin_cos();
    let (sin_lon, cos_lon) = g_pos.lon_rad.sin_cos();

    // Rotation from ECEF to ENU is:
    // R = [[-sin_lon,          cos_lon,          0     ],
    //      [-sin_lat*cos_lon, -sin_lat*sin_lon,  cos_lat],
    //      [ cos_lat*cos_lon,  cos_lat*sin_lon,  sin_lat]]
    let r = [
        [-sin_lon, cos_lon, 0.0],
        [-sin_lat * cos_lon, -sin_lat * sin_lon, cos_lat],
        [cos_lat * cos_lon, cos_lat * sin_lon, sin_lat],
    ];

    // Q_enu = R * Q_pos * Rᵀ, so trace = sum over i of (R Q_pos Rᵀ)_ii.
    // We only need diagonal entries of Q_enu for HDOP/VDOP.
    let q_pos = [
        [q[0][0], q[0][1], q[0][2]],
        [q[1][0], q[1][1], q[1][2]],
        [q[2][0], q[2][1], q[2][2]],
    ];
    let mut q_enu_diag = [0.0_f64; 3];
    for i in 0..3 {
        // Row i of R · Q_pos.
        let mut rq = [0.0_f64; 3];
        for k in 0..3 {
            rq[k] = r[i][0] * q_pos[0][k] + r[i][1] * q_pos[1][k] + r[i][2] * q_pos[2][k];
        }
        // (R Q_pos Rᵀ)_ii = rq · Rᵀ_col_i = rq · row_i_of_R.
        q_enu_diag[i] = rq[0] * r[i][0] + rq[1] * r[i][1] + rq[2] * r[i][2];
    }

    let east2 = q_enu_diag[0].max(0.0);
    let north2 = q_enu_diag[1].max(0.0);
    let up2 = q_enu_diag[2].max(0.0);

    let hdop = (east2 + north2).sqrt();
    let vdop = up2.sqrt();

    Some(DopValues {
        gdop,
        pdop,
        hdop,
        vdop,
        tdop,
    })
}

// ---------------------------------------------------------------------------
// 4x4 matrix inversion (adjugate method)
// ---------------------------------------------------------------------------

/// Invert a 4x4 matrix using the adjugate/determinant approach.  Returns
/// `None` when the determinant is smaller than `1e-14` in absolute value.
#[must_use]
fn invert_4x4(m: [[f64; 4]; 4]) -> Option<[[f64; 4]; 4]> {
    let mut inv = [[0.0_f64; 4]; 4];

    inv[0][0] =
        m[1][1] * m[2][2] * m[3][3] - m[1][1] * m[2][3] * m[3][2] - m[2][1] * m[1][2] * m[3][3]
            + m[2][1] * m[1][3] * m[3][2]
            + m[3][1] * m[1][2] * m[2][3]
            - m[3][1] * m[1][3] * m[2][2];
    inv[0][1] =
        -m[0][1] * m[2][2] * m[3][3] + m[0][1] * m[2][3] * m[3][2] + m[2][1] * m[0][2] * m[3][3]
            - m[2][1] * m[0][3] * m[3][2]
            - m[3][1] * m[0][2] * m[2][3]
            + m[3][1] * m[0][3] * m[2][2];
    inv[0][2] =
        m[0][1] * m[1][2] * m[3][3] - m[0][1] * m[1][3] * m[3][2] - m[1][1] * m[0][2] * m[3][3]
            + m[1][1] * m[0][3] * m[3][2]
            + m[3][1] * m[0][2] * m[1][3]
            - m[3][1] * m[0][3] * m[1][2];
    inv[0][3] =
        -m[0][1] * m[1][2] * m[2][3] + m[0][1] * m[1][3] * m[2][2] + m[1][1] * m[0][2] * m[2][3]
            - m[1][1] * m[0][3] * m[2][2]
            - m[2][1] * m[0][2] * m[1][3]
            + m[2][1] * m[0][3] * m[1][2];

    inv[1][0] =
        -m[1][0] * m[2][2] * m[3][3] + m[1][0] * m[2][3] * m[3][2] + m[2][0] * m[1][2] * m[3][3]
            - m[2][0] * m[1][3] * m[3][2]
            - m[3][0] * m[1][2] * m[2][3]
            + m[3][0] * m[1][3] * m[2][2];
    inv[1][1] =
        m[0][0] * m[2][2] * m[3][3] - m[0][0] * m[2][3] * m[3][2] - m[2][0] * m[0][2] * m[3][3]
            + m[2][0] * m[0][3] * m[3][2]
            + m[3][0] * m[0][2] * m[2][3]
            - m[3][0] * m[0][3] * m[2][2];
    inv[1][2] =
        -m[0][0] * m[1][2] * m[3][3] + m[0][0] * m[1][3] * m[3][2] + m[1][0] * m[0][2] * m[3][3]
            - m[1][0] * m[0][3] * m[3][2]
            - m[3][0] * m[0][2] * m[1][3]
            + m[3][0] * m[0][3] * m[1][2];
    inv[1][3] =
        m[0][0] * m[1][2] * m[2][3] - m[0][0] * m[1][3] * m[2][2] - m[1][0] * m[0][2] * m[2][3]
            + m[1][0] * m[0][3] * m[2][2]
            + m[2][0] * m[0][2] * m[1][3]
            - m[2][0] * m[0][3] * m[1][2];

    inv[2][0] =
        m[1][0] * m[2][1] * m[3][3] - m[1][0] * m[2][3] * m[3][1] - m[2][0] * m[1][1] * m[3][3]
            + m[2][0] * m[1][3] * m[3][1]
            + m[3][0] * m[1][1] * m[2][3]
            - m[3][0] * m[1][3] * m[2][1];
    inv[2][1] =
        -m[0][0] * m[2][1] * m[3][3] + m[0][0] * m[2][3] * m[3][1] + m[2][0] * m[0][1] * m[3][3]
            - m[2][0] * m[0][3] * m[3][1]
            - m[3][0] * m[0][1] * m[2][3]
            + m[3][0] * m[0][3] * m[2][1];
    inv[2][2] =
        m[0][0] * m[1][1] * m[3][3] - m[0][0] * m[1][3] * m[3][1] - m[1][0] * m[0][1] * m[3][3]
            + m[1][0] * m[0][3] * m[3][1]
            + m[3][0] * m[0][1] * m[1][3]
            - m[3][0] * m[0][3] * m[1][1];
    inv[2][3] =
        -m[0][0] * m[1][1] * m[2][3] + m[0][0] * m[1][3] * m[2][1] + m[1][0] * m[0][1] * m[2][3]
            - m[1][0] * m[0][3] * m[2][1]
            - m[2][0] * m[0][1] * m[1][3]
            + m[2][0] * m[0][3] * m[1][1];

    inv[3][0] =
        -m[1][0] * m[2][1] * m[3][2] + m[1][0] * m[2][2] * m[3][1] + m[2][0] * m[1][1] * m[3][2]
            - m[2][0] * m[1][2] * m[3][1]
            - m[3][0] * m[1][1] * m[2][2]
            + m[3][0] * m[1][2] * m[2][1];
    inv[3][1] =
        m[0][0] * m[2][1] * m[3][2] - m[0][0] * m[2][2] * m[3][1] - m[2][0] * m[0][1] * m[3][2]
            + m[2][0] * m[0][2] * m[3][1]
            + m[3][0] * m[0][1] * m[2][2]
            - m[3][0] * m[0][2] * m[2][1];
    inv[3][2] =
        -m[0][0] * m[1][1] * m[3][2] + m[0][0] * m[1][2] * m[3][1] + m[1][0] * m[0][1] * m[3][2]
            - m[1][0] * m[0][2] * m[3][1]
            - m[3][0] * m[0][1] * m[1][2]
            + m[3][0] * m[0][2] * m[1][1];
    inv[3][3] =
        m[0][0] * m[1][1] * m[2][2] - m[0][0] * m[1][2] * m[2][1] - m[1][0] * m[0][1] * m[2][2]
            + m[1][0] * m[0][2] * m[2][1]
            + m[2][0] * m[0][1] * m[1][2]
            - m[2][0] * m[0][2] * m[1][1];

    let det = m[0][0] * inv[0][0] + m[0][1] * inv[1][0] + m[0][2] * inv[2][0] + m[0][3] * inv[3][0];
    if det.abs() < 1e-14 {
        return None;
    }
    let inv_det = 1.0 / det;
    for row in &mut inv {
        for cell in row.iter_mut() {
            *cell *= inv_det;
        }
    }
    Some(inv)
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geodetic::{geodetic_to_ecef, Geodetic};

    fn tokyo() -> Ecef {
        geodetic_to_ecef(Geodetic::from_degrees(35.6895, 139.6917, 0.0))
    }

    /// Place four satellites in an approximately symmetric tetrahedron
    /// around the Tokyo receiver at 20_000 km range.
    fn tetrahedron_sats() -> Vec<Ecef> {
        let r = tokyo();
        let d = 20_000_000.0;
        vec![
            Ecef {
                x_m: r.x_m + d,
                y_m: r.y_m,
                z_m: r.z_m + d,
            },
            Ecef {
                x_m: r.x_m - d,
                y_m: r.y_m + d,
                z_m: r.z_m + d,
            },
            Ecef {
                x_m: r.x_m - d,
                y_m: r.y_m - d,
                z_m: r.z_m + d,
            },
            Ecef {
                x_m: r.x_m + d,
                y_m: r.y_m + d,
                z_m: r.z_m - d,
            },
        ]
    }

    #[test]
    fn returns_none_with_fewer_than_four_sats() {
        let receiver = tokyo();
        let sats = vec![
            Ecef {
                x_m: 20_000_000.0,
                y_m: 0.0,
                z_m: 0.0,
            };
            3
        ];
        assert!(compute_dop(receiver, &sats).is_none());
    }

    #[test]
    fn tetrahedron_produces_finite_dop() {
        let receiver = tokyo();
        let dop = compute_dop(receiver, &tetrahedron_sats()).expect("finite geometry");
        assert!(dop.pdop.is_finite());
        assert!(dop.hdop.is_finite());
        assert!(dop.vdop.is_finite());
        assert!(dop.tdop.is_finite());
        assert!(dop.gdop.is_finite());
        // A tetrahedron should yield reasonable (< 10) DOP.
        assert!(dop.pdop < 10.0, "pdop = {}", dop.pdop);
    }

    #[test]
    fn more_satellites_improve_geometry() {
        // Adding two extra satellites to the tetrahedron cannot make PDOP worse
        // than the 4-satellite baseline (least-squares residual is monotonic
        // in observations for full-rank systems).
        let receiver = tokyo();
        let mut sats = tetrahedron_sats();
        let dop_4 = compute_dop(receiver, &sats).unwrap();
        let d = 20_000_000.0;
        sats.push(Ecef {
            x_m: receiver.x_m,
            y_m: receiver.y_m + d,
            z_m: receiver.z_m + d,
        });
        sats.push(Ecef {
            x_m: receiver.x_m,
            y_m: receiver.y_m - d,
            z_m: receiver.z_m + d,
        });
        let dop_6 = compute_dop(receiver, &sats).unwrap();
        assert!(
            dop_6.pdop <= dop_4.pdop + 1e-9,
            "pdop should not increase: 6-sat={}, 4-sat={}",
            dop_6.pdop,
            dop_4.pdop
        );
    }

    #[test]
    fn gdop_is_composed_of_pdop_and_tdop() {
        let receiver = tokyo();
        let dop = compute_dop(receiver, &tetrahedron_sats()).unwrap();
        let sum = (dop.pdop * dop.pdop + dop.tdop * dop.tdop).sqrt();
        assert!((sum - dop.gdop).abs() < 1e-9);
    }

    #[test]
    fn pdop_is_composed_of_hdop_and_vdop() {
        let receiver = tokyo();
        let dop = compute_dop(receiver, &tetrahedron_sats()).unwrap();
        let sum = (dop.hdop * dop.hdop + dop.vdop * dop.vdop).sqrt();
        assert!((sum - dop.pdop).abs() < 1e-9);
    }

    #[test]
    fn coincident_satellite_returns_none() {
        let receiver = tokyo();
        let sats = vec![receiver; 4];
        assert!(compute_dop(receiver, &sats).is_none());
    }

    #[test]
    fn identity_inverse() {
        let m = [
            [1.0, 0.0, 0.0, 0.0],
            [0.0, 1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0, 0.0],
            [0.0, 0.0, 0.0, 1.0],
        ];
        let inv = invert_4x4(m).unwrap();
        for i in 0..4 {
            for j in 0..4 {
                let expected = if i == j { 1.0 } else { 0.0 };
                assert!((inv[i][j] - expected).abs() < 1e-12);
            }
        }
    }

    #[test]
    fn singular_matrix_returns_none() {
        let singular = [
            [1.0, 2.0, 3.0, 4.0],
            [2.0, 4.0, 6.0, 8.0],
            [3.0, 6.0, 9.0, 12.0],
            [4.0, 8.0, 12.0, 16.0],
        ];
        assert!(invert_4x4(singular).is_none());
    }
}
