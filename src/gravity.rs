//! Zonal gravity: the axisymmetric part of a body's gravity field.
//!
//! The potential is U = (μ/r) (1 − Σₙ Jₙ (R/r)ⁿ Pₙ(z/r)) for n = 2, 3, …,
//! with Pₙ the Legendre polynomials and z along the body's axis of symmetry.
//! The acceleration is its gradient,
//!
//! a = −(μ/r²) [ r̂ (1 − Σ Jₙ ρⁿ ((n+1) Pₙ + s Pₙ′)) + ẑ Σ Jₙ ρⁿ Pₙ′ ],
//!
//! with ρ = R/r and s = z/r; Pₙ and Pₙ′ come from Bonnet's recursion
//! n Pₙ = (2n − 1) s Pₙ₋₁ − (n − 1) Pₙ₋₂ and Pₙ′ = n Pₙ₋₁ + s Pₙ₋₁′.
//! A zonal field does not depend on longitude, so the position may be given
//! in any frame whose z axis is the body's axis (for the Earth, the celestial
//! intermediate pole: the GCRS, TEME or a true-of-date frame to the accuracy of
//! precession and nutation over the arc).
//!
//! Author: Moroya Sakamoto

use alice_det_math::sqrt64;

/// An axisymmetric gravity field: μ, the reference radius and the zonal
/// coefficients J2, J3, … (unnormalized).
#[derive(Debug, Clone, Copy)]
pub struct ZonalGravity<'a> {
    /// Gravitational parameter μ (km³/s²).
    pub mu: f64,
    /// Reference radius R of the coefficients (km).
    pub radius: f64,
    /// J2, J3, … in order (an empty slice is the point mass).
    pub zonals: &'a [f64],
}

impl ZonalGravity<'_> {
    /// Acceleration (km/s²) at position `r` (km). NaN at the centre, where
    /// the field is singular.
    #[must_use]
    pub fn acceleration(&self, r: &[f64; 3]) -> [f64; 3] {
        let r2 = r[2] * r[2] + (r[0] * r[0] + r[1] * r[1]);
        let rn = sqrt64(r2);
        if rn == 0.0 {
            return [f64::NAN; 3];
        }
        let s = r[2] / rn;
        let rho = self.radius / rn;
        let (mut p_prev, mut p) = (1.0, s);
        let mut dp = 1.0;
        let mut rho_n = rho;
        let mut radial = 1.0;
        let mut axial = 0.0;
        for (k, &j) in self.zonals.iter().enumerate() {
            let n = (k + 2) as f64;
            let p_next = ((2.0 * n - 1.0) * s * p - (n - 1.0) * p_prev) / n;
            dp = n * p + s * dp;
            p_prev = p;
            p = p_next;
            rho_n *= rho;
            radial -= j * rho_n * ((n + 1.0) * p + s * dp);
            axial += j * rho_n * dp;
        }
        let f = -self.mu / r2;
        let g = radial / rn;
        [f * (g * r[0]), f * (g * r[1]), f * (g * r[2] + axial)]
    }

    /// Potential U (km²/s²) at position `r` (km), positive, with the sign
    /// convention a = ∇U. NaN at the centre.
    #[must_use]
    pub fn potential(&self, r: &[f64; 3]) -> f64 {
        let rn = sqrt64(r[2] * r[2] + (r[0] * r[0] + r[1] * r[1]));
        if rn == 0.0 {
            return f64::NAN;
        }
        let s = r[2] / rn;
        let rho = self.radius / rn;
        let (mut p_prev, mut p) = (1.0, s);
        let mut rho_n = rho;
        let mut sum = 1.0;
        for (k, &j) in self.zonals.iter().enumerate() {
            let n = (k + 2) as f64;
            let p_next = ((2.0 * n - 1.0) * s * p - (n - 1.0) * p_prev) / n;
            p_prev = p;
            p = p_next;
            rho_n *= rho;
            sum -= j * rho_n * p;
        }
        self.mu / rn * sum
    }
}

/// The perturbing acceleration (km/s²) of a third body of gravitational
/// parameter `mu` (km³/s²) at position `s` (km) on a satellite at `r` (km),
/// both from the central body's centre:
/// a = −μ ((r − s)/|r − s|³ + s/|s|³).
///
/// The two terms nearly cancel when |r| ≪ |s|; Battin's form avoids the
/// cancellation: with q = r·(r − 2s)/|s|², |r − s|² = |s|²(1 + q) and
/// (1 + q)^{3/2} − 1 = f(q) = q (3 + 3q + q²) / (1 + (1 + q)^{3/2}), so
/// a = −μ/|r − s|³ (r + f(q) s). NaN when the body is at the centre or at the
/// satellite.
#[must_use]
pub fn third_body(r: [f64; 3], s: [f64; 3], mu: f64) -> [f64; 3] {
    let ss = s[0] * s[0] + s[1] * s[1] + s[2] * s[2];
    let r2s = [r[0] - 2.0 * s[0], r[1] - 2.0 * s[1], r[2] - 2.0 * s[2]];
    let q = (r[0] * r2s[0] + r[1] * r2s[1] + r[2] * r2s[2]) / ss;
    let opq = 1.0 + q;
    // the body at the satellite gives r + f(−1) s = 0 against k = −∞, and the
    // body at the centre q = ∞ and d2 = 0·∞: both NaN without a branch
    let d2 = ss * opq;
    let opq32 = opq * sqrt64(opq);
    let f = q * (3.0 + q * (3.0 + q)) / (1.0 + opq32);
    let k = -mu / (d2 * sqrt64(d2));
    [
        k * (r[0] + f * s[0]),
        k * (r[1] + f * s[1]),
        k * (r[2] + f * s[2]),
    ]
}
