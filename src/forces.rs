//! A time-dependent force model for an Earth satellite, and its RK4
//! propagation.
//!
//! The state is in the GCRS (km, km/s) with time counted in TT seconds from
//! the model's epoch. At every evaluation the model takes, at that instant:
//!
//! - the zonal field ([`ZonalGravity`]) and the drag ([`atmospheric_drag`]) in
//!   the celestial intermediate frame (CIRS, z along the celestial intermediate
//!   pole) from `alice_ephemeris::frames::c2i06a`, so precession and nutation
//!   turn the field's axis with the Earth's, and the rotating atmosphere
//!   turns about that pole;
//! - the Sun (the negated heliocentric Earth of ERFA `eraEpv00`) and the Moon
//!   (ERFA `eraMoon98`) from `alice_ephemeris::sunmoon`, for their
//!   perturbations ([`third_body`]) and for the radiation pressure
//!   ([`solar_radiation_pressure`]) with the Earth's and the Moon's conical
//!   shadows.
//!
//! TT stands in for TDB in the Sun and Moon models (|TDB − TT| ≤ 1.7 ms).
//!
//! Author: Moroya Sakamoto

use alice_ephemeris::frames::{c2i06a, mtxv, mxv};
use alice_ephemeris::sunmoon::{epv00, moon98};
use alice_ephemeris::Jd2;

use crate::drag::{atmospheric_drag, DragError, EARTH_ROTATION_RATE};
use crate::eclipse::{EARTH_EQUATORIAL_RADIUS_KM, MOON_RADIUS_KM};
use crate::gravity::{third_body, ZonalGravity};
use crate::orbit::SpacecraftState;
use crate::srp::{solar_radiation_pressure, Occulter, SrpError};

/// The astronomical unit (IAU 2012), km.
const AU_KM: f64 = 149_597_870.7;

/// Radiation pressure parameters.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct SrpParams {
    /// Reflectivity coefficient Cr.
    pub cr: f64,
    /// Area-to-mass ratio, m²/kg.
    pub area_to_mass: f64,
}

/// Drag parameters, with the caller's density model.
#[derive(Clone, Copy)]
pub struct DragParams<'a> {
    /// Drag coefficient Cd.
    pub cd: f64,
    /// Area-to-mass ratio, m²/kg.
    pub area_to_mass: f64,
    /// Density (kg/m³) at a TT instant and a CIRS position (km).
    pub density: &'a dyn Fn(Jd2, [f64; 3]) -> f64,
}

impl core::fmt::Debug for DragParams<'_> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("DragParams")
            .field("cd", &self.cd)
            .field("area_to_mass", &self.area_to_mass)
            .finish_non_exhaustive()
    }
}

/// The forces on an Earth satellite.
#[derive(Debug, Clone, Copy)]
pub struct ForceModel<'a> {
    /// TT Julian date of t = 0.
    pub epoch: Jd2,
    /// The Earth's field (μ, radius, zonal coefficients).
    pub gravity: ZonalGravity<'a>,
    /// The Sun's gravitational parameter (km³/s²), or no solar perturbation.
    pub sun_mu: Option<f64>,
    /// The Moon's gravitational parameter (km³/s²), or no lunar perturbation.
    pub moon_mu: Option<f64>,
    /// Radiation pressure, or none.
    pub srp: Option<SrpParams>,
    /// Drag, or none.
    pub drag: Option<DragParams<'a>>,
}

/// A body whose position the model takes from an ephemeris.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Body {
    /// The Sun (`alice_ephemeris::sunmoon::epv00`, 1900–2100).
    Sun,
    /// The Moon (`alice_ephemeris::sunmoon::moon98`, 1950–2050).
    Moon,
}

/// Why a force could not be evaluated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ForceError {
    /// The radiation pressure (an input, or the satellite inside a body).
    Srp(SrpError),
    /// The drag (an input or the density).
    Drag(DragError),
    /// The instant is outside the span of the ephemeris of a body the model
    /// uses (the position there would be an extrapolation of the series).
    OutOfEphemerisRange(Body),
}

impl core::fmt::Display for ForceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Srp(e) => write!(f, "radiation pressure: {e}"),
            Self::Drag(e) => write!(f, "drag: {e}"),
            Self::OutOfEphemerisRange(Body::Sun) => {
                f.write_str("the instant is outside the Sun's ephemeris (1900–2100)")
            }
            Self::OutOfEphemerisRange(Body::Moon) => {
                f.write_str("the instant is outside the Moon's ephemeris (1950–2050)")
            }
        }
    }
}

impl std::error::Error for ForceError {}

fn add(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale(a: [f64; 3], k: f64) -> [f64; 3] {
    [a[0] * k, a[1] * k, a[2] * k]
}

impl ForceModel<'_> {
    /// The TT instant `t_s` seconds after the epoch.
    #[must_use]
    pub fn instant(&self, t_s: f64) -> Jd2 {
        Jd2 {
            jd1: self.epoch.jd1,
            jd2: self.epoch.jd2 + t_s / 86_400.0,
        }
    }

    /// The acceleration (km/s², GCRS) at `t_s` seconds after the epoch, at
    /// GCRS position `r` (km) and velocity `v` (km/s).
    pub fn acceleration(&self, t_s: f64, r: [f64; 3], v: [f64; 3]) -> Result<[f64; 3], ForceError> {
        let tt = self.instant(t_s);
        let m = c2i06a(tt);
        let rc = mxv(&m, r);
        let mut ac = self.gravity.acceleration(&rc);
        if let Some(d) = &self.drag {
            let vc = mxv(&m, v);
            let rho = (d.density)(tt, rc);
            let ad = atmospheric_drag(rc, vc, rho, d.cd, d.area_to_mass, EARTH_ROTATION_RATE)
                .map_err(ForceError::Drag)?;
            ac = add(ac, ad);
        }
        let mut a = mtxv(&m, ac);
        let sun = if self.sun_mu.is_some() || self.srp.is_some() {
            let e = epv00(tt);
            if !e.in_range {
                return Err(ForceError::OutOfEphemerisRange(Body::Sun));
            }
            Some(scale(e.sun_geocentric()[0], AU_KM))
        } else {
            None
        };
        let moon = if self.moon_mu.is_some() || self.srp.is_some() {
            let m = moon98(tt);
            if !m.in_range {
                return Err(ForceError::OutOfEphemerisRange(Body::Moon));
            }
            Some(scale(m.geocentric[0], AU_KM))
        } else {
            None
        };
        if let (Some(mu), Some(s)) = (self.sun_mu, sun) {
            a = add(a, third_body(r, s, mu));
        }
        if let (Some(mu), Some(mo)) = (self.moon_mu, moon) {
            a = add(a, third_body(r, mo, mu));
        }
        if let (Some(p), Some(s), Some(mo)) = (self.srp, sun, moon) {
            let occulters = [
                Occulter {
                    position: [0.0; 3],
                    radius: EARTH_EQUATORIAL_RADIUS_KM,
                },
                Occulter {
                    position: mo,
                    radius: MOON_RADIUS_KM,
                },
            ];
            let ap = solar_radiation_pressure(r, s, &occulters, p.cr, p.area_to_mass)
                .map_err(ForceError::Srp)?;
            a = add(a, ap);
        }
        Ok(a)
    }

    /// Propagate `initial` (GCRS, at t = 0) by `steps` RK4 steps of `dt_s`
    /// seconds. Returns the states at each step, the initial one first.
    pub fn propagate(
        &self,
        initial: &SpacecraftState,
        dt_s: f64,
        steps: usize,
    ) -> Result<Vec<SpacecraftState>, ForceError> {
        let deriv = |t: f64, s: &[f64; 6]| -> Result<[f64; 6], ForceError> {
            let a = self.acceleration(t, [s[0], s[1], s[2]], [s[3], s[4], s[5]])?;
            Ok([s[3], s[4], s[5], a[0], a[1], a[2]])
        };
        let axpy = |s: &[f64; 6], k: &[f64; 6], h: f64| -> [f64; 6] {
            core::array::from_fn(|i| s[i] + k[i] * h)
        };
        let p = initial.position_km;
        let v = initial.velocity_km_s;
        let mut s = [p[0], p[1], p[2], v[0], v[1], v[2]];
        let dt_ns = (dt_s * 1e9) as u64;
        let mut out = Vec::with_capacity(steps + 1);
        out.extend([initial.clone()]);
        for k in 0..steps {
            let t = k as f64 * dt_s;
            let k1 = deriv(t, &s)?;
            let k2 = deriv(t + dt_s * 0.5, &axpy(&s, &k1, dt_s * 0.5))?;
            let k3 = deriv(t + dt_s * 0.5, &axpy(&s, &k2, dt_s * 0.5))?;
            let k4 = deriv(t + dt_s, &axpy(&s, &k3, dt_s))?;
            let h6 = dt_s / 6.0;
            s = core::array::from_fn(|i| {
                h6 * ((2.0 * k3[i] + (2.0 * k2[i] + k1[i])) + k4[i]) + s[i]
            });
            out.extend([SpacecraftState {
                position_km: [s[0], s[1], s[2]],
                velocity_km_s: [s[3], s[4], s[5]],
                timestamp_ns: initial.timestamp_ns + (k as u64 + 1) * dt_ns,
                fuel_kg: initial.fuel_kg,
            }]);
        }
        Ok(out)
    }
}
