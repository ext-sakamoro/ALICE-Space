//! SGP4 / SDP4 analytical propagation of NORAD two-line element sets.
//!
//! The algorithm is the 2006 revision by Vallado, Crawford, Hujsak and Kelso
//! ("Revisiting Spacetrack Report #3", AIAA 2006-6753), the reference model for
//! TLE / OMM mean elements: SGP4 for periods under 225 minutes, SDP4 (lunar /
//! solar perturbations and the 12 h / 24 h resonances) above. The code follows
//! the structure of Vallado's `sgp4unit` through its MIT-licensed Python
//! translation `python-sgp4` (Brandon Rhodes, commit 8126f773, see
//! THIRD-PARTY-NOTICES), keeping the operation order of every expression.
//!
//! # Determinism
//!
//! Every transcendental (sin, cos, atan2, pow) is evaluated by `alice-det-math`,
//! which composes them from IEEE basic operations only, so a propagation does not
//! depend on the platform's math library (`tests/sgp4_vallado.rs` pins the bits
//! of the whole verification walk). [`Sgp4::propagate`] is a pure function of the
//! element set and the time: the resonance integrator of SDP4 restarts from the
//! epoch on every call instead of keeping the last step between calls (the C++
//! reference keeps it), so the result does not depend on the order of earlier
//! calls.
//!
//! # Output frame and units
//!
//! Positions are km and velocities km/s in TEME (true equator, mean equinox of
//! date), the frame the element sets are fitted in. The time argument is
//! minutes since the element-set epoch.
//!
//! Author: Moroya Sakamoto

use alice_det_math::{atan2_64, cos64, powf64, sin64, sqrt64};

use crate::tle::Tle;

const PI: f64 = core::f64::consts::PI;
const TWOPI: f64 = 2.0 * PI;
const X2O3: f64 = 2.0 / 3.0;
const TEMP4: f64 = 1.5e-12;

/// The gravity constants a set of elements was fitted with.
///
/// NORAD fits TLEs with WGS-72; [`GravityModel::Wgs72`] is the model to use for
/// published element sets.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum GravityModel {
    /// WGS-72 with the legacy rounded `xke` of the 1980 report.
    Wgs72Old,
    /// WGS-72 (`mu` = 398600.8 km³/s², Earth radius 6378.135 km).
    Wgs72,
    /// WGS-84 (`mu` = 398600.5 km³/s², Earth radius 6378.137 km).
    Wgs84,
}

/// `(tumin, mu, radiusearthkm, xke, j2, j3, j4, j3oj2)`
fn gravity_constants(model: GravityModel) -> (f64, f64, f64, f64, f64, f64, f64, f64) {
    let (mu, radius, xke, j2, j3, j4) = match model {
        GravityModel::Wgs72Old => (
            398_600.799_64,
            6378.135,
            0.074_366_916_1,
            0.001_082_616,
            -2.538_81e-6,
            -1.655_97e-6,
        ),
        GravityModel::Wgs72 => {
            let (mu, r) = (398_600.8, 6378.135);
            (
                mu,
                r,
                60.0 / sqrt64(r * r * r / mu),
                0.001_082_616,
                -2.538_81e-6,
                -1.655_97e-6,
            )
        }
        GravityModel::Wgs84 => {
            let (mu, r) = (398_600.5, 6378.137);
            (
                mu,
                r,
                60.0 / sqrt64(r * r * r / mu),
                0.001_082_629_989_05,
                -2.532_153_06e-6,
                -1.610_987_61e-6,
            )
        }
    };
    (1.0 / xke, mu, radius, xke, j2, j3, j4, j3 / j2)
}

/// The two variants of the reference model.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OpsMode {
    /// As run by AFSPC: the sidereal time of the epoch from the 1970 series and
    /// the node of low-inclination deep-space orbits kept in `[0, 2π)`.
    Afspc,
    /// Vallado's improved mode: the IAU-1982 sidereal time of the epoch.
    Improved,
}

/// SGP4 mean elements in the units the model uses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MeanElements {
    /// Days since 1949 December 31 00:00 UT (Julian date − 2433281.5).
    pub epoch_days_since_1950: f64,
    /// Drag term B* (1 / Earth radii).
    pub bstar: f64,
    /// First derivative of the mean motion divided by two (rad / min²).
    pub ndot: f64,
    /// Second derivative of the mean motion divided by six (rad / min³).
    pub nddot: f64,
    /// Eccentricity.
    pub eccentricity: f64,
    /// Argument of perigee (rad).
    pub arg_perigee: f64,
    /// Inclination (rad).
    pub inclination: f64,
    /// Mean anomaly (rad).
    pub mean_anomaly: f64,
    /// Kozai mean motion (rad / min).
    pub mean_motion_kozai: f64,
    /// Right ascension of the ascending node (rad).
    pub raan: f64,
}

impl MeanElements {
    /// The mean elements of a TLE (angles in rad, mean motion in rev/day and
    /// its derivatives in rev/day² and rev/day³ as the TLE prints them).
    #[must_use]
    pub fn from_tle(tle: &Tle) -> Self {
        let xpdotp = 1440.0 / TWOPI;
        Self {
            epoch_days_since_1950: epoch_days_since_1950(tle.epoch_year, tle.epoch_day),
            bstar: tle.bstar,
            ndot: tle.mean_motion_dot / (xpdotp * 1440.0),
            nddot: tle.mean_motion_ddot / (xpdotp * 1440.0 * 1440.0),
            eccentricity: tle.eccentricity,
            arg_perigee: tle.arg_perigee_rad,
            inclination: tle.inclination_rad,
            mean_anomaly: tle.mean_anomaly_rad,
            mean_motion_kozai: tle.mean_motion_rev_per_day / xpdotp,
            raan: tle.raan_rad,
        }
    }
}

/// Days since 1949 December 31 00:00 of day `day_of_year` (1.0 = January 1
/// 00:00) of `year`.
#[must_use]
pub fn epoch_days_since_1950(year: u32, day_of_year: f64) -> f64 {
    let y = f64::from(year);
    // Julian date of January 1 00:00 of `year` (Vallado's jday with m = 1, d = 1)
    let jd_jan1 =
        367.0 * y - (7.0 * y / 4.0).floor() + (275.0 / 9.0f64).floor() + 1.0 + 1_721_013.5;
    jd_jan1 + (day_of_year - 1.0) - 2_433_281.5
}

/// Why a propagation has no state, with the value that failed the check (the
/// error codes 1 to 6 of the reference model).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Sgp4Error {
    /// Code 1: the mean eccentricity left `[-0.001, 1)` (or is NaN).
    MeanEccentricity(f64),
    /// Code 2: the mean motion is not positive (or is NaN, e.g. from a zero or
    /// non-finite input mean motion).
    MeanMotion(f64),
    /// Code 3: the perturbed eccentricity of SDP4 left `[0, 1]`.
    PerturbedEccentricity(f64),
    /// Code 4: the semi-latus rectum is negative (or NaN).
    SemiLatusRectum(f64),
    /// Code 6: the orbit radius is below one Earth radius (the satellite has
    /// decayed), or NaN.
    Decayed(f64),
}

impl Sgp4Error {
    /// The error code of the reference model.
    #[must_use]
    pub const fn code(&self) -> u8 {
        match self {
            Self::MeanEccentricity(_) => 1,
            Self::MeanMotion(_) => 2,
            Self::PerturbedEccentricity(_) => 3,
            Self::SemiLatusRectum(_) => 4,
            Self::Decayed(_) => 6,
        }
    }
}

impl core::fmt::Display for Sgp4Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::MeanEccentricity(e) => {
                write!(f, "mean eccentricity {e} not within range 0.0 <= e < 1.0")
            }
            Self::MeanMotion(n) => write!(f, "mean motion {n} is not positive"),
            Self::PerturbedEccentricity(e) => {
                write!(
                    f,
                    "perturbed eccentricity {e} not within range 0.0 <= e <= 1.0"
                )
            }
            Self::SemiLatusRectum(p) => write!(f, "semilatus rectum {p} is less than zero"),
            Self::Decayed(r) => {
                write!(f, "orbit radius {r} Earth radii: the satellite has decayed")
            }
        }
    }
}

impl std::error::Error for Sgp4Error {}

/// Why a set of mean elements cannot initialise the model (the reference code
/// divides by zero or takes the root of a negative number for these).
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum Sgp4InitError {
    /// The eccentricity is not in `[0, 1)`.
    Eccentricity(f64),
    /// The Kozai mean motion is not positive.
    MeanMotion(f64),
    /// A field is NaN or infinite (the field's name).
    NonFinite(&'static str),
}

impl core::fmt::Display for Sgp4InitError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Eccentricity(e) => write!(f, "eccentricity {e} not within range 0.0 <= e < 1.0"),
            Self::MeanMotion(n) => write!(f, "mean motion {n} rad/min is not positive"),
            Self::NonFinite(name) => write!(f, "{name} is not finite"),
        }
    }
}

impl std::error::Error for Sgp4InitError {}

/// A TEME state: position (km) and velocity (km/s).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct TemeState {
    /// Position (km).
    pub position_km: [f64; 3],
    /// Velocity (km/s).
    pub velocity_km_s: [f64; 3],
}

/// Lunar / solar periodic coefficients of SDP4 (`dscom` / `dpper`).
#[derive(Clone, Copy, Debug, Default)]
struct DeepPeriodics {
    e3: f64,
    ee2: f64,
    peo: f64,
    pgho: f64,
    pho: f64,
    pinco: f64,
    plo: f64,
    se2: f64,
    se3: f64,
    sgh2: f64,
    sgh3: f64,
    sgh4: f64,
    sh2: f64,
    sh3: f64,
    si2: f64,
    si3: f64,
    sl2: f64,
    sl3: f64,
    sl4: f64,
    xgh2: f64,
    xgh3: f64,
    xgh4: f64,
    xh2: f64,
    xh3: f64,
    xi2: f64,
    xi3: f64,
    xl2: f64,
    xl3: f64,
    xl4: f64,
    zmol: f64,
    zmos: f64,
}

/// Secular and resonance coefficients of SDP4 (`dsinit` / `dspace`).
#[derive(Clone, Copy, Debug, Default)]
struct DeepSecular {
    irez: u8,
    d2201: f64,
    d2211: f64,
    d3210: f64,
    d3222: f64,
    d4410: f64,
    d4422: f64,
    d5220: f64,
    d5232: f64,
    d5421: f64,
    d5433: f64,
    dedt: f64,
    del1: f64,
    del2: f64,
    del3: f64,
    didt: f64,
    dmdt: f64,
    dnodt: f64,
    domdt: f64,
    xfact: f64,
    xlamo: f64,
}

/// An initialised SGP4 / SDP4 propagator for one element set.
#[derive(Clone, Debug)]
pub struct Sgp4 {
    opsmode: OpsMode,
    radiusearthkm: f64,
    xke: f64,
    j2: f64,
    j3oj2: f64,
    bstar: f64,
    ecco: f64,
    argpo: f64,
    inclo: f64,
    mo: f64,
    nodeo: f64,
    no_unkozai: f64,
    deep: bool,
    isimp: bool,
    aycof: f64,
    con41: f64,
    cc1: f64,
    cc4: f64,
    cc5: f64,
    d2: f64,
    d3: f64,
    d4: f64,
    delmo: f64,
    eta: f64,
    argpdot: f64,
    omgcof: f64,
    sinmao: f64,
    t2cof: f64,
    t3cof: f64,
    t4cof: f64,
    t5cof: f64,
    x1mth2: f64,
    x7thm1: f64,
    mdot: f64,
    nodedot: f64,
    xlcof: f64,
    xmcof: f64,
    nodecf: f64,
    gsto: f64,
    periodics: DeepPeriodics,
    secular: DeepSecular,
}

/// IAU-1982 Greenwich mean sidereal time (rad) of the UT1 Julian date `jdut1`.
#[must_use]
pub fn gstime(jdut1: f64) -> f64 {
    let tut1 = (jdut1 - 2_451_545.0) / 36525.0;
    let mut temp = -6.2e-6 * tut1 * tut1 * tut1
        + 0.093_104 * tut1 * tut1
        + (876_600.0 * 3600.0 + 8_640_184.812_866) * tut1
        + 67_310.548_41;
    temp = temp * (PI / 180.0) / 240.0 % TWOPI;
    if temp < 0.0 {
        temp += TWOPI;
    }
    temp
}

/// `x % 2π` with the sign handling of the reference code (truncation towards
/// zero, like C's `fmod`).
fn fmod_twopi(x: f64) -> f64 {
    x % TWOPI
}

/// Python's `%` with a positive modulus (the result takes the sign of the
/// modulus): used where the reference applies it to a value that may be
/// negative.
fn mod_floor_twopi(x: f64) -> f64 {
    let r = x % TWOPI;
    if r < 0.0 {
        r + TWOPI
    } else {
        r
    }
}

impl DeepPeriodics {
    /// `dpper`: the lunar / solar long-period periodics at `t` (zero at the
    /// epoch when `init`) applied to `(ep, inclp, nodep, argpp, mp)`.
    #[allow(clippy::too_many_arguments)]
    fn apply(
        &self,
        t: f64,
        init: bool,
        opsmode: OpsMode,
        mut ep: f64,
        mut inclp: f64,
        mut nodep: f64,
        mut argpp: f64,
        mut mp: f64,
    ) -> (f64, f64, f64, f64, f64) {
        let zns = 1.194_59e-5;
        let zes = 0.016_75;
        let znl = 1.583_521_8e-4;
        let zel = 0.0549;

        let mut zm = self.zmos + zns * t;
        if init {
            zm = self.zmos;
        }
        let mut zf = zm + 2.0 * zes * sin64(zm);
        let mut sinzf = sin64(zf);
        let mut f2 = 0.5 * sinzf * sinzf - 0.25;
        let mut f3 = -0.5 * sinzf * cos64(zf);
        let ses = self.se2 * f2 + self.se3 * f3;
        let sis = self.si2 * f2 + self.si3 * f3;
        let sls = self.sl2 * f2 + self.sl3 * f3 + self.sl4 * sinzf;
        let sghs = self.sgh2 * f2 + self.sgh3 * f3 + self.sgh4 * sinzf;
        let shs = self.sh2 * f2 + self.sh3 * f3;
        zm = self.zmol + znl * t;
        if init {
            zm = self.zmol;
        }
        zf = zm + 2.0 * zel * sin64(zm);
        sinzf = sin64(zf);
        f2 = 0.5 * sinzf * sinzf - 0.25;
        f3 = -0.5 * sinzf * cos64(zf);
        let sel = self.ee2 * f2 + self.e3 * f3;
        let sil = self.xi2 * f2 + self.xi3 * f3;
        let sll = self.xl2 * f2 + self.xl3 * f3 + self.xl4 * sinzf;
        let sghl = self.xgh2 * f2 + self.xgh3 * f3 + self.xgh4 * sinzf;
        let shll = self.xh2 * f2 + self.xh3 * f3;
        let mut pe = ses + sel;
        let mut pinc = sis + sil;
        let mut pl = sls + sll;
        let mut pgh = sghs + sghl;
        let mut ph = shs + shll;

        if !init {
            pe -= self.peo;
            pinc -= self.pinco;
            pl -= self.plo;
            pgh -= self.pgho;
            ph -= self.pho;
            inclp += pinc;
            ep += pe;
            let sinip = sin64(inclp);
            let cosip = cos64(inclp);
            if inclp >= 0.2 {
                ph /= sinip;
                pgh -= cosip * ph;
                argpp += pgh;
                nodep += ph;
                mp += pl;
            } else {
                // Lyddane modification for low inclinations
                let sinop = sin64(nodep);
                let cosop = cos64(nodep);
                let mut alfdp = sinip * sinop;
                let mut betdp = sinip * cosop;
                let dalf = ph * cosop + pinc * cosip * sinop;
                let dbet = -ph * sinop + pinc * cosip * cosop;
                alfdp += dalf;
                betdp += dbet;
                nodep = fmod_twopi(nodep);
                if nodep < 0.0 && opsmode == OpsMode::Afspc {
                    nodep += TWOPI;
                }
                let xls = mp + argpp + pl + pgh + (cosip - pinc * sinip) * nodep;
                let xnoh = nodep;
                nodep = atan2_64(alfdp, betdp);
                if nodep < 0.0 && opsmode == OpsMode::Afspc {
                    nodep += TWOPI;
                }
                if (xnoh - nodep).abs() > PI {
                    if nodep < xnoh {
                        nodep += TWOPI;
                    } else {
                        nodep -= TWOPI;
                    }
                }
                mp += pl;
                argpp = xls - mp - cosip * nodep;
            }
        }
        (ep, inclp, nodep, argpp, mp)
    }
}

/// The intermediate quantities of `dscom` that `dsinit` reads.
#[derive(Clone, Copy, Debug, Default)]
struct Dscom {
    sinim: f64,
    cosim: f64,
    emsq: f64,
    s1: f64,
    s2: f64,
    s3: f64,
    s4: f64,
    s5: f64,
    ss1: f64,
    ss2: f64,
    ss3: f64,
    ss4: f64,
    ss5: f64,
    sz1: f64,
    sz3: f64,
    sz11: f64,
    sz13: f64,
    sz21: f64,
    sz23: f64,
    sz31: f64,
    sz33: f64,
    z1: f64,
    z3: f64,
    z11: f64,
    z13: f64,
    z21: f64,
    z23: f64,
    z31: f64,
    z33: f64,
    em: f64,
    nm: f64,
}

/// `dscom`: the lunar / solar terms of SDP4 at the epoch.
#[allow(clippy::too_many_lines, clippy::similar_names)]
fn dscom(
    epoch: f64,
    ep: f64,
    argpp: f64,
    tc: f64,
    inclp: f64,
    nodep: f64,
    np: f64,
) -> (Dscom, DeepPeriodics) {
    let zes = 0.016_75;
    let zel = 0.0549;
    let c1ss = 2.986_479_7e-6;
    let c1l = 4.796_806_5e-7;
    let zsinis = 0.397_854_16;
    let zcosis = 0.917_448_67;
    let zcosgs = 0.194_590_5;
    let zsings = -0.980_884_58;

    let nm = np;
    let em = ep;
    let snodm = sin64(nodep);
    let cnodm = cos64(nodep);
    let sinomm = sin64(argpp);
    let cosomm = cos64(argpp);
    let sinim = sin64(inclp);
    let cosim = cos64(inclp);
    let emsq = em * em;
    let betasq = 1.0 - emsq;
    let rtemsq = sqrt64(betasq);

    let day = epoch + 18261.5 + tc / 1440.0;
    let xnodce = mod_floor_twopi(4.523_602_0 - 9.242_202_9e-4 * day);
    let stem = sin64(xnodce);
    let ctem = cos64(xnodce);
    let zcosil = 0.913_751_64 - 0.035_680_96 * ctem;
    let zsinil = sqrt64(1.0 - zcosil * zcosil);
    let zsinhl = 0.089_683_511 * stem / zsinil;
    let zcoshl = sqrt64(1.0 - zsinhl * zsinhl);
    let gam = 5.835_151_4 + 0.001_944_368 * day;
    let mut zx = 0.397_854_16 * stem / zsinil;
    let zy = zcoshl * ctem + 0.917_448_67 * zsinhl * stem;
    zx = atan2_64(zx, zy);
    zx = gam + zx - xnodce;
    let zcosgl = cos64(zx);
    let zsingl = sin64(zx);

    let mut zcosg = zcosgs;
    let mut zsing = zsings;
    let mut zcosi = zcosis;
    let mut zsini = zsinis;
    let mut zcosh = cnodm;
    let mut zsinh = snodm;
    let mut cc = c1ss;
    let xnoi = 1.0 / nm;

    let mut d = Dscom::default();
    let (mut s1, mut s2, mut s3, mut s4, mut s5, mut s6, mut s7) =
        (0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut z1, mut z2, mut z3) = (0.0, 0.0, 0.0);
    let (mut z11, mut z12, mut z13, mut z21, mut z22, mut z23) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let (mut z31, mut z32, mut z33) = (0.0, 0.0, 0.0);
    let (mut ss6, mut ss7) = (0.0, 0.0);
    let (mut sz2, mut sz12, mut sz22, mut sz32) = (0.0, 0.0, 0.0, 0.0);

    for lsflg in 1..=2 {
        let a1 = zcosg * zcosh + zsing * zcosi * zsinh;
        let a3 = -zsing * zcosh + zcosg * zcosi * zsinh;
        let a7 = -zcosg * zsinh + zsing * zcosi * zcosh;
        let a8 = zsing * zsini;
        let a9 = zsing * zsinh + zcosg * zcosi * zcosh;
        let a10 = zcosg * zsini;
        let a2 = cosim * a7 + sinim * a8;
        let a4 = cosim * a9 + sinim * a10;
        let a5 = -sinim * a7 + cosim * a8;
        let a6 = -sinim * a9 + cosim * a10;

        let x1 = a1 * cosomm + a2 * sinomm;
        let x2 = a3 * cosomm + a4 * sinomm;
        let x3 = -a1 * sinomm + a2 * cosomm;
        let x4 = -a3 * sinomm + a4 * cosomm;
        let x5 = a5 * sinomm;
        let x6 = a6 * sinomm;
        let x7 = a5 * cosomm;
        let x8 = a6 * cosomm;

        z31 = 12.0 * x1 * x1 - 3.0 * x3 * x3;
        z32 = 24.0 * x1 * x2 - 6.0 * x3 * x4;
        z33 = 12.0 * x2 * x2 - 3.0 * x4 * x4;
        z1 = 3.0 * (a1 * a1 + a2 * a2) + z31 * emsq;
        z2 = 6.0 * (a1 * a3 + a2 * a4) + z32 * emsq;
        z3 = 3.0 * (a3 * a3 + a4 * a4) + z33 * emsq;
        z11 = -6.0 * a1 * a5 + emsq * (-24.0 * x1 * x7 - 6.0 * x3 * x5);
        z12 = -6.0 * (a1 * a6 + a3 * a5)
            + emsq * (-24.0 * (x2 * x7 + x1 * x8) - 6.0 * (x3 * x6 + x4 * x5));
        z13 = -6.0 * a3 * a6 + emsq * (-24.0 * x2 * x8 - 6.0 * x4 * x6);
        z21 = 6.0 * a2 * a5 + emsq * (24.0 * x1 * x5 - 6.0 * x3 * x7);
        z22 = 6.0 * (a4 * a5 + a2 * a6)
            + emsq * (24.0 * (x2 * x5 + x1 * x6) - 6.0 * (x4 * x7 + x3 * x8));
        z23 = 6.0 * a4 * a6 + emsq * (24.0 * x2 * x6 - 6.0 * x4 * x8);
        z1 = z1 + z1 + betasq * z31;
        z2 = z2 + z2 + betasq * z32;
        z3 = z3 + z3 + betasq * z33;
        s3 = cc * xnoi;
        s2 = -0.5 * s3 / rtemsq;
        s4 = s3 * rtemsq;
        s1 = -15.0 * em * s4;
        s5 = x1 * x3 + x2 * x4;
        s6 = x2 * x3 + x1 * x4;
        s7 = x2 * x4 - x1 * x3;

        if lsflg == 1 {
            d.ss1 = s1;
            d.ss2 = s2;
            d.ss3 = s3;
            d.ss4 = s4;
            d.ss5 = s5;
            ss6 = s6;
            ss7 = s7;
            d.sz1 = z1;
            sz2 = z2;
            d.sz3 = z3;
            d.sz11 = z11;
            sz12 = z12;
            d.sz13 = z13;
            d.sz21 = z21;
            sz22 = z22;
            d.sz23 = z23;
            d.sz31 = z31;
            sz32 = z32;
            d.sz33 = z33;
            zcosg = zcosgl;
            zsing = zsingl;
            zcosi = zcosil;
            zsini = zsinil;
            zcosh = zcoshl * cnodm + zsinhl * snodm;
            zsinh = snodm * zcoshl - cnodm * zsinhl;
            cc = c1l;
        }
    }

    let p = DeepPeriodics {
        zmol: mod_floor_twopi(4.719_967_2 + 0.229_971_50 * day - gam),
        zmos: mod_floor_twopi(6.256_583_7 + 0.017_201_977 * day),
        se2: 2.0 * d.ss1 * ss6,
        se3: 2.0 * d.ss1 * ss7,
        si2: 2.0 * d.ss2 * sz12,
        si3: 2.0 * d.ss2 * (d.sz13 - d.sz11),
        sl2: -2.0 * d.ss3 * sz2,
        sl3: -2.0 * d.ss3 * (d.sz3 - d.sz1),
        sl4: -2.0 * d.ss3 * (-21.0 - 9.0 * emsq) * zes,
        sgh2: 2.0 * d.ss4 * sz32,
        sgh3: 2.0 * d.ss4 * (d.sz33 - d.sz31),
        sgh4: -18.0 * d.ss4 * zes,
        sh2: -2.0 * d.ss2 * sz22,
        sh3: -2.0 * d.ss2 * (d.sz23 - d.sz21),
        ee2: 2.0 * s1 * s6,
        e3: 2.0 * s1 * s7,
        xi2: 2.0 * s2 * z12,
        xi3: 2.0 * s2 * (z13 - z11),
        xl2: -2.0 * s3 * z2,
        xl3: -2.0 * s3 * (z3 - z1),
        xl4: -2.0 * s3 * (-21.0 - 9.0 * emsq) * zel,
        xgh2: 2.0 * s4 * z32,
        xgh3: 2.0 * s4 * (z33 - z31),
        xgh4: -18.0 * s4 * zel,
        xh2: -2.0 * s2 * z22,
        xh3: -2.0 * s2 * (z23 - z21),
        peo: 0.0,
        pgho: 0.0,
        pho: 0.0,
        pinco: 0.0,
        plo: 0.0,
    };

    d.sinim = sinim;
    d.cosim = cosim;
    d.emsq = emsq;
    d.s1 = s1;
    d.s2 = s2;
    d.s3 = s3;
    d.s4 = s4;
    d.s5 = s5;
    d.z1 = z1;
    d.z3 = z3;
    d.z11 = z11;
    d.z13 = z13;
    d.z21 = z21;
    d.z23 = z23;
    d.z31 = z31;
    d.z33 = z33;
    d.em = em;
    d.nm = nm;
    (d, p)
}

impl Sgp4 {
    /// Initialise the propagator for `elements` (`sgp4init`).
    ///
    /// An element set that initialises but that the model cannot propagate at
    /// some time reports that from [`Self::propagate`] (also at the epoch).
    ///
    /// # Errors
    ///
    /// [`Sgp4InitError`] for a non-finite field, an eccentricity outside
    /// `[0, 1)` or a mean motion that is not positive.
    #[allow(clippy::too_many_lines, clippy::similar_names)]
    pub fn new(
        elements: &MeanElements,
        gravity: GravityModel,
        opsmode: OpsMode,
    ) -> Result<Self, Sgp4InitError> {
        let fields = [
            ("epoch", elements.epoch_days_since_1950),
            ("bstar", elements.bstar),
            ("ndot", elements.ndot),
            ("nddot", elements.nddot),
            ("eccentricity", elements.eccentricity),
            ("argument of perigee", elements.arg_perigee),
            ("inclination", elements.inclination),
            ("mean anomaly", elements.mean_anomaly),
            ("mean motion", elements.mean_motion_kozai),
            ("right ascension of the node", elements.raan),
        ];
        if let Some((name, _)) = fields.iter().find(|(_, x)| !x.is_finite()) {
            return Err(Sgp4InitError::NonFinite(name));
        }
        if !(0.0..1.0).contains(&elements.eccentricity) {
            return Err(Sgp4InitError::Eccentricity(elements.eccentricity));
        }
        if elements.mean_motion_kozai <= 0.0 {
            return Err(Sgp4InitError::MeanMotion(elements.mean_motion_kozai));
        }
        let (tumin, _mu, radiusearthkm, xke, j2, _j3, j4, j3oj2) = gravity_constants(gravity);
        let epoch = elements.epoch_days_since_1950;
        let mut s = Self {
            opsmode,
            radiusearthkm,
            xke,
            j2,
            j3oj2,
            bstar: elements.bstar,
            ecco: elements.eccentricity,
            argpo: elements.arg_perigee,
            inclo: elements.inclination,
            mo: elements.mean_anomaly,
            nodeo: elements.raan,
            no_unkozai: 0.0,
            deep: false,
            isimp: false,
            aycof: 0.0,
            con41: 0.0,
            cc1: 0.0,
            cc4: 0.0,
            cc5: 0.0,
            d2: 0.0,
            d3: 0.0,
            d4: 0.0,
            delmo: 0.0,
            eta: 0.0,
            argpdot: 0.0,
            omgcof: 0.0,
            sinmao: 0.0,
            t2cof: 0.0,
            t3cof: 0.0,
            t4cof: 0.0,
            t5cof: 0.0,
            x1mth2: 0.0,
            x7thm1: 0.0,
            mdot: 0.0,
            nodedot: 0.0,
            xlcof: 0.0,
            xmcof: 0.0,
            nodecf: 0.0,
            gsto: 0.0,
            periodics: DeepPeriodics::default(),
            secular: DeepSecular::default(),
        };
        let _ = tumin;

        let ss = 78.0 / radiusearthkm + 1.0;
        let qzms2ttemp = (120.0 - 78.0) / radiusearthkm;
        let qzms2t = qzms2ttemp * qzms2ttemp * qzms2ttemp * qzms2ttemp;

        // initl
        let eccsq = s.ecco * s.ecco;
        let omeosq = 1.0 - eccsq;
        let rteosq = sqrt64(omeosq);
        let cosio = cos64(s.inclo);
        let cosio2 = cosio * cosio;
        let ak = powf64(xke / elements.mean_motion_kozai, X2O3);
        let d1 = 0.75 * j2 * (3.0 * cosio2 - 1.0) / (rteosq * omeosq);
        let mut del = d1 / (ak * ak);
        let adel = ak * (1.0 - del * del - del * (1.0 / 3.0 + 134.0 * del * del / 81.0));
        del = d1 / (adel * adel);
        s.no_unkozai = elements.mean_motion_kozai / (1.0 + del);
        let ao = powf64(xke / s.no_unkozai, X2O3);
        let sinio = sin64(s.inclo);
        let po = ao * omeosq;
        let con42 = 1.0 - 5.0 * cosio2;
        s.con41 = -con42 - cosio2 - cosio2;
        let posq = po * po;
        let rp = ao * (1.0 - s.ecco);
        s.gsto = match opsmode {
            OpsMode::Afspc => {
                let ts70 = epoch - 7305.0;
                let ds70 = (ts70 + 1.0e-8).floor();
                let tfrac = ts70 - ds70;
                let c1 = 1.720_279_169_407_036_2e-2;
                let thgr70 = 1.732_134_385_650_937_5;
                let fk5r = 5.075_514_194_322_695e-15;
                let c1p2p = c1 + TWOPI;
                let mut g =
                    mod_floor_twopi(thgr70 + c1 * ds70 + c1p2p * tfrac + ts70 * ts70 * fk5r);
                if g < 0.0 {
                    g += TWOPI;
                }
                g
            }
            OpsMode::Improved => gstime(epoch + 2_433_281.5),
        };

        if omeosq >= 0.0 || s.no_unkozai >= 0.0 {
            s.isimp = rp < 220.0 / radiusearthkm + 1.0;
            let mut sfour = ss;
            let mut qzms24 = qzms2t;
            let perige = (rp - 1.0) * radiusearthkm;
            if perige < 156.0 {
                sfour = perige - 78.0;
                if perige < 98.0 {
                    sfour = 20.0;
                }
                let qzms24temp = (120.0 - sfour) / radiusearthkm;
                qzms24 = qzms24temp * qzms24temp * qzms24temp * qzms24temp;
                sfour = sfour / radiusearthkm + 1.0;
            }
            let pinvsq = 1.0 / posq;
            let tsi = 1.0 / (ao - sfour);
            s.eta = ao * s.ecco * tsi;
            let etasq = s.eta * s.eta;
            let eeta = s.ecco * s.eta;
            let psisq = (1.0 - etasq).abs();
            let coef = qzms24 * powf64(tsi, 4.0);
            let coef1 = coef / powf64(psisq, 3.5);
            let cc2 = coef1
                * s.no_unkozai
                * (ao * (1.0 + 1.5 * etasq + eeta * (4.0 + etasq))
                    + 0.375 * j2 * tsi / psisq * s.con41 * (8.0 + 3.0 * etasq * (8.0 + etasq)));
            s.cc1 = s.bstar * cc2;
            let mut cc3 = 0.0;
            if s.ecco > 1.0e-4 {
                cc3 = -2.0 * coef * tsi * j3oj2 * s.no_unkozai * sinio / s.ecco;
            }
            s.x1mth2 = 1.0 - cosio2;
            s.cc4 = 2.0
                * s.no_unkozai
                * coef1
                * ao
                * omeosq
                * (s.eta * (2.0 + 0.5 * etasq) + s.ecco * (0.5 + 2.0 * etasq)
                    - j2 * tsi / (ao * psisq)
                        * (-3.0 * s.con41 * (1.0 - 2.0 * eeta + etasq * (1.5 - 0.5 * eeta))
                            + 0.75
                                * s.x1mth2
                                * (2.0 * etasq - eeta * (1.0 + etasq))
                                * cos64(2.0 * s.argpo)));
            s.cc5 = 2.0 * coef1 * ao * omeosq * (1.0 + 2.75 * (etasq + eeta) + eeta * etasq);
            let cosio4 = cosio2 * cosio2;
            let temp1 = 1.5 * j2 * pinvsq * s.no_unkozai;
            let temp2 = 0.5 * temp1 * j2 * pinvsq;
            let temp3 = -0.468_75 * j4 * pinvsq * pinvsq * s.no_unkozai;
            s.mdot = s.no_unkozai
                + 0.5 * temp1 * rteosq * s.con41
                + 0.0625 * temp2 * rteosq * (13.0 - 78.0 * cosio2 + 137.0 * cosio4);
            s.argpdot = -0.5 * temp1 * con42
                + 0.0625 * temp2 * (7.0 - 114.0 * cosio2 + 395.0 * cosio4)
                + temp3 * (3.0 - 36.0 * cosio2 + 49.0 * cosio4);
            let xhdot1 = -temp1 * cosio;
            s.nodedot = xhdot1
                + (0.5 * temp2 * (4.0 - 19.0 * cosio2) + 2.0 * temp3 * (3.0 - 7.0 * cosio2))
                    * cosio;
            let xpidot = s.argpdot + s.nodedot;
            s.omgcof = s.bstar * cc3 * cos64(s.argpo);
            s.xmcof = 0.0;
            if s.ecco > 1.0e-4 {
                s.xmcof = -X2O3 * coef * s.bstar / eeta;
            }
            s.nodecf = 3.5 * omeosq * xhdot1 * s.cc1;
            s.t2cof = 1.5 * s.cc1;
            s.xlcof = if (cosio + 1.0).abs() > 1.5e-12 {
                -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / (1.0 + cosio)
            } else {
                -0.25 * j3oj2 * sinio * (3.0 + 5.0 * cosio) / TEMP4
            };
            s.aycof = -0.5 * j3oj2 * sinio;
            let delmotemp = 1.0 + s.eta * cos64(s.mo);
            s.delmo = delmotemp * delmotemp * delmotemp;
            s.sinmao = sin64(s.mo);
            s.x7thm1 = 7.0 * cosio2 - 1.0;

            // deep space initialisation
            if 2.0 * PI / s.no_unkozai >= 225.0 {
                s.deep = true;
                s.isimp = true;
                let tc = 0.0;
                let inclm = s.inclo;
                let (d, p) = dscom(epoch, s.ecco, s.argpo, tc, s.inclo, s.nodeo, s.no_unkozai);
                s.periodics = p;
                let (ecco, inclo, nodeo, argpo, mo) = s
                    .periodics
                    .apply(0.0, true, opsmode, s.ecco, s.inclo, s.nodeo, s.argpo, s.mo);
                s.ecco = ecco;
                s.inclo = inclo;
                s.nodeo = nodeo;
                s.argpo = argpo;
                s.mo = mo;
                s.secular = s.dsinit(&d, tc, xpidot, eccsq, inclm);
            }

            if !s.isimp {
                let cc1sq = s.cc1 * s.cc1;
                s.d2 = 4.0 * ao * tsi * cc1sq;
                let temp = s.d2 * tsi * s.cc1 / 3.0;
                s.d3 = (17.0 * ao + sfour) * temp;
                s.d4 = 0.5 * temp * ao * tsi * (221.0 * ao + 31.0 * sfour) * s.cc1;
                s.t3cof = s.d2 + 2.0 * cc1sq;
                s.t4cof = 0.25 * (3.0 * s.d3 + s.cc1 * (12.0 * s.d2 + 10.0 * cc1sq));
                s.t5cof = 0.2
                    * (3.0 * s.d4
                        + 12.0 * s.cc1 * s.d3
                        + 6.0 * s.d2 * s.d2
                        + 15.0 * cc1sq * (2.0 * s.d2 + cc1sq));
            }
        }
        Ok(s)
    }

    /// `dsinit`: the secular rates and resonance terms of SDP4 (the terms the
    /// reference also applies at `t` = 0 change nothing there and are left out).
    #[allow(clippy::too_many_lines, clippy::similar_names)]
    fn dsinit(&self, d: &Dscom, tc: f64, xpidot: f64, eccsq: f64, inclm: f64) -> DeepSecular {
        let q22 = 1.789_167_9e-6;
        let q31 = 2.146_074_8e-6;
        let q33 = 2.212_301_5e-7;
        let root22 = 1.789_167_9e-6;
        let root44 = 7.363_695_3e-9;
        let root54 = 2.176_580_3e-9;
        let rptim = 4.375_269_088_011_3e-3;
        let root32 = 3.739_379_2e-7;
        let root52 = 1.142_863_9e-7;
        let znl = 1.583_521_8e-4;
        let zns = 1.194_59e-5;

        let nm = d.nm;
        let em = d.em;
        let mut emsq = d.emsq;
        let (sinim, cosim) = (d.sinim, d.cosim);

        let mut out = DeepSecular::default();
        if 0.003_490_658_5 < nm && nm < 0.005_235_987_7 {
            out.irez = 1;
        }
        if (0.008_26..=0.009_24).contains(&nm) && em >= 0.5 {
            out.irez = 2;
        }

        let ses = d.ss1 * zns * d.ss5;
        let sis = d.ss2 * zns * (d.sz11 + d.sz13);
        let sls = -zns * d.ss3 * (d.sz1 + d.sz3 - 14.0 - 6.0 * emsq);
        let sghs = d.ss4 * zns * (d.sz31 + d.sz33 - 6.0);
        let mut shs = -zns * d.ss2 * (d.sz21 + d.sz23);
        if !(5.235_987_7e-2..=PI - 5.235_987_7e-2).contains(&inclm) {
            shs = 0.0;
        }
        if sinim != 0.0 {
            shs /= sinim;
        }
        let sgs = sghs - cosim * shs;

        out.dedt = ses + d.s1 * znl * d.s5;
        out.didt = sis + d.s2 * znl * (d.z11 + d.z13);
        out.dmdt = sls - znl * d.s3 * (d.z1 + d.z3 - 14.0 - 6.0 * emsq);
        let sghl = d.s4 * znl * (d.z31 + d.z33 - 6.0);
        let mut shll = -znl * d.s2 * (d.z21 + d.z23);
        if !(5.235_987_7e-2..=PI - 5.235_987_7e-2).contains(&inclm) {
            shll = 0.0;
        }
        out.domdt = sgs + sghl;
        out.dnodt = shs;
        if sinim != 0.0 {
            out.domdt -= cosim / sinim * shll;
            out.dnodt += shll / sinim;
        }

        let theta = mod_floor_twopi(self.gsto + tc * rptim);

        if out.irez != 0 {
            let aonv = powf64(nm / self.xke, X2O3);
            if out.irez == 2 {
                // geopotential resonance for 12 hour orbits
                let cosisq = cosim * cosim;
                let em = self.ecco;
                emsq = eccsq;
                let eoc = em * emsq;
                let g201 = -0.306 - (em - 0.64) * 0.440;
                let (g211, g310, g322, g410, g422, g520);
                if em <= 0.65 {
                    g211 = 3.616 - 13.2470 * em + 16.2900 * emsq;
                    g310 = -19.302 + 117.3900 * em - 228.4190 * emsq + 156.5910 * eoc;
                    g322 = -18.9068 + 109.7927 * em - 214.6334 * emsq + 146.5816 * eoc;
                    g410 = -41.122 + 242.6940 * em - 471.0940 * emsq + 313.9530 * eoc;
                    g422 = -146.407 + 841.8800 * em - 1629.014 * emsq + 1083.4350 * eoc;
                    g520 = -532.114 + 3017.977 * em - 5740.032 * emsq + 3708.2760 * eoc;
                } else {
                    g211 = -72.099 + 331.819 * em - 508.738 * emsq + 266.724 * eoc;
                    g310 = -346.844 + 1582.851 * em - 2415.925 * emsq + 1246.113 * eoc;
                    g322 = -342.585 + 1554.908 * em - 2366.899 * emsq + 1215.972 * eoc;
                    g410 = -1052.797 + 4758.686 * em - 7193.992 * emsq + 3651.957 * eoc;
                    g422 = -3581.690 + 16178.110 * em - 24462.770 * emsq + 12422.520 * eoc;
                    g520 = if em > 0.715 {
                        -5149.66 + 29936.92 * em - 54087.36 * emsq + 31324.56 * eoc
                    } else {
                        1464.74 - 4664.75 * em + 3763.64 * emsq
                    };
                }
                let (g533, g521, g532);
                if em < 0.7 {
                    g533 = -919.22770 + 4988.6100 * em - 9064.7700 * emsq + 5542.21 * eoc;
                    g521 = -822.71072 + 4568.6173 * em - 8491.4146 * emsq + 5337.524 * eoc;
                    g532 = -853.66600 + 4690.2500 * em - 8624.7700 * emsq + 5341.4 * eoc;
                } else {
                    g533 = -37995.780 + 161616.52 * em - 229838.20 * emsq + 109377.94 * eoc;
                    g521 = -51752.104 + 218913.95 * em - 309468.16 * emsq + 146349.42 * eoc;
                    g532 = -40023.880 + 170470.89 * em - 242699.48 * emsq + 115605.82 * eoc;
                }
                let sini2 = sinim * sinim;
                let f220 = 0.75 * (1.0 + 2.0 * cosim + cosisq);
                let f221 = 1.5 * sini2;
                let f321 = 1.875 * sinim * (1.0 - 2.0 * cosim - 3.0 * cosisq);
                let f322 = -1.875 * sinim * (1.0 + 2.0 * cosim - 3.0 * cosisq);
                let f441 = 35.0 * sini2 * f220;
                let f442 = 39.3750 * sini2 * sini2;
                let f522 = 9.84375
                    * sinim
                    * (sini2 * (1.0 - 2.0 * cosim - 5.0 * cosisq)
                        + 0.333_333_33 * (-2.0 + 4.0 * cosim + 6.0 * cosisq));
                let f523 = sinim
                    * (4.921_875_12 * sini2 * (-2.0 - 4.0 * cosim + 10.0 * cosisq)
                        + 6.562_500_12 * (1.0 + 2.0 * cosim - 3.0 * cosisq));
                let f542 = 29.53125
                    * sinim
                    * (2.0 - 8.0 * cosim + cosisq * (-12.0 + 8.0 * cosim + 10.0 * cosisq));
                let f543 = 29.53125
                    * sinim
                    * (-2.0 - 8.0 * cosim + cosisq * (12.0 + 8.0 * cosim - 10.0 * cosisq));
                let xno2 = nm * nm;
                let ainv2 = aonv * aonv;
                let mut temp1 = 3.0 * xno2 * ainv2;
                let mut temp = temp1 * root22;
                out.d2201 = temp * f220 * g201;
                out.d2211 = temp * f221 * g211;
                temp1 *= aonv;
                temp = temp1 * root32;
                out.d3210 = temp * f321 * g310;
                out.d3222 = temp * f322 * g322;
                temp1 *= aonv;
                temp = 2.0 * temp1 * root44;
                out.d4410 = temp * f441 * g410;
                out.d4422 = temp * f442 * g422;
                temp1 *= aonv;
                temp = temp1 * root52;
                out.d5220 = temp * f522 * g520;
                out.d5232 = temp * f523 * g532;
                temp = 2.0 * temp1 * root54;
                out.d5421 = temp * f542 * g521;
                out.d5433 = temp * f543 * g533;
                out.xlamo = mod_floor_twopi(self.mo + self.nodeo + self.nodeo - theta - theta);
                out.xfact = self.mdot + out.dmdt + 2.0 * (self.nodedot + out.dnodt - rptim)
                    - self.no_unkozai;
                emsq = d.emsq;
            }
            if out.irez == 1 {
                // synchronous resonance terms
                let g200 = 1.0 + emsq * (-2.5 + 0.8125 * emsq);
                let g310 = 1.0 + 2.0 * emsq;
                let g300 = 1.0 + emsq * (-6.0 + 6.609_37 * emsq);
                let f220 = 0.75 * (1.0 + cosim) * (1.0 + cosim);
                let f311 = 0.9375 * sinim * sinim * (1.0 + 3.0 * cosim) - 0.75 * (1.0 + cosim);
                let mut f330 = 1.0 + cosim;
                f330 = 1.875 * f330 * f330 * f330;
                let mut del1 = 3.0 * nm * nm * aonv * aonv;
                out.del2 = 2.0 * del1 * f220 * g200 * q22;
                out.del3 = 3.0 * del1 * f330 * g300 * q33 * aonv;
                del1 = del1 * f311 * g310 * q31 * aonv;
                out.del1 = del1;
                out.xlamo = mod_floor_twopi(self.mo + self.nodeo + self.argpo - theta);
                out.xfact =
                    self.mdot + xpidot - rptim + out.dmdt + out.domdt + out.dnodt - self.no_unkozai;
            }
        }
        let _ = emsq;
        out
    }

    /// `dspace`: the deep-space secular effects and resonance integration at
    /// `t`, from the epoch. Returns `(em, argpm, inclm, mm, nodem, nm)`.
    #[allow(clippy::too_many_arguments, clippy::similar_names)]
    fn dspace(
        &self,
        t: f64,
        mut em: f64,
        mut argpm: f64,
        mut inclm: f64,
        mut mm: f64,
        mut nodem: f64,
        mut nm: f64,
    ) -> (f64, f64, f64, f64, f64, f64) {
        let ds = &self.secular;
        let fasx2 = 0.131_309_08;
        let fasx4 = 2.884_319_8;
        let fasx6 = 0.374_480_87;
        let g22 = 5.768_639_6;
        let g32 = 0.952_408_98;
        let g44 = 1.801_499_8;
        let g52 = 1.050_833_0;
        let g54 = 4.410_889_8;
        let rptim = 4.375_269_088_011_3e-3;
        let stepp = 720.0;
        let stepn = -720.0;
        let step2 = 259_200.0;

        let theta = mod_floor_twopi(self.gsto + t * rptim);
        em += ds.dedt * t;
        inclm += ds.didt * t;
        argpm += ds.domdt * t;
        nodem += ds.dnodt * t;
        mm += ds.dmdt * t;

        if ds.irez != 0 {
            // the integrator starts from the epoch on every call
            let mut atime = 0.0;
            let mut xni = self.no_unkozai;
            let mut xli = ds.xlamo;
            let delt = if t > 0.0 { stepp } else { stepn };
            let ft;
            let (mut xndt, mut xldot, mut xnddt);
            loop {
                if ds.irez != 2 {
                    xndt = ds.del1 * sin64(xli - fasx2)
                        + ds.del2 * sin64(2.0 * (xli - fasx4))
                        + ds.del3 * sin64(3.0 * (xli - fasx6));
                    xldot = xni + ds.xfact;
                    xnddt = ds.del1 * cos64(xli - fasx2)
                        + 2.0 * ds.del2 * cos64(2.0 * (xli - fasx4))
                        + 3.0 * ds.del3 * cos64(3.0 * (xli - fasx6));
                    xnddt *= xldot;
                } else {
                    let xomi = self.argpo + self.argpdot * atime;
                    let x2omi = xomi + xomi;
                    let x2li = xli + xli;
                    xndt = ds.d2201 * sin64(x2omi + xli - g22)
                        + ds.d2211 * sin64(xli - g22)
                        + ds.d3210 * sin64(xomi + xli - g32)
                        + ds.d3222 * sin64(-xomi + xli - g32)
                        + ds.d4410 * sin64(x2omi + x2li - g44)
                        + ds.d4422 * sin64(x2li - g44)
                        + ds.d5220 * sin64(xomi + xli - g52)
                        + ds.d5232 * sin64(-xomi + xli - g52)
                        + ds.d5421 * sin64(xomi + x2li - g54)
                        + ds.d5433 * sin64(-xomi + x2li - g54);
                    xldot = xni + ds.xfact;
                    xnddt = ds.d2201 * cos64(x2omi + xli - g22)
                        + ds.d2211 * cos64(xli - g22)
                        + ds.d3210 * cos64(xomi + xli - g32)
                        + ds.d3222 * cos64(-xomi + xli - g32)
                        + ds.d5220 * cos64(xomi + xli - g52)
                        + ds.d5232 * cos64(-xomi + xli - g52)
                        + 2.0
                            * (ds.d4410 * cos64(x2omi + x2li - g44)
                                + ds.d4422 * cos64(x2li - g44)
                                + ds.d5421 * cos64(xomi + x2li - g54)
                                + ds.d5433 * cos64(-xomi + x2li - g54));
                    xnddt *= xldot;
                }
                if (t - atime).abs() >= stepp {
                    xli = xli + xldot * delt + xndt * step2;
                    xni = xni + xndt * delt + xnddt * step2;
                    atime += delt;
                } else {
                    ft = t - atime;
                    break;
                }
            }
            nm = xni + xndt * ft + xnddt * ft * ft * 0.5;
            let xl = xli + xldot * ft + xndt * ft * ft * 0.5;
            if ds.irez != 1 {
                mm = xl - 2.0 * nodem + 2.0 * theta;
            } else {
                mm = xl - nodem - argpm + theta;
            }
            let dndt = nm - self.no_unkozai;
            nm = self.no_unkozai + dndt;
        }
        (em, argpm, inclm, mm, nodem, nm)
    }

    /// The TEME state `tsince` minutes after the epoch.
    ///
    /// # Errors
    ///
    /// [`Sgp4Error`] when the model has no state there (codes 1 to 4), and
    /// [`Sgp4Error::Decayed`] when the computed radius is below one Earth
    /// radius (code 6; the reference computes a state and flags it, here no
    /// state is returned).
    #[allow(clippy::too_many_lines, clippy::similar_names)]
    pub fn propagate(&self, tsince: f64) -> Result<TemeState, Sgp4Error> {
        let vkmpersec = self.radiusearthkm * self.xke / 60.0;
        let t = tsince;

        let xmdf = self.mo + self.mdot * t;
        let argpdf = self.argpo + self.argpdot * t;
        let nodedf = self.nodeo + self.nodedot * t;
        let mut argpm = argpdf;
        let mut mm = xmdf;
        let t2 = t * t;
        let mut nodem = nodedf + self.nodecf * t2;
        let mut tempa = 1.0 - self.cc1 * t;
        let mut tempe = self.bstar * self.cc4 * t;
        let mut templ = self.t2cof * t2;

        if !self.isimp {
            let delomg = self.omgcof * t;
            let delmtemp = 1.0 + self.eta * cos64(xmdf);
            let delm = self.xmcof * (delmtemp * delmtemp * delmtemp - self.delmo);
            let temp = delomg + delm;
            mm = xmdf + temp;
            argpm = argpdf - temp;
            let t3 = t2 * t;
            let t4 = t3 * t;
            tempa = tempa - self.d2 * t2 - self.d3 * t3 - self.d4 * t4;
            tempe += self.bstar * self.cc5 * (sin64(mm) - self.sinmao);
            templ = templ + self.t3cof * t3 + t4 * (self.t4cof + t * self.t5cof);
        }

        let mut nm = self.no_unkozai;
        let mut em = self.ecco;
        let mut inclm = self.inclo;
        if self.deep {
            let r = self.dspace(t, em, argpm, inclm, mm, nodem, nm);
            em = r.0;
            argpm = r.1;
            inclm = r.2;
            mm = r.3;
            nodem = r.4;
            nm = r.5;
        }

        // written so that NaN fails the check (the reference would go on with NaN)
        if nm.is_nan() || nm <= 0.0 {
            return Err(Sgp4Error::MeanMotion(nm));
        }
        let am = powf64(self.xke / nm, X2O3) * tempa * tempa;
        nm = self.xke / powf64(am, 1.5);
        em -= tempe;

        if !(-0.001..1.0).contains(&em) {
            return Err(Sgp4Error::MeanEccentricity(em));
        }
        if em < 1.0e-6 {
            em = 1.0e-6;
        }
        mm += self.no_unkozai * templ;
        let mut xlm = mm + argpm + nodem;
        nodem = fmod_twopi(nodem);
        argpm = mod_floor_twopi(argpm);
        xlm = mod_floor_twopi(xlm);
        mm = mod_floor_twopi(xlm - argpm - nodem);

        let sinim = sin64(inclm);
        let cosim = cos64(inclm);

        let mut ep = em;
        let mut xincp = inclm;
        let mut argpp = argpm;
        let mut nodep = nodem;
        let mut mp = mm;
        let mut sinip = sinim;
        let mut cosip = cosim;
        let mut aycof = self.aycof;
        let mut xlcof = self.xlcof;
        let mut con41 = self.con41;
        let mut x1mth2 = self.x1mth2;
        let mut x7thm1 = self.x7thm1;

        if self.deep {
            let r = self
                .periodics
                .apply(t, false, self.opsmode, ep, xincp, nodep, argpp, mp);
            ep = r.0;
            xincp = r.1;
            nodep = r.2;
            argpp = r.3;
            mp = r.4;
            if xincp < 0.0 {
                xincp = -xincp;
                nodep += PI;
                argpp -= PI;
            }
            if !(0.0..=1.0).contains(&ep) {
                return Err(Sgp4Error::PerturbedEccentricity(ep));
            }
            sinip = sin64(xincp);
            cosip = cos64(xincp);
            aycof = -0.5 * self.j3oj2 * sinip;
            xlcof = if (cosip + 1.0).abs() > 1.5e-12 {
                -0.25 * self.j3oj2 * sinip * (3.0 + 5.0 * cosip) / (1.0 + cosip)
            } else {
                -0.25 * self.j3oj2 * sinip * (3.0 + 5.0 * cosip) / TEMP4
            };
        }

        // long period periodics
        let axnl = ep * cos64(argpp);
        let mut temp = 1.0 / (am * (1.0 - ep * ep));
        let aynl = ep * sin64(argpp) + temp * aycof;
        let xl = mp + argpp + nodep + temp * xlcof * axnl;

        // Kepler's equation
        let u = mod_floor_twopi(xl - nodep);
        let mut eo1 = u;
        let mut tem5: f64 = 9999.9;
        let mut ktr = 1;
        let (mut sineo1, mut coseo1) = (0.0, 0.0);
        while tem5.abs() >= 1.0e-12 && ktr <= 10 {
            sineo1 = sin64(eo1);
            coseo1 = cos64(eo1);
            tem5 = 1.0 - coseo1 * axnl - sineo1 * aynl;
            tem5 = (u - aynl * coseo1 + axnl * sineo1 - eo1) / tem5;
            if tem5.abs() >= 0.95 {
                tem5 = if tem5 > 0.0 { 0.95 } else { -0.95 };
            }
            eo1 += tem5;
            ktr += 1;
        }

        // short period preliminary quantities
        let ecose = axnl * coseo1 + aynl * sineo1;
        let esine = axnl * sineo1 - aynl * coseo1;
        let el2 = axnl * axnl + aynl * aynl;
        let pl = am * (1.0 - el2);
        if pl.is_nan() || pl < 0.0 {
            return Err(Sgp4Error::SemiLatusRectum(pl));
        }
        let rl = am * (1.0 - ecose);
        let rdotl = sqrt64(am) * esine / rl;
        let rvdotl = sqrt64(pl) / rl;
        let betal = sqrt64(1.0 - el2);
        temp = esine / (1.0 + betal);
        let sinu = am / rl * (sineo1 - aynl - axnl * temp);
        let cosu = am / rl * (coseo1 - axnl + aynl * temp);
        let mut su = atan2_64(sinu, cosu);
        let sin2u = (cosu + cosu) * sinu;
        let cos2u = 1.0 - 2.0 * sinu * sinu;
        temp = 1.0 / pl;
        let temp1 = 0.5 * self.j2 * temp;
        let temp2 = temp1 * temp;

        if self.deep {
            let cosisq = cosip * cosip;
            con41 = 3.0 * cosisq - 1.0;
            x1mth2 = 1.0 - cosisq;
            x7thm1 = 7.0 * cosisq - 1.0;
        }
        let mrt = rl * (1.0 - 1.5 * temp2 * betal * con41) + 0.5 * temp1 * x1mth2 * cos2u;
        su -= 0.25 * temp2 * x7thm1 * sin2u;
        let xnode = nodep + 1.5 * temp2 * cosip * sin2u;
        let xinc = xincp + 1.5 * temp2 * cosip * sinip * cos2u;
        let mvt = rdotl - nm * temp1 * x1mth2 * sin2u / self.xke;
        let rvdot = rvdotl + nm * temp1 * (x1mth2 * cos2u + 1.5 * con41) / self.xke;

        // orientation vectors
        let sinsu = sin64(su);
        let cossu = cos64(su);
        let snod = sin64(xnode);
        let cnod = cos64(xnode);
        let sini = sin64(xinc);
        let cosi = cos64(xinc);
        let xmx = -snod * cosi;
        let xmy = cnod * cosi;
        let ux = xmx * sinsu + cnod * cossu;
        let uy = xmy * sinsu + snod * cossu;
        let uz = sini * sinsu;
        let vx = xmx * cossu - cnod * sinsu;
        let vy = xmy * cossu - snod * sinsu;
        let vz = sini * cossu;

        if mrt.is_nan() || mrt < 1.0 {
            return Err(Sgp4Error::Decayed(mrt));
        }
        let mr = mrt * self.radiusearthkm;
        Ok(TemeState {
            position_km: [mr * ux, mr * uy, mr * uz],
            velocity_km_s: [
                (mvt * ux + rvdot * vx) * vkmpersec,
                (mvt * uy + rvdot * vy) * vkmpersec,
                (mvt * uz + rvdot * vz) * vkmpersec,
            ],
        })
    }
}
