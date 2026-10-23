//! Degenerate inputs of the SGP4 propagator: each case states the expected
//! result (an initialisation error, a propagation error code, or a finite state)
//! and asserts it. No case may return a state with a non-finite component.
//! Inputs for which the reference code divides by zero or takes the root of a
//! negative number at initialisation (e outside [0, 1), n <= 0) are rejected by
//! `new`.
//!
//! Author: Moroya Sakamoto

use alice_det_math::sqrt64;

use alice_space::sgp4::{GravityModel, MeanElements, OpsMode, Sgp4, Sgp4Error, Sgp4InitError};

/// A near-Earth set (ISS-like, 92 min, e = 0.0005) to vary one field at a time.
fn base() -> MeanElements {
    MeanElements {
        epoch_days_since_1950: 27_000.5,
        bstar: 3.0e-5,
        ndot: 0.0,
        nddot: 0.0,
        eccentricity: 5.0e-4,
        arg_perigee: 1.0,
        inclination: 0.9,
        mean_anomaly: 2.0,
        mean_motion_kozai: 2.0 * core::f64::consts::PI / 92.0,
        raan: 3.0,
    }
}

fn init(e: &MeanElements) -> Result<Sgp4, Sgp4InitError> {
    Sgp4::new(e, GravityModel::Wgs72, OpsMode::Improved)
}

fn propagate(e: &MeanElements, t: f64) -> Result<[f64; 6], Sgp4Error> {
    let s = init(e).expect("initialises").propagate(t)?;
    let r = s.position_km;
    let v = s.velocity_km_s;
    let out = [r[0], r[1], r[2], v[0], v[1], v[2]];
    assert!(
        out.iter().all(|x| x.is_finite()),
        "a returned state is finite: {out:?}"
    );
    Ok(out)
}

/// The base set propagates (the cases below differ from it in one field).
#[test]
fn the_base_set_propagates_to_a_low_earth_orbit_radius() {
    let s = propagate(&base(), 0.0).expect("propagates");
    let r = sqrt64(s[0] * s[0] + s[1] * s[1] + s[2] * s[2]);
    assert!((6600.0..6900.0).contains(&r), "radius {r} km");
}

/// e = 0 exactly: the model floors the mean eccentricity at 1e-6 and returns a
/// finite state.
#[test]
fn a_circular_orbit_propagates() {
    let e = MeanElements {
        eccentricity: 0.0,
        ..base()
    };
    for t in [0.0, 100.0, -100.0] {
        propagate(&e, t).expect("circular orbit propagates");
    }
}

/// e ≥ 1, e < 0: rejected at initialisation.
#[test]
fn an_eccentricity_outside_zero_to_one_is_rejected() {
    for ecc in [1.0, 1.5, -0.1] {
        let e = MeanElements {
            eccentricity: ecc,
            ..base()
        };
        assert_eq!(
            init(&e).err(),
            Some(Sgp4InitError::Eccentricity(ecc)),
            "e = {ecc}"
        );
    }
}

/// A zero or negative mean motion: rejected at initialisation.
#[test]
fn a_non_positive_mean_motion_is_rejected() {
    for n in [0.0, -0.01] {
        let e = MeanElements {
            mean_motion_kozai: n,
            ..base()
        };
        assert_eq!(
            init(&e).err(),
            Some(Sgp4InitError::MeanMotion(n)),
            "n = {n}"
        );
    }
}

/// A NaN or infinite field: rejected at initialisation, naming the field.
#[test]
fn a_non_finite_field_is_rejected() {
    let cases = [
        (
            MeanElements {
                eccentricity: f64::NAN,
                ..base()
            },
            "eccentricity",
        ),
        (
            MeanElements {
                mean_motion_kozai: f64::NAN,
                ..base()
            },
            "mean motion",
        ),
        (
            MeanElements {
                bstar: f64::INFINITY,
                ..base()
            },
            "bstar",
        ),
        (
            MeanElements {
                raan: f64::NAN,
                ..base()
            },
            "right ascension of the node",
        ),
    ];
    for (e, name) in cases {
        assert_eq!(init(&e).err(), Some(Sgp4InitError::NonFinite(name)));
    }
}

/// A NaN time: an error, not a NaN state.
#[test]
fn a_nan_time_is_an_error() {
    assert!(propagate(&base(), f64::NAN).is_err());
}

/// A huge drag term decays the orbit within days: error code 1 or 6 (the
/// eccentricity or the radius leaves its range first), never a state below the
/// surface.
#[test]
fn a_huge_drag_term_ends_in_decay() {
    let e = MeanElements {
        bstar: 0.5,
        ..base()
    };
    let mut ended = None;
    for k in 0..200 {
        if let Err(x) = propagate(&e, f64::from(k) * 60.0) {
            ended = Some(x.code());
            break;
        }
    }
    assert!(matches!(ended, Some(1 | 6)), "ended with {ended:?}");
}
