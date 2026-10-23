//! The Laws' semantics identifiers (src/law.rs).
//!
//! Each identifier is the SHA-256 of the Law's name, its data, the
//! alice-det-math `SEMANTICS_ID` where it evaluates `f64` transcendentals, and
//! the digest of its outputs over a **frozen reference grid**. The rule: the
//! identifier moves if and only if the Law's inputs or its outputs on the
//! reference grid move. The grids' inputs are pinned by their own SHA-256 and
//! count (`reference_grids_are_frozen`), so editing a grid is a visible failure.
//!
//! All digests are the same bits on every target. Re-recording after an
//! intentional change of meaning, in this order: the `REFERENCE_*` digests from
//! the failure message, then the identifiers in src/law.rs from the next one.
//!
//! Author: Moroya Sakamoto

use alice_space::eclipse::{
    conical_illumination, Eclipse, EARTH_EQUATORIAL_RADIUS_KM, MOON_RADIUS_KM, SUN_RADIUS_KM,
};
use alice_space::interp;
use alice_space::law::{self, Law};
use alice_space::sgp4::{GravityModel, MeanElements, OpsMode, Sgp4};
use alice_space::tle;
use sha2::{Digest, Sha256};

const TLE_FILE: &str = include_str!("data/sgp4/SGP4-VER.TLE");

/// SHA-256 and count of the frozen reference inputs.
const REFERENCE_SGP4_INPUTS: (&str, usize) = (
    "4325242aab1a88f81cf291428f419b90cd1a56a12f544e7074c70507235d50a9",
    33,
);
const REFERENCE_INTERP_INPUTS: (&str, usize) = (
    "c961a51b59284a2dd0b246c71155a9467fdc7fe62654b55f76628180dd50bf40",
    65,
);
const REFERENCE_ECLIPSE_INPUTS: (&str, usize) = (
    "861ec7ff24119095b759d0addd620f776dbc325b368e0a65ce17de9c9cf19025",
    459,
);

/// Output digests over the reference grids (folded into the identifiers).
const REFERENCE_SGP4_AFSPC_WGS72: &str =
    "7b7fc3e229099a579e331bbf2168ddcd9b8491a82bd6617bc6b80d8d7ab26398";
const REFERENCE_ODM_INTERP_LAGRANGE: &str =
    "f39757dfd592d4cf3432e82a504ec2a0fa9b984630bb0e78af3877e93bf2ae09";
const REFERENCE_ODM_INTERP_HERMITE: &str =
    "a784541638c556fbb44843db44e33cd71bc10465aeb0ff7ccaa2034421cdd9a2";
const REFERENCE_ECLIPSE_CONICAL: &str =
    "79fb3bf85b47c6a1c2901395be9c8a7f5578aa8c9b5256fecac44ffc4903bc8b";

fn hex(b: &[u8]) -> String {
    b.iter().map(|x| format!("{x:02x}")).collect()
}

fn digest(bytes: &[u8]) -> String {
    hex(&Sha256::digest(bytes))
}

fn bits(out: &mut Vec<u8>, v: &[f64]) {
    for x in v {
        out.extend(x.to_bits().to_le_bytes());
    }
}

// ---- SGP4 ------------------------------------------------------------------

/// Minutes from the epoch at which every verification element set is
/// propagated (both sides of the epoch, a day, three days, a week).
const SGP4_TIMES: [f64; 10] = [
    -1440.0, -720.0, -1.0, 0.0, 1.0, 90.0, 360.0, 1440.0, 4320.0, 10080.0,
];

/// The verification element sets of Vallado et al. (2006): line 1 and line 2
/// cut to 68 columns with their check digits recomputed.
fn sgp4_sets() -> Vec<String> {
    let check = |line: &str| {
        let payload = &line[..68];
        let sum: u32 = payload
            .chars()
            .map(|c| c.to_digit(10).unwrap_or(u32::from(c == '-')))
            .sum();
        format!("{payload}{}", sum % 10)
    };
    let mut out = Vec::new();
    let mut lines = TLE_FILE.lines().map(|l| l.trim_end_matches('\r'));
    while let Some(l1) = lines.next() {
        if l1.starts_with("1 ") {
            let l2 = lines.next().unwrap();
            out.push(format!("{}\n{}", check(l1), check(l2)));
        }
    }
    out
}

fn sgp4_input_bytes(sets: &[String]) -> Vec<u8> {
    let mut out = Vec::new();
    for s in sets {
        out.extend(s.as_bytes());
        out.push(0);
    }
    bits(&mut out, &SGP4_TIMES);
    out
}

fn sgp4_outputs(sets: &[String]) -> Vec<u8> {
    let mut out = Vec::new();
    for s in sets {
        let elements = MeanElements::from_tle(&tle::parse(s).unwrap());
        for mode in [OpsMode::Afspc, OpsMode::Improved] {
            match Sgp4::new(&elements, GravityModel::Wgs72, mode) {
                Err(_) => out.push(0xff),
                Ok(sat) => {
                    out.push(0);
                    for &t in &SGP4_TIMES {
                        match sat.propagate(t) {
                            Ok(st) => {
                                out.push(0);
                                bits(&mut out, &st.position_km);
                                bits(&mut out, &st.velocity_km_s);
                            }
                            Err(e) => out.push(e.code()),
                        }
                    }
                }
            }
        }
    }
    out
}

// ---- interpolation ---------------------------------------------------------

/// A 41-node table, 60 s apart from a far epoch offset, with smooth
/// positions and velocities, and query times across it.
/// Node times, positions, velocities and query times.
type InterpTable = (Vec<f64>, Vec<[f64; 3]>, Vec<[f64; 3]>, Vec<f64>);

fn interp_table() -> InterpTable {
    let times: Vec<f64> = (0..41).map(|i| 60.0 * f64::from(i)).collect();
    let pos = |t: f64| {
        [
            7000.0 + 0.0005 * t * t,
            -1200.0 + 3.0 * t,
            50.0 * t / 60.0 - t * t * t / 1.0e7,
        ]
    };
    let vel = |t: f64| [0.001 * t, 3.0, 50.0 / 60.0 - 3.0 * t * t / 1.0e7];
    let values = times.iter().map(|&t| pos(t)).collect();
    let velocities = times.iter().map(|&t| vel(t)).collect();
    let queries = (0..=64).map(|k| 37.3 * f64::from(k)).collect();
    (times, values, velocities, queries)
}

fn interp_input_bytes() -> Vec<u8> {
    let (times, values, velocities, queries) = interp_table();
    let mut out = Vec::new();
    bits(&mut out, &times);
    for v in values.iter().chain(&velocities) {
        bits(&mut out, v);
    }
    bits(&mut out, &queries);
    out
}

fn state_bytes(out: &mut Vec<u8>, r: Result<interp::State, interp::InterpError>) {
    match r {
        Ok((p, v)) => {
            out.push(0);
            bits(out, &p);
            bits(out, &v);
        }
        Err(_) => out.push(0xff),
    }
}

fn lagrange_outputs() -> Vec<u8> {
    let (times, values, _, queries) = interp_table();
    let mut out = Vec::new();
    for degree in 1..=9 {
        for &t in &queries {
            state_bytes(&mut out, interp::lagrange(&times, &values, degree, t));
        }
    }
    out
}

fn hermite_outputs() -> Vec<u8> {
    let (times, values, velocities, queries) = interp_table();
    let mut out = Vec::new();
    for degree in [1, 3, 5, 7, 9, 11] {
        for &t in &queries {
            state_bytes(
                &mut out,
                interp::hermite(&times, &values, &velocities, degree, t),
            );
        }
    }
    out
}

// ---- eclipse ---------------------------------------------------------------

type EclipseCase = ([f64; 3], [f64; 3], [f64; 3], f64);

/// Satellites behind the Earth (Sun at 1 au on +x) from LEO to beyond the
/// umbra tip, and behind the Moon, offset across the shadow edges.
fn eclipse_cases() -> Vec<EclipseCase> {
    let sun = [149_597_870.7, 0.0, 0.0];
    let mut v = Vec::new();
    for dist in [6778.0, 7178.0, 26_560.0, 42_164.0, 400_000.0, 1_500_000.0] {
        for k in -24..=24 {
            let off = EARTH_EQUATORIAL_RADIUS_KM * f64::from(k) / 20.0;
            v.push((
                [-dist, off, off / 3.0],
                sun,
                [0.0; 3],
                EARTH_EQUATORIAL_RADIUS_KM,
            ));
        }
    }
    let moon = [384_400.0, 0.0, 0.0];
    for back in [5000.0, 50_000.0, 374_000.0, 380_000.0, 420_000.0] {
        for k in -16..=16 {
            let off = MOON_RADIUS_KM * f64::from(k) / 6.0;
            v.push(([moon[0] - back, off, -off / 2.0], sun, moon, MOON_RADIUS_KM));
        }
    }
    v
}

fn eclipse_input_bytes(cases: &[EclipseCase]) -> Vec<u8> {
    let mut out = Vec::new();
    for (sat, sun, body, r) in cases {
        bits(&mut out, sat);
        bits(&mut out, sun);
        bits(&mut out, body);
        bits(&mut out, &[*r]);
    }
    out
}

fn eclipse_outputs(cases: &[EclipseCase]) -> Vec<u8> {
    let mut out = Vec::new();
    for &(sat, sun, body, r) in cases {
        match conical_illumination(sat, sun, body, r) {
            Ok(i) => {
                out.push(match i.eclipse {
                    Eclipse::None => 0,
                    Eclipse::Partial => 1,
                    Eclipse::Annular => 2,
                    Eclipse::Total => 3,
                });
                bits(
                    &mut out,
                    &[
                        i.fraction,
                        i.sun_angular_radius,
                        i.body_angular_radius,
                        i.separation,
                    ],
                );
            }
            Err(_) => out.push(0xff),
        }
    }
    out
}

// ---- identifiers -----------------------------------------------------------

fn law_id(law: &Law, data: &[u8], grid: (&str, usize), reference: &str) -> [u8; 32] {
    let mut h = Sha256::new();
    h.update(b"alice-space law v1\0");
    h.update(law.name.as_bytes());
    h.update([0]);
    h.update((data.len() as u64).to_le_bytes());
    h.update(data);
    h.update(grid.0.as_bytes());
    h.update((grid.1 as u64).to_le_bytes());
    h.update(reference.as_bytes());
    h.finalize().into()
}

#[test]
fn reference_grids_are_frozen() {
    let sets = sgp4_sets();
    let (_, _, _, queries) = interp_table();
    let cases = eclipse_cases();
    let got = [
        (
            "SGP4_INPUTS",
            digest(&sgp4_input_bytes(&sets)),
            sets.len(),
            REFERENCE_SGP4_INPUTS,
        ),
        (
            "INTERP_INPUTS",
            digest(&interp_input_bytes()),
            queries.len(),
            REFERENCE_INTERP_INPUTS,
        ),
        (
            "ECLIPSE_INPUTS",
            digest(&eclipse_input_bytes(&cases)),
            cases.len(),
            REFERENCE_ECLIPSE_INPUTS,
        ),
    ];
    let mut bad = Vec::new();
    for (name, h, n, want) in got {
        assert!(n > 30, "{name}: {n}");
        if (h.as_str(), n) != want {
            bad.push(format!(
                "const REFERENCE_{name}: (&str, usize) = (\"{h}\", {n});"
            ));
        }
    }
    assert!(
        bad.is_empty(),
        "a reference grid changed (it is part of the Law identifiers):\n{}",
        bad.join("\n")
    );
}

#[test]
fn reference_digests() {
    let sets = sgp4_sets();
    let cases = eclipse_cases();
    let got = [
        (
            "SGP4_AFSPC_WGS72",
            sgp4_outputs(&sets),
            REFERENCE_SGP4_AFSPC_WGS72,
        ),
        (
            "ODM_INTERP_LAGRANGE",
            lagrange_outputs(),
            REFERENCE_ODM_INTERP_LAGRANGE,
        ),
        (
            "ODM_INTERP_HERMITE",
            hermite_outputs(),
            REFERENCE_ODM_INTERP_HERMITE,
        ),
        (
            "ECLIPSE_CONICAL",
            eclipse_outputs(&cases),
            REFERENCE_ECLIPSE_CONICAL,
        ),
    ];
    let mut bad = Vec::new();
    for (name, bytes, want) in &got {
        assert!(bytes.len() > 1000, "{name}: {} bytes", bytes.len());
        let h = digest(bytes);
        if h != *want {
            bad.push(format!("const REFERENCE_{name}: &str = \"{h}\";"));
        }
    }
    assert!(
        bad.is_empty(),
        "reference digests changed:\n{}",
        bad.join("\n")
    );
}

#[test]
fn semantics_ids() {
    // WGS-72 (mu, Earth radius, J2, J3, J4) as the propagator uses them
    let mut sgp4 = Vec::new();
    bits(
        &mut sgp4,
        &[
            398_600.8,
            6378.135,
            0.001_082_616,
            -2.538_81e-6,
            -1.655_97e-6,
        ],
    );
    sgp4.extend(alice_det_math::SEMANTICS_ID);
    // the window rule (no transcendental: + − × ÷ only)
    let window = b"centred window: even counts split evenly, odd counts on the nearest node; times relative to the first node".to_vec();
    let mut eclipse = Vec::new();
    bits(&mut eclipse, &[SUN_RADIUS_KM]);
    eclipse.extend(alice_det_math::SEMANTICS_ID);
    let cases = [
        (
            law::SGP4_AFSPC_WGS72,
            sgp4,
            REFERENCE_SGP4_INPUTS,
            REFERENCE_SGP4_AFSPC_WGS72,
        ),
        (
            law::ODM_INTERP_LAGRANGE,
            window.clone(),
            REFERENCE_INTERP_INPUTS,
            REFERENCE_ODM_INTERP_LAGRANGE,
        ),
        (
            law::ODM_INTERP_HERMITE,
            window,
            REFERENCE_INTERP_INPUTS,
            REFERENCE_ODM_INTERP_HERMITE,
        ),
        (
            law::ECLIPSE_CONICAL,
            eclipse,
            REFERENCE_ECLIPSE_INPUTS,
            REFERENCE_ECLIPSE_CONICAL,
        ),
    ];
    let mut bad = Vec::new();
    for (l, data, grid, reference) in &cases {
        let id = law_id(l, data, *grid, reference);
        if id != l.semantics_id {
            let bytes: Vec<String> = id.iter().map(|b| format!("0x{b:02x}")).collect();
            bad.push(format!("{}: [{}]", l.name, bytes.join(", ")));
        }
    }
    assert!(bad.is_empty(), "semantics ids changed:\n{}", bad.join("\n"));
    let mut ids: Vec<[u8; 32]> = law::ALL.iter().map(|l| l.semantics_id).collect();
    ids.sort_unstable();
    ids.dedup();
    assert_eq!(ids.len(), law::ALL.len(), "two Laws share an identifier");
}
