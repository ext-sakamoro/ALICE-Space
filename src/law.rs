//! The Laws of this crate, their v2 law identifiers and their law references.
//!
//! A Law is one evaluation rule whose output a stored value may come from. Each Law is a v2 law
//! reference (`alice_zip::law_ref::LawRef`): its kind is the Law's name, its semantics is the
//! semantics it is evaluated under (alice-det-math's `SEMANTICS_ID` where it evaluates `f64`
//! transcendentals, [`BASIC_OPS_SEMANTICS_ID`] where it uses + − × ÷ only), and its rows are its
//! constants (by their bits), its rule, the Laws it depends on, and the digest and count of a
//! frozen reference grid of inputs with the digest of its outputs over that grid. Its
//! `semantics_id` is the v2 law identifier of those (`alice_zip::law_id::LawIdHasher`), so any
//! crate that holds the reference ([`law_ref`]) recomputes the same identifier.
//!
//! **The identifier moves if and only if the Law's semantics, its constants or rule, or its
//! inputs or outputs on the reference grid move.** `tests/law_ids.rs` recomputes the grids'
//! digests from the code, `tests/law_ref_v2.rs` checks the references against
//! `tests/data/law/v2_rows.txt`, and `scripts/law_ref_check.py` recomputes the identifiers
//! from that file with the standard library only. The identifiers before v2 are kept in
//! [`LEGACY_V1_IDS`].
//!
//! Author: Moroya Sakamoto

use alice_zip::law_ref::{LawRef, Value};

/// A Law: a name and a 32-byte semantics identifier (its v2 law identifier).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Law {
    /// Stable name, `alice-space/<law>`; the kind of its law reference.
    pub name: &'static str,
    /// The v2 law identifier.
    pub semantics_id: [u8; 32],
}

const fn hex32(s: &str) -> [u8; 32] {
    const fn nibble(c: u8) -> u8 {
        match c {
            b'0'..=b'9' => c - b'0',
            b'a'..=b'f' => c - b'a' + 10,
            _ => panic!("not lower-case hex"),
        }
    }
    let b = s.as_bytes();
    assert!(b.len() == 64, "32 bytes of hex");
    let mut out = [0u8; 32];
    let mut i = 0;
    while i < 32 {
        out[i] = (nibble(b[2 * i]) << 4) | nibble(b[2 * i + 1]);
        i += 1;
    }
    out
}

/// The semantics of a Law evaluated with the IEEE 754 binary64 operations + − × ÷ only (no
/// transcendental, so not alice-det-math's): SHA-256 of
/// `"alice-space/semantics/ieee754-basic-ops/v1"`.
pub const BASIC_OPS_SEMANTICS_ID: [u8; 32] =
    hex32("7765c4eee3f2248f8755040b2746a74741e608e3245dbcd63d3f4415ac292080");

/// The SGP4 / SDP4 propagator (Vallado et al. 2006 revision of the AFSPC
/// model) with the WGS-72 constants, in both operation modes.
pub const SGP4_AFSPC_WGS72_ID: [u8; 32] =
    hex32("353b40dc831c6313db05c93fcf4b9947a516799612caf514c9d2021defa41834");
/// [`SGP4_AFSPC_WGS72_ID`] with its name.
pub const SGP4_AFSPC_WGS72: Law = Law {
    name: "alice-space/sgp4-afspc-wgs72",
    semantics_id: SGP4_AFSPC_WGS72_ID,
};

/// Lagrange interpolation of tabulated trajectories (centred window,
/// relative times), as used for OEM segments and SP3 orbits.
pub const ODM_INTERP_LAGRANGE_ID: [u8; 32] =
    hex32("6e82199cf10b8c0a9601ea6aff46a8f7587026fb7ef8c08c7cfc85c362ae66a2");
/// [`ODM_INTERP_LAGRANGE_ID`] with its name.
pub const ODM_INTERP_LAGRANGE: Law = Law {
    name: "alice-space/odm-interp-lagrange",
    semantics_id: ODM_INTERP_LAGRANGE_ID,
};

/// Hermite interpolation of tabulated positions with velocities (centred
/// window, relative times), as used for OEM segments.
pub const ODM_INTERP_HERMITE_ID: [u8; 32] =
    hex32("220c655452c843551e9cdfa0ba349a281cd0070cf04a621dd1559c730bdf33bc");
/// [`ODM_INTERP_HERMITE_ID`] with its name.
pub const ODM_INTERP_HERMITE: Law = Law {
    name: "alice-space/odm-interp-hermite",
    semantics_id: ODM_INTERP_HERMITE_ID,
};

/// The conical eclipse of the Sun by a spherical body, from positions.
pub const ECLIPSE_CONICAL_ID: [u8; 32] =
    hex32("44bb4ddc6cf7cb17605756e44986e8361bf9b4bd9a15ac681363cf742bb8c1cd");
/// [`ECLIPSE_CONICAL_ID`] with its name.
pub const ECLIPSE_CONICAL: Law = Law {
    name: "alice-space/eclipse-conical",
    semantics_id: ECLIPSE_CONICAL_ID,
};

/// The zonal gravity field: the acceleration and potential of
/// `gravity::ZonalGravity` (Bonnet's recursion, terms summed from J2 up).
pub const ZONAL_GRAVITY_ID: [u8; 32] =
    hex32("8fe6a972775cac678cf59e8fd94c7af6361a2fbcdfb78f7c0fafd4504a0973c0");
/// [`ZONAL_GRAVITY_ID`] with its name.
pub const ZONAL_GRAVITY: Law = Law {
    name: "alice-space/zonal-gravity",
    semantics_id: ZONAL_GRAVITY_ID,
};

/// The third-body perturbation `gravity::third_body` (Battin's f(q)).
pub const THIRD_BODY_ID: [u8; 32] =
    hex32("2da76d1e57048bf6d1466dfd9fe20fda19bbe7241992e49e6ca0b55b43f07630");
/// [`THIRD_BODY_ID`] with its name.
pub const THIRD_BODY: Law = Law {
    name: "alice-space/third-body-battin",
    semantics_id: THIRD_BODY_ID,
};

/// Solar radiation pressure `srp::solar_radiation_pressure`: the cannonball
/// model scaled by the conical shadow's visible fraction of the Sun.
pub const SOLAR_RADIATION_PRESSURE_ID: [u8; 32] =
    hex32("d890f70de121079472fcfda71e18e4a9f042c6a11e9c5ffd6dfe911ed2bb8849");
/// [`SOLAR_RADIATION_PRESSURE_ID`] with its name.
pub const SOLAR_RADIATION_PRESSURE: Law = Law {
    name: "alice-space/srp-cannonball-conical",
    semantics_id: SOLAR_RADIATION_PRESSURE_ID,
};

/// Atmospheric drag `drag::atmospheric_drag` relative to a co-rotating
/// atmosphere, with the caller's density.
pub const DRAG_ID: [u8; 32] =
    hex32("5324b239e9f7f3319cab6001a0f15f0f8bc323d58c9a9c9cd6ac71d0ff80067b");
/// [`DRAG_ID`] with its name.
pub const DRAG: Law = Law {
    name: "alice-space/drag-corotating",
    semantics_id: DRAG_ID,
};

/// Every Law of the crate.
pub const ALL: [Law; 8] = [
    SGP4_AFSPC_WGS72,
    ODM_INTERP_LAGRANGE,
    ODM_INTERP_HERMITE,
    ECLIPSE_CONICAL,
    ZONAL_GRAVITY,
    THIRD_BODY,
    SOLAR_RADIATION_PRESSURE,
    DRAG,
];

/// The identifiers before v2 (SHA-256 of `"alice-space law v1\0"`, the name, the data, the grid
/// digest and count and the output digest), by Law name. A value recorded with one of these was
/// produced by the same evaluation as the Law of that name; `tests/law_ids.rs` recomputes them.
pub const LEGACY_V1_IDS: [(&str, [u8; 32]); 8] = [
    (
        "alice-space/sgp4-afspc-wgs72",
        hex32("ec4269262615d98fe2a14d36d0c25f752d60728fca59d9fb69086ddb7fc023ee"),
    ),
    (
        "alice-space/odm-interp-lagrange",
        hex32("2e049cf3be395de570dd439fe13ceaca84eef6492eb19d078b769061c5665ecc"),
    ),
    (
        "alice-space/odm-interp-hermite",
        hex32("dbddb1a3d38ad76dd175caf574f5c622310f0eb08a00a0e45ee72c440b7ffc67"),
    ),
    (
        "alice-space/eclipse-conical",
        hex32("7d91ae6c4072152ee609929f72d62bcca74f15e784748e0cd2e9c3fd49642b90"),
    ),
    (
        "alice-space/zonal-gravity",
        hex32("3225b4dea0c0d750af4630c7b1e2e9d7a676d03d0c729489fed5f1b62e0f1d07"),
    ),
    (
        "alice-space/third-body-battin",
        hex32("22dba1cb06c0f8a738004d3fd86d9a618190d30553a106c152cd48dc5199a8d1"),
    ),
    (
        "alice-space/srp-cannonball-conical",
        hex32("21f54f94edace7d6c3faa55cd7eb468bd044d4dc31cf74079e857de04029489b"),
    ),
    (
        "alice-space/drag-corotating",
        hex32("8aa774af42c1463a5cc7bd9a769194bb3383cc7a57f97c609047e952ce36f81f"),
    ),
];

/// The frozen reference grid of a Law: the SHA-256 and count of its inputs and the SHA-256 of
/// its outputs over them (`tests/law_ids.rs` recomputes all three from the code).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct ReferenceGrid {
    /// SHA-256 of the grid's inputs.
    pub inputs: [u8; 32],
    /// The number of inputs.
    pub input_count: u64,
    /// SHA-256 of the Law's outputs over the inputs.
    pub outputs: [u8; 32],
}

const SGP4_INPUTS: [u8; 32] =
    hex32("4325242aab1a88f81cf291428f419b90cd1a56a12f544e7074c70507235d50a9");
const INTERP_INPUTS: [u8; 32] =
    hex32("c961a51b59284a2dd0b246c71155a9467fdc7fe62654b55f76628180dd50bf40");
const THIRD_INPUTS: [u8; 32] =
    hex32("a0de1e1d27ad2de07ade34495fd9e480a0346d7939b11fb36ef49afb80a7c57c");
const SRP_INPUTS: [u8; 32] =
    hex32("6e4a193fffce69850f1e788ffe06f0ad67efba84b0047577fe50c98cd46efcba");
const DRAG_INPUTS: [u8; 32] =
    hex32("ca686d751f7bbfab54fcac0055aaaa35f9193b9d83490fec17426eea9c2d37d9");
const ZONAL_INPUTS: [u8; 32] =
    hex32("be322080b9a1056ade3b4d56c494049f661d6b157b5c4c43195d5f3692e37e34");
const ECLIPSE_INPUTS: [u8; 32] =
    hex32("861ec7ff24119095b759d0addd620f776dbc325b368e0a65ce17de9c9cf19025");

/// The reference grid of `law`, or `None` for a Law that is not one of [`ALL`].
#[must_use]
pub fn reference_grid(law: &Law) -> Option<ReferenceGrid> {
    let (inputs, input_count, outputs) = match law.name {
        "alice-space/sgp4-afspc-wgs72" => (
            SGP4_INPUTS,
            33,
            "7b7fc3e229099a579e331bbf2168ddcd9b8491a82bd6617bc6b80d8d7ab26398",
        ),
        "alice-space/odm-interp-lagrange" => (
            INTERP_INPUTS,
            65,
            "f39757dfd592d4cf3432e82a504ec2a0fa9b984630bb0e78af3877e93bf2ae09",
        ),
        "alice-space/odm-interp-hermite" => (
            INTERP_INPUTS,
            65,
            "a784541638c556fbb44843db44e33cd71bc10465aeb0ff7ccaa2034421cdd9a2",
        ),
        "alice-space/eclipse-conical" => (
            ECLIPSE_INPUTS,
            459,
            "79fb3bf85b47c6a1c2901395be9c8a7f5578aa8c9b5256fecac44ffc4903bc8b",
        ),
        "alice-space/zonal-gravity" => (
            ZONAL_INPUTS,
            121,
            "eca9b9a0460ae56e5032586e19bfd245382c0cbc1305c444f95013dcb0b8ff0b",
        ),
        "alice-space/third-body-battin" => (
            THIRD_INPUTS,
            98,
            "ba8a14c951f28c516b40981a4f456494924c302542e7b1a76683b059acca8041",
        ),
        "alice-space/srp-cannonball-conical" => (
            SRP_INPUTS,
            462,
            "37ee5039d3a8b685596de241667839e6ac188f8f6dfb8c6858aacb9b034ff248",
        ),
        "alice-space/drag-corotating" => (
            DRAG_INPUTS,
            120,
            "66055564a695ea39a0dda533ad3b011c72fcd03bdb6ce0568c1dc7ddb7f5f6d6",
        ),
        _ => return None,
    };
    Some(ReferenceGrid {
        inputs,
        input_count,
        outputs: hex32(outputs),
    })
}

const WINDOW_RULE: &str = "centred window: even counts split evenly, odd counts on the nearest node; times relative to the first node";

/// The reference of `law`: its kind, semantics and rows, with its v2 law identifier, or `None`
/// for a Law that is not one of [`ALL`]. The identifier of the result is `law.semantics_id`.
#[must_use]
pub fn law_ref(law: &Law) -> Option<LawRef> {
    let grid = reference_grid(law)?;
    let f = |name: &str, v: f64| (String::from(name), Value::F64(v));
    let rule = |text: &str| (String::from("rule"), Value::Str(String::from(text)));
    let (semantics, mut rows): ([u8; 32], Vec<(String, Value)>) = match law.name {
        "alice-space/sgp4-afspc-wgs72" => (
            alice_det_math::SEMANTICS_ID,
            std::vec![
                f("mu_km3_s2", 398_600.8),
                f("earth_radius_km", 6378.135),
                f("j2", 0.001_082_616),
                f("j3", -2.538_81e-6),
                f("j4", -1.655_97e-6),
            ],
        ),
        "alice-space/odm-interp-lagrange" | "alice-space/odm-interp-hermite" => {
            (BASIC_OPS_SEMANTICS_ID, std::vec![rule(WINDOW_RULE)])
        }
        "alice-space/eclipse-conical" => (
            alice_det_math::SEMANTICS_ID,
            std::vec![f("sun_radius_km", crate::eclipse::SUN_RADIUS_KM)],
        ),
        "alice-space/zonal-gravity" => (
            alice_det_math::SEMANTICS_ID,
            std::vec![rule(
                "Bonnet recursion; terms summed from J2 up; a = -(mu/r^2)(radial r^ + axial z^)"
            )],
        ),
        "alice-space/third-body-battin" => (
            alice_det_math::SEMANTICS_ID,
            std::vec![rule(
                "Battin: q = r.(r - 2s)/|s|^2, f(q) = q(3 + 3q + q^2)/(1 + (1+q)^(3/2)), a = -mu (r + f s)/|r - s|^3"
            )],
        ),
        "alice-space/srp-cannonball-conical" => (
            alice_det_math::SEMANTICS_ID,
            std::vec![
                f("solar_irradiance_w_m2", crate::srp::SOLAR_IRRADIANCE_W_M2),
                f("speed_of_light_m_s", crate::srp::SPEED_OF_LIGHT_M_S),
                f("au_km", crate::srp::AU_KM),
                (String::from("shadow_law"), Value::Id(ECLIPSE_CONICAL_ID)),
            ],
        ),
        "alice-space/drag-corotating" => (
            alice_det_math::SEMANTICS_ID,
            std::vec![rule(
                "a = -1/2 rho Cd (A/m) |v - w x r| (v - w x r), km and km/s"
            )],
        ),
        _ => return None,
    };
    rows.push((String::from("reference_inputs"), Value::Id(grid.inputs)));
    rows.push((
        String::from("reference_input_count"),
        Value::U64(grid.input_count),
    ));
    rows.push((String::from("reference_outputs"), Value::Id(grid.outputs)));
    // every row name is distinct and short and every float is finite, so the reference is valid
    Some(
        LawRef::new(law.name.as_bytes(), semantics, rows)
            .expect("a Law's rows form a valid law reference"),
    )
}

/// The references of [`ALL`], in that order.
#[must_use]
pub fn law_refs() -> Vec<LawRef> {
    ALL.iter()
        .map(|l| law_ref(l).expect("every Law of ALL has a reference"))
        .collect()
}
