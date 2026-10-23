//! The Laws of this crate and their semantics identifiers.
//!
//! A Law is one evaluation rule whose output a stored value may come from. Its
//! `semantics_id` is the SHA-256 of its name, its data (constants as their bit
//! patterns), the alice-det-math `SEMANTICS_ID` where it evaluates `f64`
//! transcendentals, and the digest of its outputs over a frozen reference grid
//! of inputs (`tests/law_ids.rs` recomputes all of them).
//!
//! **The identifier moves if and only if the Law's inputs or its outputs on
//! the reference grid move.** The reference grids are pinned by their own
//! input hash, so they change only together with a change of meaning. A value
//! recorded with its identifier says which Law produced it.
//!
//! Author: Moroya Sakamoto

/// A Law: a name and a 32-byte semantics identifier.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct Law {
    /// Stable name, `alice-space/<law>`.
    pub name: &'static str,
    /// SHA-256 semantics identifier.
    pub semantics_id: [u8; 32],
}

/// The SGP4 / SDP4 propagator (Vallado et al. 2006 revision of the AFSPC
/// model) with the WGS-72 constants, in both operation modes.
pub const SGP4_AFSPC_WGS72_ID: [u8; 32] = [
    0xec, 0x42, 0x69, 0x26, 0x26, 0x15, 0xd9, 0x8f, 0xe2, 0xa1, 0x4d, 0x36, 0xd0, 0xc2, 0x5f, 0x75,
    0x2d, 0x60, 0x72, 0x8f, 0xca, 0x59, 0xd9, 0xfb, 0x69, 0x08, 0x6d, 0xdb, 0x7f, 0xc0, 0x23, 0xee,
];
/// [`SGP4_AFSPC_WGS72_ID`] with its name.
pub const SGP4_AFSPC_WGS72: Law = Law {
    name: "alice-space/sgp4-afspc-wgs72",
    semantics_id: SGP4_AFSPC_WGS72_ID,
};

/// Lagrange interpolation of tabulated trajectories (centred window,
/// relative times), as used for OEM segments and SP3 orbits.
pub const ODM_INTERP_LAGRANGE_ID: [u8; 32] = [
    0x2e, 0x04, 0x9c, 0xf3, 0xbe, 0x39, 0x5d, 0xe5, 0x70, 0xdd, 0x43, 0x9f, 0xe1, 0x3c, 0xea, 0xca,
    0x84, 0xee, 0xf6, 0x49, 0x2e, 0xb1, 0x9d, 0x07, 0x8b, 0x76, 0x90, 0x61, 0xc5, 0x66, 0x5e, 0xcc,
];
/// [`ODM_INTERP_LAGRANGE_ID`] with its name.
pub const ODM_INTERP_LAGRANGE: Law = Law {
    name: "alice-space/odm-interp-lagrange",
    semantics_id: ODM_INTERP_LAGRANGE_ID,
};

/// Hermite interpolation of tabulated positions with velocities (centred
/// window, relative times), as used for OEM segments.
pub const ODM_INTERP_HERMITE_ID: [u8; 32] = [
    0xdb, 0xdd, 0xb1, 0xa3, 0xd3, 0x8a, 0xd7, 0x6d, 0xd1, 0x75, 0xca, 0xf5, 0x74, 0xf5, 0xc6, 0x22,
    0x31, 0x0f, 0x0e, 0xb0, 0x8a, 0x00, 0xa0, 0xe4, 0x5e, 0xe7, 0x2c, 0x44, 0x0b, 0x7f, 0xfc, 0x67,
];
/// [`ODM_INTERP_HERMITE_ID`] with its name.
pub const ODM_INTERP_HERMITE: Law = Law {
    name: "alice-space/odm-interp-hermite",
    semantics_id: ODM_INTERP_HERMITE_ID,
};

/// The conical eclipse of the Sun by a spherical body, from positions.
pub const ECLIPSE_CONICAL_ID: [u8; 32] = [
    0x7d, 0x91, 0xae, 0x6c, 0x40, 0x72, 0x15, 0x2e, 0xe6, 0x09, 0x92, 0x9f, 0x72, 0xd6, 0x2b, 0xcc,
    0xa7, 0x4f, 0x15, 0xe7, 0x84, 0x74, 0x8e, 0x0c, 0xd2, 0xe9, 0xc3, 0xfd, 0x49, 0x64, 0x2b, 0x90,
];
/// [`ECLIPSE_CONICAL_ID`] with its name.
pub const ECLIPSE_CONICAL: Law = Law {
    name: "alice-space/eclipse-conical",
    semantics_id: ECLIPSE_CONICAL_ID,
};

/// The zonal gravity field: the acceleration and potential of
/// `gravity::ZonalGravity` (Bonnet's recursion, terms summed from J2 up).
pub const ZONAL_GRAVITY_ID: [u8; 32] = [
    0x32, 0x25, 0xb4, 0xde, 0xa0, 0xc0, 0xd7, 0x50, 0xaf, 0x46, 0x30, 0xc7, 0xb1, 0xe2, 0xe9, 0xd7,
    0xa6, 0x76, 0xd0, 0x3d, 0x0c, 0x72, 0x94, 0x89, 0xfe, 0xd5, 0xf1, 0xb6, 0x2e, 0x0f, 0x1d, 0x07,
];
/// [`ZONAL_GRAVITY_ID`] with its name.
pub const ZONAL_GRAVITY: Law = Law {
    name: "alice-space/zonal-gravity",
    semantics_id: ZONAL_GRAVITY_ID,
};

/// Every Law of the crate.
pub const ALL: [Law; 5] = [
    SGP4_AFSPC_WGS72,
    ODM_INTERP_LAGRANGE,
    ODM_INTERP_HERMITE,
    ECLIPSE_CONICAL,
    ZONAL_GRAVITY,
];
