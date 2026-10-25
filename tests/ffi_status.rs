//! The C entry points that take buffers: a null pointer, a buffer of the wrong length, a non-finite
//! input or a non-finite result is a nonzero status, and on any nonzero status nothing is written
//! to the output buffer (it keeps the caller's values).
#![cfg(feature = "ffi")]

use alice_space::ffi::{
    alice_space_delta_v_hohmann, alice_space_rk4_single, alice_space_two_body_accel,
    ALICE_SPACE_ERR_LENGTH, ALICE_SPACE_ERR_NON_FINITE_INPUT, ALICE_SPACE_ERR_NON_FINITE_RESULT,
    ALICE_SPACE_ERR_NULL, ALICE_SPACE_OK,
};

const MU: f64 = 398_600.441_8;
/// a value no computation writes: an output left with it was not written
const UNTOUCHED: f64 = -12_345.0;

#[test]
fn the_statuses_are_distinct_and_ok_is_zero() {
    let all = [
        ALICE_SPACE_OK,
        ALICE_SPACE_ERR_NULL,
        ALICE_SPACE_ERR_LENGTH,
        ALICE_SPACE_ERR_NON_FINITE_INPUT,
        ALICE_SPACE_ERR_NON_FINITE_RESULT,
    ];
    assert_eq!(ALICE_SPACE_OK, 0);
    for (i, a) in all.iter().enumerate() {
        for b in &all[i + 1..] {
            assert_ne!(a, b);
        }
    }
}

#[test]
fn hohmann_writes_on_success_only() {
    let mut out = [UNTOUCHED; 2];
    let s = unsafe { alice_space_delta_v_hohmann(6578.0, 42_164.0, MU, out.as_mut_ptr(), 2) };
    assert_eq!(s, ALICE_SPACE_OK);
    assert!(out.iter().all(|v| v.is_finite() && *v > 0.0), "{out:?}");

    let s = unsafe { alice_space_delta_v_hohmann(6578.0, 42_164.0, MU, core::ptr::null_mut(), 2) };
    assert_eq!(s, ALICE_SPACE_ERR_NULL);

    let mut out = [UNTOUCHED; 2];
    let s = unsafe { alice_space_delta_v_hohmann(6578.0, 42_164.0, MU, out.as_mut_ptr(), 1) };
    assert_eq!(s, ALICE_SPACE_ERR_LENGTH);
    assert_eq!(out, [UNTOUCHED; 2]);

    for (r1, r2, mu) in [
        (f64::NAN, 42_164.0, MU),
        (6578.0, f64::INFINITY, MU),
        (6578.0, 42_164.0, f64::NAN),
    ] {
        let mut out = [UNTOUCHED; 2];
        let s = unsafe { alice_space_delta_v_hohmann(r1, r2, mu, out.as_mut_ptr(), 2) };
        assert_eq!(s, ALICE_SPACE_ERR_NON_FINITE_INPUT, "{r1} {r2} {mu}");
        assert_eq!(out, [UNTOUCHED; 2]);
    }
    // finite inputs whose result is not finite (a zero radius)
    let mut out = [UNTOUCHED; 2];
    let s = unsafe { alice_space_delta_v_hohmann(0.0, 42_164.0, MU, out.as_mut_ptr(), 2) };
    assert_eq!(s, ALICE_SPACE_ERR_NON_FINITE_RESULT);
    assert_eq!(out, [UNTOUCHED; 2]);
}

#[test]
fn rk4_reads_and_writes_six_values_on_success_only() {
    let state = [6778.0, 0.0, 0.0, 0.0, 7.67, 0.0];
    let mut out = [UNTOUCHED; 6];
    let s = unsafe { alice_space_rk4_single(state.as_ptr(), 6, MU, 60.0, out.as_mut_ptr(), 6) };
    assert_eq!(s, ALICE_SPACE_OK);
    assert!(out.iter().all(|v| v.is_finite() && *v != UNTOUCHED));

    let mut out = [UNTOUCHED; 6];
    let s = unsafe { alice_space_rk4_single(core::ptr::null(), 6, MU, 60.0, out.as_mut_ptr(), 6) };
    assert_eq!(s, ALICE_SPACE_ERR_NULL);
    assert_eq!(out, [UNTOUCHED; 6]);
    let s =
        unsafe { alice_space_rk4_single(state.as_ptr(), 6, MU, 60.0, core::ptr::null_mut(), 6) };
    assert_eq!(s, ALICE_SPACE_ERR_NULL);

    for (lin, lout) in [(5, 6), (6, 5), (0, 6), (7, 6)] {
        let mut out = [UNTOUCHED; 6];
        let s = unsafe {
            alice_space_rk4_single(state.as_ptr(), lin, MU, 60.0, out.as_mut_ptr(), lout)
        };
        assert_eq!(s, ALICE_SPACE_ERR_LENGTH, "{lin} {lout}");
        assert_eq!(out, [UNTOUCHED; 6]);
    }

    let mut bad = state;
    bad[4] = f64::NAN;
    let mut out = [UNTOUCHED; 6];
    let s = unsafe { alice_space_rk4_single(bad.as_ptr(), 6, MU, 60.0, out.as_mut_ptr(), 6) };
    assert_eq!(s, ALICE_SPACE_ERR_NON_FINITE_INPUT);
    assert_eq!(out, [UNTOUCHED; 6]);
    let s = unsafe {
        alice_space_rk4_single(state.as_ptr(), 6, MU, f64::INFINITY, out.as_mut_ptr(), 6)
    };
    assert_eq!(s, ALICE_SPACE_ERR_NON_FINITE_INPUT);
    assert_eq!(out, [UNTOUCHED; 6]);

    // finite inputs whose result overflows: μ = 1e308 at 1e-10 km (the acceleration is −∞).
    // The centre itself is not one: the two-body acceleration is 0 there by its own contract
    let near = [1e-10, 0.0, 0.0, 0.0, 0.0, 0.0];
    let s = unsafe { alice_space_rk4_single(near.as_ptr(), 6, 1e308, 60.0, out.as_mut_ptr(), 6) };
    assert_eq!(s, ALICE_SPACE_ERR_NON_FINITE_RESULT);
    assert_eq!(out, [UNTOUCHED; 6]);
}

#[test]
fn rk4_may_read_and_write_the_same_buffer() {
    let mut state = [6778.0, 0.0, 0.0, 0.0, 7.67, 0.0];
    let mut separate = [0.0; 6];
    let s =
        unsafe { alice_space_rk4_single(state.as_ptr(), 6, MU, 60.0, separate.as_mut_ptr(), 6) };
    assert_eq!(s, ALICE_SPACE_OK);
    let p = state.as_mut_ptr();
    let s = unsafe { alice_space_rk4_single(p, 6, MU, 60.0, p, 6) };
    assert_eq!(s, ALICE_SPACE_OK);
    assert_eq!(state.map(f64::to_bits), separate.map(f64::to_bits));
}

#[test]
fn two_body_writes_on_success_only() {
    let r = [6778.0, 0.0, 0.0];
    let mut out = [UNTOUCHED; 3];
    let s = unsafe { alice_space_two_body_accel(r.as_ptr(), 3, MU, out.as_mut_ptr(), 3) };
    assert_eq!(s, ALICE_SPACE_OK);
    assert!(out[0] < 0.0);

    let mut out = [UNTOUCHED; 3];
    for s in [
        unsafe { alice_space_two_body_accel(core::ptr::null(), 3, MU, out.as_mut_ptr(), 3) },
        unsafe { alice_space_two_body_accel(r.as_ptr(), 3, MU, core::ptr::null_mut(), 3) },
    ] {
        assert_eq!(s, ALICE_SPACE_ERR_NULL);
    }
    assert_eq!(out, [UNTOUCHED; 3]);
    let s = unsafe { alice_space_two_body_accel(r.as_ptr(), 2, MU, out.as_mut_ptr(), 3) };
    assert_eq!(s, ALICE_SPACE_ERR_LENGTH);
    let s = unsafe { alice_space_two_body_accel(r.as_ptr(), 3, f64::NAN, out.as_mut_ptr(), 3) };
    assert_eq!(s, ALICE_SPACE_ERR_NON_FINITE_INPUT);
    // μ = 1e308 at 1e-10 km overflows to −∞ (the centre itself gives 0, a finite result)
    let s = unsafe {
        alice_space_two_body_accel([1e-10, 0.0, 0.0].as_ptr(), 3, 1e308, out.as_mut_ptr(), 3)
    };
    assert_eq!(s, ALICE_SPACE_ERR_NON_FINITE_RESULT);
    assert_eq!(out, [UNTOUCHED; 3]);
    let s = unsafe { alice_space_two_body_accel([0.0; 3].as_ptr(), 3, MU, out.as_mut_ptr(), 3) };
    assert_eq!(s, ALICE_SPACE_OK);
    assert_eq!(out, [0.0; 3]);
}
