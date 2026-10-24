//! The Laws as v2 law references (src/law.rs, alice-zip `law_ref::LawRef`).
//!
//! `tests/data/law/v2_rows.txt` holds each Law's kind, semantics and rows, and the identifier and
//! bytes that `scripts/law_ref_check.py` (an independent standard-library implementation)
//! computes from them. These tests check that alice-zip gives the same identifier and bytes, that
//! the Laws of src/law.rs carry those identifiers, and that the identifier moves with the
//! semantics and with every reference input but not with the order the rows are written in.
//!
//! Author: Moroya Sakamoto

use alice_space::law;
use alice_zip::law_ref::{LawRef, Value};

const ROWS: &str = include_str!("data/law/v2_rows.txt");

/// One Law of the file.
struct Spec {
    kind: String,
    semantics: [u8; 32],
    rows: Vec<(String, Value)>,
    id: [u8; 32],
    bytes: Vec<u8>,
}

fn unhex(s: &str) -> Vec<u8> {
    assert!(s.len() % 2 == 0, "odd hex length: {s}");
    (0..s.len())
        .step_by(2)
        .map(|i| u8::from_str_radix(&s[i..i + 2], 16).expect("hex"))
        .collect()
}

fn id32(s: &str) -> [u8; 32] {
    unhex(s).try_into().expect("32 bytes")
}

fn value(ty: &str, text: &str) -> Value {
    match ty {
        "f64" => Value::F64(f64::from_bits(u64::from_be_bytes(
            unhex(text).try_into().expect("8 bytes"),
        ))),
        "u64" => Value::U64(text.parse().expect("u64")),
        "id" => Value::Id(id32(text)),
        "str" => Value::Str(String::from_utf8(unhex(text)).expect("UTF-8")),
        other => panic!("unknown value type {other}"),
    }
}

fn specs() -> Vec<Spec> {
    let mut out = Vec::new();
    let mut cur: Option<Spec> = None;
    for line in ROWS.lines() {
        let (word, rest) = line.split_once(' ').unwrap_or((line, ""));
        match word {
            "law" => {
                cur = Some(Spec {
                    kind: rest.to_owned(),
                    semantics: [0; 32],
                    rows: Vec::new(),
                    id: [0; 32],
                    bytes: Vec::new(),
                });
            }
            "semantics" => cur.as_mut().expect("law first").semantics = id32(rest),
            "row" => {
                let mut f = rest.splitn(3, ' ');
                let (name, ty, text) = (
                    f.next().expect("name"),
                    f.next().expect("type"),
                    f.next().expect("value"),
                );
                cur.as_mut()
                    .expect("law first")
                    .rows
                    .push((name.to_owned(), value(ty, text)));
            }
            "expect-id" => cur.as_mut().expect("law first").id = id32(rest),
            "expect-ref" => {
                let mut s = cur.take().expect("law first");
                s.bytes = unhex(rest);
                out.push(s);
            }
            _ => {}
        }
    }
    out
}

fn reference(s: &Spec) -> LawRef {
    LawRef::new(s.kind.as_bytes(), s.semantics, s.rows.clone()).expect("a valid law reference")
}

/// oracle: scripts/law_ref_check.py (standard-library Python, not alice-zip) computed the
/// identifier and the bytes of every Law in the file; alice-zip gives the same, and reads its
/// own bytes back to the same reference.
#[test]
fn the_independent_reference_and_alice_zip_agree() {
    let specs = specs();
    assert_eq!(specs.len(), law::ALL.len(), "one entry per Law");
    for s in &specs {
        let r = reference(s);
        assert_eq!(r.law_id(), &s.id, "{}: law id", s.kind);
        assert_eq!(r.to_bytes(), s.bytes, "{}: bytes", s.kind);
        assert_eq!(LawRef::from_bytes(&s.bytes).expect("reads back"), r, "{}", s.kind);
    }
}

/// The Laws of src/law.rs carry the v2 identifiers of the file, in the same order.
#[test]
fn the_laws_carry_their_v2_identifiers() {
    let specs = specs();
    assert_eq!(specs.len(), law::ALL.len());
    for (l, s) in law::ALL.iter().zip(&specs) {
        assert_eq!(l.name, s.kind);
        assert_eq!(l.semantics_id, s.id, "{}: semantics_id is not the v2 law id", l.name);
    }
}

/// The identifier moves when only the semantics moves.
#[test]
fn the_identifier_moves_with_the_semantics() {
    for s in &specs() {
        let mut sem = s.semantics;
        sem[31] ^= 1;
        let moved = LawRef::new(s.kind.as_bytes(), sem, s.rows.clone()).expect("valid");
        assert_ne!(moved.law_id(), &s.id, "{}", s.kind);
    }
}

/// The identifier moves when any reference row moves: the input digest, the input count or the
/// output digest (one bit or one count each).
#[test]
fn the_identifier_moves_with_every_reference_row() {
    for s in &specs() {
        let mut compared = 0;
        for (i, (name, v)) in s.rows.iter().enumerate() {
            let moved = match (name.as_str(), v) {
                ("reference_inputs" | "reference_outputs", Value::Id(id)) => {
                    let mut id = *id;
                    id[0] ^= 1;
                    Value::Id(id)
                }
                ("reference_input_count", Value::U64(n)) => Value::U64(n + 1),
                _ => continue,
            };
            let mut rows = s.rows.clone();
            rows[i].1 = moved;
            let r = LawRef::new(s.kind.as_bytes(), s.semantics, rows).expect("valid");
            assert_ne!(r.law_id(), &s.id, "{}: {name}", s.kind);
            compared += 1;
        }
        assert_eq!(compared, 3, "{}: three reference rows", s.kind);
    }
}

/// The order the rows are written in does not move the identifier (the rows are sorted).
#[test]
fn the_row_order_does_not_move_the_identifier() {
    for s in &specs() {
        let mut rows = s.rows.clone();
        rows.reverse();
        let r = LawRef::new(s.kind.as_bytes(), s.semantics, rows).expect("valid");
        assert_eq!(r.law_id(), &s.id, "{}", s.kind);
    }
}

/// Every float row is finite, so no Law can carry a NaN (refused by alice-zip, and by the
/// Python reference).
#[test]
fn every_float_row_is_finite() {
    let mut floats = 0;
    for s in &specs() {
        for (name, v) in &s.rows {
            if let Value::F64(x) = v {
                assert!(x.is_finite(), "{}: {name} = {x}", s.kind);
                floats += 1;
            }
        }
    }
    assert!(floats > 0, "the file has float rows");
}

/// The two interpolation Laws are evaluated with + − × ÷ only; their semantics is the named
/// constant, not alice-det-math's, and the others are alice-det-math's.
#[test]
fn the_semantics_of_each_law() {
    let basic: [u8; 32] = {
        use sha2::{Digest, Sha256};
        Sha256::digest(b"alice-space/semantics/ieee754-basic-ops/v1").into()
    };
    for s in &specs() {
        let want = if s.kind.starts_with("alice-space/odm-interp-") {
            basic
        } else {
            alice_det_math::SEMANTICS_ID
        };
        assert_eq!(s.semantics, want, "{}", s.kind);
    }
}
