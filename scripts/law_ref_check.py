#!/usr/bin/env python3
"""An independent check of the Laws' v2 identifiers and law references (src/law.rs).

Each Law of this crate is a v2 law reference (alice-zip `law_ref::LawRef`): a kind (the Law's
name), the semantics it is evaluated under and named rows. Its identifier is the v2 law
identifier of the kind and the semantics with, for each row in ascending name order, the row's
name as a string value and then the row's value. This script computes both from the rows in
tests/data/law/v2_rows.txt with the standard library only; tests/law_ref_v2.rs reads the same
file and checks that src/law.rs gives the same identifiers and the same bytes.

The law identifier, every length a big-endian u64:

    len(b"alice-law-id/v2") b"alice-law-id/v2" | len(kind) kind | semantics id (32)
    then each value with a one-byte tag: 0x01 bytes (len, data) 0x02 u32 (4 BE) 0x03 u64 (8 BE)
    0x04 f64 (bits, 8 BE) 0x05 str (len, UTF-8) 0x06 id (32) 0x07 f32 (bits, 4 BE)

The law reference, every integer big-endian:

    "LOLR" | u16 2 | u8 len(kind) kind | semantics id (32) | law id (32) | u8 row count
    each row: u8 len(name) name | u8 tag | value (bytes and str: u64 len, then the data)

    python3 scripts/law_ref_check.py --check   # every expect line equals this computation
    python3 scripts/law_ref_check.py --write   # recompute the expect lines (a deliberate move)

Fails when the file holds no Law, so it never passes vacuously.

Author: Moroya Sakamoto
"""

from __future__ import annotations

import hashlib
import pathlib
import struct
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
ROWS = ROOT / "tests" / "data" / "law" / "v2_rows.txt"
DOMAIN = b"alice-law-id/v2"
TAGS = {"bytes": 0x01, "u32": 0x02, "u64": 0x03, "f64": 0x04, "str": 0x05, "id": 0x06, "f32": 0x07}
CANONICAL_NAN_F64 = 0x7FF8000000000000


def value_bytes(kind: str, raw: bytes) -> bytes:
    """The value as the law identifier and the law reference write it after the tag."""
    if kind in ("bytes", "str"):
        return struct.pack(">Q", len(raw)) + raw
    return raw


def parse_value(kind: str, text: str) -> bytes:
    if kind == "u64":
        return struct.pack(">Q", int(text))
    if kind == "u32":
        return struct.pack(">I", int(text))
    if kind == "f64":
        raw = bytes.fromhex(text)
        bits = struct.unpack(">Q", raw)[0]
        if (bits & 0x7FF0000000000000) == 0x7FF0000000000000 and bits & 0x000FFFFFFFFFFFFF:
            if bits != CANONICAL_NAN_F64:
                raise ValueError(f"non-canonical NaN {text}")
        return raw
    if kind in ("id",):
        raw = bytes.fromhex(text)
        if len(raw) != 32:
            raise ValueError(f"an id is 32 bytes, not {len(raw)}")
        return raw
    if kind in ("str", "bytes"):
        return bytes.fromhex(text)
    raise ValueError(f"unknown value type {kind}")


def law_id(kind: bytes, semantics: bytes, rows: list[tuple[str, str, bytes]]) -> bytes:
    h = hashlib.sha256()
    h.update(struct.pack(">Q", len(DOMAIN)) + DOMAIN)
    h.update(struct.pack(">Q", len(kind)) + kind)
    h.update(semantics)
    for name, vkind, raw in rows:
        n = name.encode()
        h.update(bytes([TAGS["str"]]) + struct.pack(">Q", len(n)) + n)
        h.update(bytes([TAGS[vkind]]) + value_bytes(vkind, raw))
    return h.digest()


def law_ref(kind: bytes, semantics: bytes, ident: bytes, rows: list[tuple[str, str, bytes]]) -> bytes:
    out = b"LOLR" + struct.pack(">H", 2) + bytes([len(kind)]) + kind + semantics + ident
    out += bytes([len(rows)])
    for name, vkind, raw in rows:
        n = name.encode()
        out += bytes([len(n)]) + n + bytes([TAGS[vkind]]) + value_bytes(vkind, raw)
    return out


def main() -> int:
    if len(sys.argv) != 2 or sys.argv[1] not in ("--check", "--write"):
        print(__doc__)
        return 2
    write = sys.argv[1] == "--write"
    lines = ROWS.read_text(encoding="utf-8").split("\n")
    out, bad, laws = [], [], 0
    kind = semantics = None
    rows: list[tuple[str, str, bytes]] = []
    for line in lines:
        word = line.split(" ", 1)[0]
        if word == "law":
            kind, semantics, rows = line.split(" ", 1)[1].encode(), None, []
        elif word == "semantics":
            semantics = bytes.fromhex(line.split(" ", 1)[1])
        elif word == "row":
            _, name, vkind, text = line.split(" ", 3)
            rows.append((name, vkind, parse_value(vkind, text)))
        elif word in ("expect-id", "expect-ref"):
            if kind is None or semantics is None:
                raise ValueError(f"{word} before law / semantics")
            names = [r[0].encode() for r in rows]
            if names != sorted(set(names)):
                raise ValueError(f"{kind.decode()}: rows are not in strictly ascending name order")
            ident = law_id(kind, semantics, rows)
            got = ident.hex() if word == "expect-id" else law_ref(kind, semantics, ident, rows).hex()
            if word == "expect-id":
                laws += 1
            if write:
                line = f"{word} {got}"
            elif line.split(" ", 1)[1] != got:
                bad.append(f"{kind.decode()}: {word} {line.split(' ', 1)[1]} != computed {got}")
        out.append(line)
    if laws == 0:
        print("law_ref_check: no Law in the file", file=sys.stderr)
        return 1
    if write:
        ROWS.write_text("\n".join(out), encoding="utf-8")
        print(f"law_ref_check: wrote {laws} Laws")
        return 0
    for b in bad:
        print(b, file=sys.stderr)
    print(f"law_ref_check: {laws} Laws compared, {len(bad)} differ")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
