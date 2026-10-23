#!/usr/bin/env python3
"""Generate tests/golden/eclipse.csv: the fraction of the Sun's disk seen by a
satellite past an occulting body, for positions given in kilometres.

Independent of the crate's closed form: the apparent radii and the
separation come from mpmath at 40 digits, and the overlap of the two disks
is integrated numerically (mpmath.quad of the overlapping chord heights,
split at the crossing abscissa), not taken from the circle-intersection
formula.

Usage (version pinned in scripts/oracle-requirements.txt):
    python3 -m venv .venv && .venv/bin/pip install -r scripts/oracle-requirements.txt
    .venv/bin/python -I scripts/gen_eclipse_golden.py           # write
    .venv/bin/python -I scripts/gen_eclipse_golden.py --check   # compare
"""

import pathlib
import sys

import mpmath

ROOT = pathlib.Path(__file__).resolve().parent.parent
OUT = ROOT / "tests" / "golden" / "eclipse.csv"
GRAZING_OUT = ROOT / "tests" / "golden" / "eclipse_grazing.csv"
PIN = "1.4.1"
mpmath.mp.dps = 40

R_SUN = mpmath.mpf(695700)
AU = mpmath.mpf("149597870.7")
BODIES = {"earth": mpmath.mpf("6378.137"), "moon": mpmath.mpf("1737.4")}


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def norm(v):
    return mpmath.sqrt(sum(x * x for x in v))


def overlap(a, b, c):
    """Area of the intersection of disks of radii a, b at centre distance c."""
    if c >= a + b:
        return mpmath.mpf(0)
    if c <= abs(a - b):
        return mpmath.pi * min(a, b) ** 2
    lo, hi = max(-a, c - b), min(a, c + b)
    x0 = (c * c + a * a - b * b) / (2 * c)

    def h(x):
        ha = mpmath.sqrt(max(a * a - x * x, 0))
        hb = mpmath.sqrt(max(b * b - (x - c) ** 2, 0))
        return 2 * min(ha, hb)

    return mpmath.quad(h, [lo, x0, hi])


def fraction(sat, sun, body, r_body):
    to_sun, to_body = sub(sun, sat), sub(body, sat)
    a = mpmath.asin(R_SUN / norm(to_sun))
    b = mpmath.asin(r_body / norm(to_body))
    cross = [to_sun[1] * to_body[2] - to_sun[2] * to_body[1],
             to_sun[2] * to_body[0] - to_sun[0] * to_body[2],
             to_sun[0] * to_body[1] - to_sun[1] * to_body[0]]
    dot = sum(x * y for x, y in zip(to_sun, to_body))
    c = mpmath.atan2(norm(cross), dot)
    if c >= a + b:
        kind = "none"
    elif c <= b - a:
        kind = "total"
    elif c <= a - b:
        kind = "annular"
    else:
        kind = "partial"
    nu = 1 - overlap(a, b, c) / (mpmath.pi * a * a)
    return kind, nu, a, b, c


def cases():
    """Satellites near the shadow of the Earth (Sun on +x at 1 au) and of the
    Moon, swept across umbra, penumbra and annular regions."""
    out = []
    sun = [AU, mpmath.mpf(0), mpmath.mpf(0)]
    re = BODIES["earth"]
    # LEO to beyond GEO, behind the Earth, offset across the shadow edge
    for dist in [mpmath.mpf(v) for v in ("6778", "7178", "26560", "42164", "400000", "1500000")]:
        for k in range(-24, 25):
            off = re * (mpmath.mpf(k) / 20)
            out.append(("earth", [-dist, off, off / 3], sun, [0, 0, 0], re))
    # the Moon at 384400 km towards the Sun: satellites behind it at
    # distances that see the Moon larger and smaller than the Sun
    rm = BODIES["moon"]
    moon = [mpmath.mpf(384400), mpmath.mpf(0), mpmath.mpf(0)]
    for back in [mpmath.mpf(v) for v in ("5000", "50000", "374000", "380000", "420000")]:
        for k in range(-16, 17):
            off = rm * (mpmath.mpf(k) / 6)
            out.append(("moon", [moon[0] - back, off, -off / 2], sun, moon, rm))
    # just inside and just outside each edge (c = a + b, the penumbra; c = b - a,
    # the umbra) found at 40 digits: the offsets 1e-12 apart in y change the
    # separation by about 1e-14 rad, far above the f64 rounding of c (1e-18 rad)
    for name, x, body, rb in [("earth", mpmath.mpf(-7000), [0, 0, 0], re),
                              ("earth", mpmath.mpf(-42164), [0, 0, 0], re),
                              ("earth", mpmath.mpf(-400000), [0, 0, 0], re),
                              ("moon", mpmath.mpf(380000), moon, rm)]:
        for edge in ("penumbra", "umbra"):
            y = edge_offset(x, sun, body, rb, edge)
            if y is None:
                continue
            for f in (1 - mpmath.mpf("1e-12"), 1 + mpmath.mpf("1e-12")):
                # the row is read back as f64: classify the f64 input itself
                out.append((name, [x, mpmath.mpf(float(y * f)), mpmath.mpf(0)], sun, body, rb))
    return out


def edge_offset(x, sun, body, rb, edge):
    """The y > 0 at which the satellite (x, y, 0) sits on the edge, by bisection
    at 40 digits; None when the edge does not exist at that x."""
    def gap(y):
        _, _, a, b, c = fraction([x, y, mpmath.mpf(0)], sun, body, rb)
        return c - {"penumbra": a + b, "umbra": b - a, "annular": a - b}[edge]
    lo, hi = mpmath.mpf(0), 4 * rb
    if gap(lo) >= 0 or gap(hi) <= 0:
        return None
    for _ in range(160):
        mid = (lo + hi) / 2
        if gap(mid) < 0:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def textbook(sat, sun, body, rb):
    """The visible fraction from the circle-intersection formula at 50 digits
    (a different formula from the crate's; at this precision its cancellation
    near tangency costs nothing), and the two segment parameters u = 2*theta."""
    with mpmath.workdps(50):
        kind, _, a, b, c = fraction(sat, sun, body, rb)
        if kind != "partial":
            return kind, (1 if kind == "none" else (1 - b * b / (a * a) if kind == "annular" else 0)), None, None
        x = (c * c + a * a - b * b) / (2 * c)
        y = mpmath.sqrt(a * a - x * x)
        ta, tb = mpmath.acos(x / a), mpmath.acos((c - x) / b)
        lens = a * a * ta + b * b * tb - c * y
        return kind, 1 - lens / (mpmath.pi * a * a), 2 * ta, 2 * tb


def grazing():
    """Satellites within 1e-13 .. 1e-1 (relative offset in y) of every edge
    (penumbra c = a + b, umbra c = b - a, annular c = a - b), and around the
    switch of the crate's segment series (u = 2*theta = 0.25) for each disk."""
    sun = [AU, mpmath.mpf(0), mpmath.mpf(0)]
    re, rm = BODIES["earth"], BODIES["moon"]
    moon = [mpmath.mpf(384400), mpmath.mpf(0), mpmath.mpf(0)]
    configs = [("earth", mpmath.mpf(-d), [0, 0, 0], re)
               for d in (6778, 7000, 8000, 12000, 26560, 42164, 100000, 400000, 1000000, 1500000)]
    configs += [("moon", moon[0] - mpmath.mpf(back), moon, rm) for back in (5000, 50000, 300000, 380000, 420000)]
    deltas = [mpmath.mpf(10) ** (-13 + 12 * k / 39) for k in range(40)]
    rows = []
    for name, x, body, rb in configs:
        for edge in ("penumbra", "umbra", "annular"):
            y = edge_offset(x, sun, body, rb, edge)
            if y is None:
                continue
            for d in deltas:
                for f in (1 - d, 1 + d):
                    rows.append((name, [x, mpmath.mpf(float(y * f)), mpmath.mpf(0)], sun, body, rb))
        for which in (2, 3):
            # the y at which u of the Sun's (2) or the body's (3) segment is 0.25
            lo, hi = edge_offset(x, sun, body, rb, "umbra") or mpmath.mpf(0), edge_offset(x, sun, body, rb, "penumbra")
            if hi is None:
                continue
            def u_at(yy):
                r = textbook([x, yy, mpmath.mpf(0)], sun, body, rb)
                return r[which] if r[which] is not None else mpmath.mpf(0)
            a_, b_ = lo, hi
            if (u_at((a_ + b_) / 2 * 0 + a_ * (1 + mpmath.mpf("1e-9"))) - mpmath.mpf("0.25")) * (u_at(b_ * (1 - mpmath.mpf("1e-9"))) - mpmath.mpf("0.25")) > 0:
                continue
            for _ in range(160):
                m = (a_ + b_) / 2
                if (u_at(m) - mpmath.mpf("0.25")) * (u_at(a_ * (1 + mpmath.mpf("1e-9"))) - mpmath.mpf("0.25")) > 0:
                    a_ = m
                else:
                    b_ = m
            ys = (a_ + b_) / 2
            for f in (1 - mpmath.mpf("1e-6"), 1 + mpmath.mpf("1e-6")):
                rows.append((name, [x, mpmath.mpf(float(ys * f)), mpmath.mpf(0)], sun, body, rb))
    return rows


def main():
    if mpmath.__version__ != PIN:
        sys.exit(f"mpmath {mpmath.__version__} != pinned {PIN}")
    lines = [
        "# generated by scripts/gen_eclipse_golden.py (mpmath 1.4.1, 40 digits, quadrature)",
        "body,sat_x,sat_y,sat_z,sun_x,sun_y,sun_z,body_x,body_y,body_z,body_radius,"
        "kind,fraction,sun_angular_radius,body_angular_radius,separation",
    ]
    for name, sat, sun, body, rb in cases():
        kind, nu, a, b, c = fraction(sat, sun, body, rb)
        nums = [*sat, *sun, *body, rb]
        lines.append(
            name + "," + ",".join(repr(float(v)) for v in nums) + f",{kind},"
            + ",".join(mpmath.nstr(v, 25) for v in (nu, a, b, c))
        )
    text = "\n".join(lines) + "\n"
    graze = [
        "# generated by scripts/gen_eclipse_golden.py: grazing configurations, the circle-intersection formula at 50 digits",
        "body,sat_x,sat_y,sat_z,sun_x,sun_y,sun_z,body_x,body_y,body_z,body_radius,kind,fraction,u_sun,u_body",
    ]
    for name, sat, sun, body, rb in grazing():
        kind, nu, ua, ub = textbook(sat, sun, body, rb)
        nums = [*sat, *sun, *body, rb]
        us = ",".join("" if u is None else mpmath.nstr(u, 20) for u in (ua, ub))
        graze.append(name + "," + ",".join(repr(float(v)) for v in nums) + f",{kind},{mpmath.nstr(nu, 25)},{us}")
    graze_text = "\n".join(graze) + "\n"
    if "--check" in sys.argv:
        if GRAZING_OUT.read_text(encoding="utf-8") != graze_text:
            sys.exit(f"{GRAZING_OUT} differs from the regenerated table")
        if OUT.read_text(encoding="utf-8") != text:
            sys.exit(f"{OUT} differs from the regenerated table")
        print(f"{OUT}: {len(lines) - 2} rows match")
    else:
        OUT.parent.mkdir(exist_ok=True)
        OUT.write_text(text, encoding="utf-8")
        GRAZING_OUT.write_text(graze_text, encoding="utf-8")
        print(f"wrote {GRAZING_OUT}: {len(graze) - 2} rows")
        print(f"wrote {OUT}: {len(lines) - 2} rows")


if __name__ == "__main__":
    main()
