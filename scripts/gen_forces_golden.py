#!/usr/bin/env python3
"""Generate tests/golden/forces.csv: RK4 propagations of an Earth satellite
under the force model of src/forces.rs, computed by an independent reference.

The reference shares no code with the crate: the Sun, the Moon and the GCRS
to CIRS matrix come from ERFA through pyerfa (eraEpv00, eraMoon98,
eraC2i06a); the forces are written here from their definitions (the zonal
potential's gradient by the Legendre recursion of the textbook form, the
third body as the direct difference, the radiation pressure with the conical
shadow of Montenbruck & Gill section 3.4.2, the drag relative to an atmosphere
rotating about the pole); the integrator is the same classical RK4 with the
same steps and the same time arithmetic, so the two propagations differ only
by the rounding of their arithmetic and of ERFA's libm.

Usage (versions pinned in scripts/oracle-requirements.txt):
    .venv/bin/python -I scripts/gen_forces_golden.py           # write
    .venv/bin/python -I scripts/gen_forces_golden.py --check   # compare
"""

import math
import pathlib
import sys

import erfa
import numpy

OUT = pathlib.Path(__file__).resolve().parent.parent / "tests" / "golden" / "forces.csv"
PINS = {"pyerfa": "2.0.1.5", "numpy": "2.5.4"}
AU = 149597870.7
R_SUN = 695700.0
R_EARTH = 6378.137
R_MOON = 1737.4
S_OVER_C = 1361.0 / 299792458.0
OMEGA = 7.292115e-5
MU = 398600.4415
RE = 6378.1363
J = [1.08262668e-3, -2.5327e-6, -1.6196e-6, -2.273e-7, 5.407e-7]
MU_SUN = 132712440041.279419
MU_MOON = 4902.800118


def density(r):
    """exponential atmosphere: 3.725e-12 kg/m^3 at 400 km, scale height 58.515 km"""
    h = math.sqrt(r[0] ** 2 + r[1] ** 2 + r[2] ** 2) - R_EARTH
    return 3.725e-12 * math.exp(-(h - 400.0) / 58.515)


def zonal(r):
    x, y, z = r
    rn = math.sqrt(x * x + y * y + z * z)
    s = z / rn
    rho = RE / rn
    # P_n and P_n' by the textbook recursions, then a = grad U with U = mu/r (1 - sum Jn rho^n Pn)
    P = [1.0, s]
    dP = [0.0, 1.0]
    for n in range(2, len(J) + 2):
        P.append(((2 * n - 1) * s * P[n - 1] - (n - 1) * P[n - 2]) / n)
        dP.append(n * P[n - 1] + s * dP[n - 1])
    radial = 1.0
    axial = 0.0
    for n, jn in enumerate(J, start=2):
        radial -= jn * rho ** n * ((n + 1) * P[n] + s * dP[n])
        axial += jn * rho ** n * dP[n]
    f = -MU / (rn * rn)
    return [f * radial * x / rn, f * radial * y / rn, f * (radial * z / rn + axial)]


def third(r, s, mu):
    d = [r[i] - s[i] for i in range(3)]
    nd = math.sqrt(sum(c * c for c in d))
    ns = math.sqrt(sum(c * c for c in s))
    return [-mu * (d[i] / nd ** 3 + s[i] / ns ** 3) for i in range(3)]


def visible(sat, sun, body, radius):
    """Montenbruck & Gill 3.4.2: apparent radii a (Sun), b (body), separation c"""
    to_sun = [sun[i] - sat[i] for i in range(3)]
    to_body = [body[i] - sat[i] for i in range(3)]
    ds = math.sqrt(sum(c * c for c in to_sun))
    db = math.sqrt(sum(c * c for c in to_body))
    a = math.asin(R_SUN / ds)
    b = math.asin(radius / db)
    c = math.acos(max(-1.0, min(1.0, sum(to_sun[i] * to_body[i] for i in range(3)) / (ds * db))))
    if c >= a + b:
        return 1.0
    if c <= b - a:
        return 0.0
    if c <= a - b:
        return 1.0 - (b * b) / (a * a)
    x = (c * c + a * a - b * b) / (2 * c)
    y = math.sqrt(max(0.0, a * a - x * x))
    area = a * a * math.acos(x / a) + b * b * math.acos((c - x) / b) - c * y
    return 1.0 - area / (math.pi * a * a)


def srp(sat, sun, moon, cr, am):
    nu = visible(sat, sun, [0.0, 0.0, 0.0], R_EARTH) * visible(sat, sun, moon, R_MOON)
    d = [sat[i] - sun[i] for i in range(3)]
    dn = math.sqrt(sum(c * c for c in d))
    k = nu * S_OVER_C * cr * am * (AU / dn) ** 2 / 1000.0 / dn
    return [k * c for c in d]


def drag(rc, vc, cd, am):
    vr = [vc[0] + OMEGA * rc[1], vc[1] - OMEGA * rc[0], vc[2]]
    sp = math.sqrt(sum(c * c for c in vr))
    k = -0.5 * density(rc) * cd * am * sp * 1000.0
    return [k * c for c in vr]


def accel(case, t, r, v):
    jd1, jd2 = case["epoch"]
    t1, t2 = jd1, jd2 + t / 86400.0
    m = numpy.asarray(erfa.c2i06a(t1, t2))
    rc = list(m @ numpy.array(r))
    ac = zonal(rc)
    if case["drag"]:
        vc = list(m @ numpy.array(v))
        ac = [ac[i] + x for i, x in enumerate(drag(rc, vc, *case["drag"]))]
    a = list(m.T @ numpy.array(ac))
    pvh, _ = erfa.epv00(t1, t2)
    sun = [-float(c) * AU for c in pvh[0]]
    moon = [float(c) * AU for c in erfa.moon98(t1, t2)[0]]
    a = [a[i] + x for i, x in enumerate(third(r, sun, MU_SUN))]
    a = [a[i] + x for i, x in enumerate(third(r, moon, MU_MOON))]
    if case["srp"]:
        a = [a[i] + x for i, x in enumerate(srp(r, sun, moon, *case["srp"]))]
    return [float(c) for c in a]


def propagate(case):
    s = list(case["r"]) + list(case["v"])
    dt = case["dt"]
    out = [s]
    for k in range(case["steps"]):
        t = k * dt

        def f(tt, st):
            a = accel(case, tt, st[:3], st[3:])
            return st[3:] + a

        k1 = f(t, s)
        k2 = f(t + dt * 0.5, [s[i] + k1[i] * (dt * 0.5) for i in range(6)])
        k3 = f(t + dt * 0.5, [s[i] + k2[i] * (dt * 0.5) for i in range(6)])
        k4 = f(t + dt, [s[i] + k3[i] * dt for i in range(6)])
        h6 = dt / 6.0
        s = [h6 * ((2.0 * k3[i] + (2.0 * k2[i] + k1[i])) + k4[i]) + s[i] for i in range(6)]
        out.append(s)
    return out


CASES = {
    # 400 km, 51.6 degrees, through the Earth's shadow; every force
    "leo": {"epoch": (2461100.5, 0.25), "r": (6778.137, 0.0, 0.0),
            "v": (0.0, 7.6686 * math.cos(math.radians(51.6)), 7.6686 * math.sin(math.radians(51.6))),
            "dt": 20.0, "steps": 300, "srp": (1.3, 0.02), "drag": (2.2, 0.02)},
    # geostationary, one day; no drag
    "geo": {"epoch": (2461100.5, 0.0), "r": (42164.0, 0.0, 0.0), "v": (0.0, 3.0746, 0.0),
            "dt": 600.0, "steps": 144, "srp": (1.5, 0.01), "drag": None},
}


def main():
    got = {"pyerfa": erfa.__version__, "numpy": numpy.__version__}
    if got != PINS:
        sys.exit(f"tool versions {got} differ from the pins {PINS}")
    lines = ["# generated by scripts/gen_forces_golden.py (pyerfa + forces written from their definitions)",
             "case,step,x,y,z,vx,vy,vz"]
    for name, case in CASES.items():
        for k, s in enumerate(propagate(case)):
            if k % 10 == 0 or k == case["steps"]:
                lines.append(f"{name},{k}," + ",".join(repr(float(c)) for c in s))
    text = "\n".join(lines) + "\n"
    if "--check" in sys.argv[1:]:
        if OUT.read_text(encoding="utf-8") != text:
            sys.exit(f"{OUT} differs from the regenerated table")
        print(f"{OUT.name}: {len(lines) - 2} rows match")
    else:
        OUT.write_text(text, encoding="utf-8")
        print(f"wrote {OUT.name}: {len(lines) - 2} rows")


if __name__ == "__main__":
    main()
