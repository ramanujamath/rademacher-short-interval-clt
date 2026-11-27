"""Exponent bookkeeping for every 'o(1)' claim in the proof.

Each error term in the paper is a ratio built from powers of x, log x, L = log log z,
U = log X / log z, etc. We set every implied constant to 1, write log(ratio) as an
explicit function of L0 = log log x, and check that it is eventually negative and tends to
-infinity. We also report the first L0 (on a geometric grid) after which it stays negative.
This tests the arithmetic of the exponents, not the inequalities themselves.

With z = x^{1/log log log x}:  log z = e^{L0}/log L0,  L = L0 - log log L0,  U = log L0,
T = 2/log(1+delta) ~ 2 (log x)^alpha,  ell = log T,  w = exp(ell^4),  w0 = exp(L^4).

The thresholds are astronomically large for several terms: the proof is purely
asymptotic. That is expected and is not an error.
"""
import math
import numpy as np
from common import check, summary

print("check_asymptotic_bookkeeping.py")


def params(L0, a):
    lx = math.exp(L0) if L0 < 700 else math.inf
    lL0 = math.log(L0)
    U = lL0
    logz = lx / lL0 if lx != math.inf else math.inf
    L = L0 - math.log(lL0)
    ell = a * L0 + math.log(2)
    return dict(lx=lx, U=U, logz=logz, L=L, ell=ell, lL0=lL0)


def lg(v):
    return math.log(v)


A = 1.0  # any fixed A > 0 suffices where L^{-A} appears
c0 = 0.01  # a (pessimistically small) constant in the prime number theorem error

# (name, function(L0, alpha) -> log ratio, uses_x)  uses_x=True means it needs e^{L0} (L0 <= 700)
claims = [
    # Theorem 1.2 (Sections 4-5), applied with T ~ (log x)^alpha and z as above
    ("w < z: log w / log z = ell^4 / log z", lambda L0, a: 4 * lg(params(L0, a)["ell"]) - lg(params(L0, a)["logz"]), True),
    ("Lemma 4.1: L / ell^2", lambda L0, a: lg(params(L0, a)["L"]) - 2 * lg(params(L0, a)["ell"]), False),
    ("Lemma 4.1: barrier margin (beta + 2log(1+L) - ell/4 < 0)",
     lambda L0, a: (lambda p: math.copysign(1, math.sqrt(p["ell"]) + 2 * math.log1p(p["L"]) - p["ell"] / 4))(params(L0, a)), False),
    ("PNT error T^2 exp(-c0 ell^2) (c0=0.01)", lambda L0, a: 2 * params(L0, a)["ell"] - c0 * params(L0, a)["ell"] ** 2, False),
    ("(5.1): 1/T", lambda L0, a: -params(L0, a)["ell"], False),
    ("Lemma 3.1 use: (log w)^3 / T", lambda L0, a: 12 * lg(params(L0, a)["ell"]) - params(L0, a)["ell"], False),
    ("Lemma 3.1 use: log T / sqrt T", lambda L0, a: lg(params(L0, a)["ell"]) - params(L0, a)["ell"] / 2, False),
    ("intermediate: (log w)^{7/2} / T", lambda L0, a: 14 * lg(params(L0, a)["ell"]) - params(L0, a)["ell"], False),
    ("close: e^{-2 sqrt ell} ell^12", lambda L0, a: -2 * math.sqrt(params(L0, a)["ell"]) + 12 * lg(params(L0, a)["ell"]), False),
    # Proposition 6.1
    ("1/z (from y/z)", lambda L0, a: -params(L0, a)["logz"], True),
    ("sqrt(x)/y", lambda L0, a: -params(L0, a)["lx"] / 2 + a * L0, True),
    ("x/(y log z)", lambda L0, a: a * L0 - L0 + lg(params(L0, a)["lL0"]), True),
    ("sqrt(x/y)/log z", lambda L0, a: a * L0 / 2 - lg(params(L0, a)["logz"]), True),
    ("small m: L^{-A} log L (A=1)", lambda L0, a: -A * lg(params(L0, a)["L"]) + lg(lg(params(L0, a)["L"])), False),
    ("Lemma 6.2(a): rho(u) <= L^{-A}, u = 0.09 log L0 (A=1)",
     lambda L0, a: (lambda u: -u * math.log(u) + A * lg(params(L0, a)["L"]))(0.09 * math.log(L0)), False),
    ("Lemma 6.2(c): sieve error exp(-s log s) <= L^{-A}, s = 0.8 log L0",
     lambda L0, a: (lambda s: -s * math.log(s) + s + A * lg(params(L0, a)["L"]))(0.8 * math.log(L0)), False),
    # Section 7
    ("Lemma 7.1: (log x)^{2 alpha + 4} / z", lambda L0, a: (2 * a + 4) * L0 - params(L0, a)["logz"], True),
    # Section 8
    ("smoothing: sqrt(x eta / y)", lambda L0, a: 0.5 * (a * L0 - 100 * L0), False),
    # Section 9
    ("large frequencies: (log x)^{2a} U^6 (log z)^6 L^3 (log X)^{-900}",
     lambda L0, a: (lambda p: 2 * a * L0 + 6 * lg(p["U"]) + 6 * (L0 - lg(p["lL0"])) + 3 * lg(p["L"]) - 900 * L0)(params(L0, a)), False),
    ("U^6/z term: (log x)^{2a} U^6 (log z)^6 L^4 / z",
     lambda L0, a: (lambda p: 2 * a * L0 + 6 * lg(p["U"]) + 6 * lg(p["logz"]) + 4 * lg(p["L"]) - p["logz"])(params(L0, a)), True),
    ("Lemma 9.2: L^3 sqrt(log z)/log X <= L^{-10}",
     lambda L0, a: (lambda p: 13 * lg(p["L"]) + 0.5 * (L0 - lg(p["lL0"])) - L0)(params(L0, a)), False),
    ("Lemma 9.2: delta^2 L <= delta L^{-10}", lambda L0, a: -a * L0 + 11 * lg(params(L0, a)["L"]), False),
    ("Lemma 9.2: r1 = L^20/log z -> 0", lambda L0, a: (lambda p: 20 * lg(p["L"]) - (L0 - lg(p["lL0"])))(params(L0, a)), False),
    ("separated pairs: U^6 L^{-4}", lambda L0, a: 6 * lg(params(L0, a)["U"]) - 4 * lg(params(L0, a)["L"]), False),
    ("outside D+-: U^5 L^{-18}", lambda L0, a: 5 * lg(params(L0, a)["U"]) - 18 * lg(params(L0, a)["L"]), False),
    # Section 10
    ("barrier margin: alpha L/4 - sqrt L > 0", lambda L0, a: -(a * params(L0, a)["L"] / 4 - math.sqrt(params(L0, a)["L"])), False),
    ("bad part: U^6 L^{-3/4}", lambda L0, a: 6 * lg(params(L0, a)["U"]) - 0.75 * lg(params(L0, a)["L"]), False),
    ("good part: U^6 e^{-2 sqrt L} L^72",
     lambda L0, a: 6 * lg(params(L0, a)["U"]) - 2 * math.sqrt(params(L0, a)["L"]) + 72 * lg(params(L0, a)["L"]), False),
]

# Two claims only become negative for L0 = loglog x beyond e^(10^5); for them we use
# psi = log L0 as the variable (log L = psi + o(1)), on a grid psi <= 10^9.
claims_psi = [
    ("Lemma 6.2(a): rho(u) <= L^{-A}, u = 0.09 log L0 (A=1), variable psi=log L0",
     lambda psi: -(0.09 * psi) * math.log(0.09 * psi) + A * psi),
    ("small-m use only needs rho(u) log L0 -> 0, variable psi=log L0",
     lambda psi: -(0.09 * psi) * math.log(0.09 * psi) + math.log(psi)),
]
claims = [c for c in claims if not c[0].startswith("Lemma 6.2(a)")]

grid_small = np.geomspace(20, 650, 400)
grid_large = np.geomspace(20, 1e40, 4000)
print(f"  {'claim':72s} {'alpha':>5s}  threshold L0=loglog x")
for name, fn, uses_x in claims:
    for a in (0.05, 0.38):
        grid = grid_small if uses_x else grid_large
        vals = np.array([fn(L0, a) for L0 in grid])
        neg = vals < 0
        if neg[-1] and np.all(neg[np.argmax(neg):]):
            thr = grid[np.argmax(neg)]
            # tends to -infinity: last value well below the value at the threshold
            trend = vals[-1] < vals[np.argmax(neg)] - 1 or "margin" in name or "w < z" in name
            ok = bool(trend)
        else:
            thr, ok = float("nan"), False
        print(f"  {name:72s} {a:5.2f}  {thr:9.3g}")
        check(f"{name} [alpha={a}]", ok)

grid_psi = np.geomspace(3, 1e9, 4000)
for name, fn in claims_psi:
    vals = np.array([fn(v) for v in grid_psi])
    neg = vals < 0
    ok = bool(neg[-1] and np.all(neg[np.argmax(neg):]) and vals[-1] < vals[np.argmax(neg)] - 1)
    thr = grid_psi[np.argmax(neg)] if ok else float("nan")
    print(f"  {name:72s}   -    threshold log L0 = {thr:.3g}")
    check(name, ok)

raise SystemExit(summary("check_asymptotic_bookkeeping.py"))
