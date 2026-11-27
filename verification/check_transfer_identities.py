"""Exact identities behind Sections 6 and 7, checked on small examples.

1. Parseval step of Proposition 6.1: for any signs and any finite set of z-smooth n,
       int_0^inf ( sum_{u<n<=(1+delta)u} f(n) )^2 du/u^2
     = e^{1/T}/(2 pi T^2) int W_T(t) |R_z(t)|^2 dt ,   T = 2/log(1+delta).
   The left side is computed exactly (piecewise constant integrand); the right side
   both via the closed form of hat W_T and by direct numerical integration.
2. Normalization: Theta_z M_z x (e^{2/T}-1) = y prod_{p<=z}(1-p^-2).
3. Smooth/rough decomposition: sum_{x<n<=X} f(n) = (smooth part) + sum_m a_m f(m)
   and S_{>z} = sum_p D_p, for random sign realizations.
4. Tower identity (eq:tower): E(Q_x | F_z) = V_f, computed exactly by expanding each
   B_p as a polynomial in the rough signs.
"""
import itertools
import numpy as np
from scipy.integrate import quad
from common import check, summary, primes_upto, W_T, What_closed

print("check_transfer_identities.py")
rng = np.random.default_rng(11)


def factor(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            e = 0
            while n % d == 0:
                n //= d; e += 1
            out.append((d, e))
        d += 1
    if n > 1:
        out.append((n, 1))
    return out


def sqfree_primes(n):
    f = factor(n)
    if any(e > 1 for _, e in f):
        return None
    return [p for p, _ in f]


# ---------- 1. Parseval identity ----------
z = 13
pz = [int(p) for p in primes_upto(z)]
smooth_sf = sorted(int(np.prod(c)) for r in range(len(pz) + 1) for c in itertools.combinations(pz, r))

for trial, delta in enumerate([0.3, 0.05, 0.01]):
    eps = {p: rng.choice([-1, 1]) for p in pz}
    fval = {n: int(np.prod([eps[p] for p in pz if n % p == 0])) for n in smooth_sf}
    T = 2 / np.log1p(delta)
    # exact LHS: integrand is constant between consecutive breakpoints n and n/(1+delta)
    bps = sorted(set([n for n in smooth_sf] + [n / (1 + delta) for n in smooth_sf]))
    lhs = 0.0
    for a, b in zip(bps[:-1], bps[1:]):
        u = 0.5 * (a + b)
        S = sum(fval[n] for n in smooth_sf if u < n <= (1 + delta) * u)
        lhs += S * S * (1 / a - 1 / b)
    # RHS via closed-form Fourier transform: int W_T |R|^2 = sum_{m,n} f f (mn)^{-1/2} hatW(log(n/m))
    arr = np.array(smooth_sf, float); fv = np.array([fval[n] for n in smooth_sf], float)
    M = np.outer(fv, fv) / np.sqrt(np.outer(arr, arr)) * What_closed(np.log(arr[None, :] / arr[:, None]), T)
    rhs_closed = np.exp(1 / T) / (2 * np.pi * T * T) * M.sum()
    check(f"Parseval identity, closed-form side (delta={delta})", abs(lhs - rhs_closed) < 1e-9 * max(1, abs(lhs)),
          f"LHS={lhs:.10f}, RHS={rhs_closed:.10f}")
    if trial == 0:  # independent numerical integration in t (slow, do once)
        def integrand(t):
            R = np.sum(fv * arr ** (-0.5 - 1j * t))
            return W_T(t, T) * abs(R) ** 2
        val = 2 * sum(quad(integrand, k, k + 5, limit=200)[0] for k in range(0, 4000, 5))
        # tail beyond 4000: W_T <= 4 T^2/t^2 and |R|^2 <= (sum n^{-1/2})^2
        tail = 2 * 4 * T * T / 4000 * np.sum(arr ** -0.5) ** 2
        rhs_num = np.exp(1 / T) / (2 * np.pi * T * T) * val
        check("Parseval identity, direct t-integration", abs(lhs - rhs_num) < 2e-3 * lhs + np.exp(1 / T) / (2 * np.pi * T * T) * tail,
              f"num RHS={rhs_num:.6f} (tail bound {np.exp(1/T)/(2*np.pi*T*T)*tail:.1e})")

# ---------- 2. normalization ----------
for x, alpha in [(1e8, 0.2), (1e12, 0.35), (1e20, 0.1)]:
    y = x / np.log(x) ** alpha
    T = 2 / np.log1p(y / x)
    zz = 1000
    pr = primes_upto(zz).astype(float)
    Theta = np.prod(1 - 1 / pr); Mz = np.prod(1 + 1 / pr)
    J = 4 * np.pi * T * T * np.sinh(1 / T)
    lhs = Theta * x * np.exp(1 / T) / (2 * np.pi * T * T) * Mz * J
    rhs = y * np.prod(1 - pr ** -2)
    check(f"normalization identity (x={x:.0e})", abs(lhs / rhs - 1) < 1e-10)
check("prod(1-p^-2) -> 6/pi^2", abs(np.prod(1 - primes_upto(10**6).astype(float) ** -2) - 6 / np.pi ** 2) < 1e-6)

# ---------- 3 & 4. decompositions and the tower identity ----------
x, y, z = 3000, 400, 7
X = x + y
primes = [int(p) for p in primes_upto(X)]
small = [p for p in primes if p <= z]
rough_primes = [p for p in primes if p > z]
sf = {n: sqfree_primes(n) for n in range(1, X + 1)}

ok3, ok4 = True, True
for trial in range(5):
    eps = {p: int(rng.choice([-1, 1])) for p in primes}
    f = lambda n: 0 if sf[n] is None else int(np.prod([eps[p] for p in sf[n]])) if sf[n] else 1
    total = sum(f(n) for n in range(x + 1, X + 1))
    smooth = sum(f(n) for n in range(x + 1, X + 1) if all(p <= z for p in (sf[n] or [1])) and sf[n] is not None)
    # a_m for rough m
    def a(m):
        return sum(f(n) for n in range(1, X // m + 1) if x < m * n <= X and sf[n] is not None and all(p <= z for p in sf[n]))
    rough_ms = [m for m in range(2, X + 1) if all(p > z for p, _ in factor(m))]
    S_gt = sum(a(m) * f(m) for m in rough_ms)
    ok3 &= (total == smooth + S_gt)
    # martingale differences D_p
    Dsum = 0
    for p in rough_primes:
        Bp = sum(f(m) for m in range(1, X // p + 1) if x < p * m <= X and (m == 1 or max(q for q, _ in factor(m)) < p))
        Dsum += eps[p] * Bp
    ok3 &= (Dsum == S_gt)
    # Tower identity: expand B_p as polynomial in rough signs given smooth signs.
    Vf = sum((a(m) ** 2) for m in rough_ms if sqfree_primes(m) is not None)
    EQ = 0
    for p in rough_primes:
        poly = {}
        for m in range(1, X // p + 1):
            if not (x < p * m <= X) or sf[m] is None:
                continue
            if m > 1 and max(sf[m]) >= p:
                continue
            sm = [q for q in sf[m] if q <= z]; ro = tuple(q for q in sf[m] if q > z)
            coef = int(np.prod([eps[q] for q in sm])) if sm else 1
            poly[ro] = poly.get(ro, 0) + coef
        EQ += sum(c * c for c in poly.values())  # orthonormality of rough monomials
    ok4 &= (EQ == Vf)
check("sum = smooth part + S_{>z}, and S_{>z} = sum_p D_p (5 realizations)", ok3)
check("E(Q_x | F_z) = V_f exactly (eq:tower), 5 realizations", ok4)

raise SystemExit(summary("check_transfer_identities.py"))
