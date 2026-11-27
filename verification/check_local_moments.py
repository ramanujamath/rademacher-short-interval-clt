"""Exact single-prime identities and Taylor expansions used in Sections 2, 4, 5, 8, 9, 10.

Each identity is checked by averaging exactly over the two values of the sign.
"""
import numpy as np
from common import check, summary, primes_upto

print("check_local_moments.py")
rng = np.random.default_rng(7)
P = primes_upto(10**5)


def fac(eps, p, t):
    return 1 + eps * p ** (-0.5 - 1j * t)


# (eq:4mom) per-prime factor:  E|f(t)|^2|f(u)|^2 / (1+1/p)^2 = 1 + 4p/(p+1)^2 cos cos
worst = 0.0
for _ in range(2000):
    p = int(rng.choice(P[:2000])); t, u = rng.uniform(-300, 300, 2)
    lhs = np.mean([abs(fac(e, p, t)) ** 2 * abs(fac(e, p, u)) ** 2 for e in (1, -1)]) / (1 + 1 / p) ** 2
    rhs = 1 + 4 * p / (p + 1) ** 2 * np.cos(t * np.log(p)) * np.cos(u * np.log(p))
    worst = max(worst, abs(lhs - rhs))
check("exact fourth-moment factor (eq:4mom)", worst < 1e-12, f"max err {worst:.1e}")

# E|R(t)|^2 = M : per-prime mean of |1+eps p^{-1/2-it}|^2 is 1+1/p
worst = max(abs(np.mean([abs(fac(e, p, t)) ** 2 for e in (1, -1)]) - (1 + 1 / p))
            for p, t in zip(rng.choice(P[:500], 500), rng.uniform(-100, 100, 500)))
check("E|1+eps p^{-1/2-it}|^2 = 1+1/p", worst < 1e-12)

# (eq:4momlog): log of factor = 2[cos(t-u)+cos(t+u)]/p + O(p^-2); check |err| p^2 bounded
C = 0.0
for p in P[P > 50][:3000]:
    t, u = rng.uniform(-300, 300, 2)
    lp = np.log(p)
    exact = np.log(1 + 4 * p / (p + 1) ** 2 * np.cos(t * lp) * np.cos(u * lp))
    approx = 2 * (np.cos((t - u) * lp) + np.cos((t + u) * lp)) / p
    C = max(C, abs(exact - approx) * p * p)
check("log-factor expansion error is O(p^-2) (eq:4momlog)", C < 20, f"sup p^2|err| = {C:.2f}")

# (eq:low4) bound: (1+1/p)^2 + (4/p)|cos cos| <= 1 + 6/p + 1/p^2
ok = all((1 + 1 / p) ** 2 + 4 / p <= 1 + 6 / p + 1 / p ** 2 + 1e-15 for p in P[:1000])
check("low-prime fourth-moment bound (eq:low4)", ok)

# Tilted measure (Lemma 4.1): E_Q eps = 2p^{-1/2}cos/(1+1/p); E_Q Lambda_p expansion (eq:tiltmean)
C1, C2, C3 = 0.0, 0.0, 0.0
for p in P[P > 20][:3000]:
    t = rng.uniform(1, 1e4)
    w = {e: abs(fac(e, p, t)) ** 2 / (2 * (1 + 1 / p)) for e in (1, -1)}  # Q-probabilities
    C1 = max(C1, abs(sum(w.values()) - 1))
    EQeps = w[1] - w[-1]
    C2 = max(C2, abs(EQeps - 2 * p ** -0.5 * np.cos(t * np.log(p)) / (1 + 1 / p)))
    Lam = {e: np.log(abs(fac(e, p, t))) for e in (1, -1)}
    EQL = w[1] * Lam[1] + w[-1] * Lam[-1]
    pred = 1 / p + np.cos(2 * t * np.log(p)) / (2 * p)
    C3 = max(C3, abs(EQL - pred) * p * p)
check("tilted weights form a probability measure", C1 < 1e-12)
check("E_Q eps_p = 2p^{-1/2}cos(t log p)/(1+1/p)", C2 < 1e-12)
check("E_Q Lambda_p = 1/p + cos(2t log p)/(2p) + O(p^-2) (eq:tiltmean)", C3 < 20, f"sup p^2|err| = {C3:.2f}")

# Lemma 9.1 expansion of log E|f(t)||f(u)|
C = 0.0
for p in P[P > 20][:3000]:
    t, u = rng.uniform(1, 1e4, 2)
    lp = np.log(p)
    exact = np.log(np.mean([abs(fac(e, p, t)) * abs(fac(e, p, u)) for e in (1, -1)]))
    pred = (1 / (2 * p) - (np.cos(2 * t * lp) + np.cos(2 * u * lp)) / (4 * p)
            + (np.cos((t - u) * lp) + np.cos((t + u) * lp)) / (2 * p))
    C = max(C, abs(exact - pred) * p * p)
check("two-point expansion in Lemma 9.1 is accurate to O(p^-2)", C < 20, f"sup p^2|err| = {C:.2f}")

# Martingale CLT (Lemma 7.2): e^{x^2/2}cos x - 1 = -x^4/12 + O(x^6); bound C_t b^4 for |b|<=sqrt2
x = np.linspace(-3, 3, 60001); x = x[np.abs(x) > 0.05]  # avoid float cancellation at x~0
r = (np.exp(x ** 2 / 2) * np.cos(x) - 1) / x ** 4
check("e^{x^2/2}cos x - 1 = -x^4/12 + O(x^6)", abs(r[np.argmin(abs(x))] + 1 / 12) < 1e-3)
check("|e^{x^2/2}cos x - 1| <= C x^4 on |x| <= 3", np.max(abs(r)) < 2, f"C = {np.max(abs(r)):.3f}")

# Relaxation in Lemma 8.2: (1+6/p)(1-1/p)^6 in (0,1] and = 1 - O(p^-2)
vals = np.array([(1 + 6 / p) * (1 - 1 / p) ** 6 for p in P])
check("(1+6/p)(1-1/p)^6 lies in (0,1]", np.all((vals > 0) & (vals <= 1)))
check("1-(1+6/p)(1-1/p)^6 = O(p^-2)", np.max((1 - vals) * P.astype(float) ** 2) < 30,
      f"sup p^2(1-.) = {np.max((1 - vals) * P.astype(float) ** 2):.2f}")

# Walsh expansion (Lemma 3.1): |1+eps p^{-1/2-it}|^2/(1+1/p) = 1 + 2sqrt(p)/(p+1) eps cos
worst = max(abs(abs(fac(e, p, t)) ** 2 / (1 + 1 / p) - (1 + 2 * np.sqrt(p) / (p + 1) * e * np.cos(t * np.log(p))))
            for e in (1, -1) for p, t in zip(rng.choice(P[:300], 300), rng.uniform(-50, 50, 300)))
check("Walsh expansion of normalized local factor", worst < 1e-12)

raise SystemExit(summary("check_local_moments.py"))
