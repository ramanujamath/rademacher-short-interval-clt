"""Numerical sanity checks of the prime-sum inputs (Lemma 2.1 and (eq:Ebound)).

1. |sum_{a<p<=b} p^{ih}/p| * |h| log a stays bounded (first term of (eq:pnt)) once
   |h| log a is large, for a moderate range of h.
2. M_{a,b} log a / log b stays between two constants (eq:mertens).
3. |E_Y(v) - 1| / (U Gamma(v)) stays bounded, where E_Y(v) = prod_{z<p<=Y}(1-p^{-1-iv})^{-1},
   U = log X/log z and Gamma(v) = min(1, 1/(|v| log X))  (eq:Ebound).
These are finite-range sanity checks of standard estimates, not proofs.
"""
import numpy as np
from common import check, summary, primes_upto

print("check_prime_sums.py")
P = primes_upto(2 * 10**7).astype(float)
logP = np.log(P)

# 1. oscillatory prime sums
worst = 0.0
for a in [1e3, 1e4, 1e5]:
    sel = (P > a) & (P <= 2e7)
    for h in np.geomspace(0.05, 50, 60):
        if h * np.log(a) < 5:
            continue
        S = np.sum(np.exp(1j * h * logP[sel]) / P[sel])
        worst = max(worst, abs(S) * h * np.log(a))
check("|sum p^{ih}/p| <= C/(|h| log a) on the tested range", worst < 5, f"empirical C = {worst:.2f}")

# 2. Mertens
r = []
for a in [10, 100, 1e3, 1e4]:
    for b in [1e5, 1e6, 2e7]:
        sel = (P > a) & (P <= b)
        r.append(np.exp(np.sum(np.log1p(1 / P[sel]))) * np.log(a) / np.log(b))
check("M_{a,b} asymp log b/log a", 0.5 < min(r) and max(r) < 2, f"range [{min(r):.3f}, {max(r):.3f}]")

# 3. Euler-product bound
# At computable z the prime-number-theorem error (1+|v|)exp(-c0 sqrt(log z)) is not small for
# large |v| (sqrt(log 1e4) ~ 3), so we test |v| <= 10, where the main term dominates. In the paper
# z = x^{1/logloglog x} and this error is negligible for all |v| <= 2H.
X = 2e7
worst = 0.0
for z in [1e3, 1e4]:
  U = np.log(X) / np.log(z)
  for Y in [2 * z, 1e5, 2e7]:
    sel = (P > z) & (P <= Y)
    for v in np.concatenate([np.geomspace(1e-4, 10, 200), -np.geomspace(1e-4, 10, 40)]):
        E = np.exp(-np.sum(np.log1p(-P[sel] ** (-1 - 1j * v))))
        Gam = min(1.0, 1 / (abs(v) * np.log(X)))
        worst = max(worst, abs(E - 1) / (U * Gam))
z = 1e3; U = np.log(X) / np.log(z)
check("|E_Y(v)-1| <= C U Gamma(v) (eq:Ebound) on the tested range", worst < 5, f"empirical C = {worst:.2f}")
check("|E_Y(v)| <= E_X(0) asymp U", np.exp(-np.sum(np.log1p(-1 / P[(P > z) & (P <= X)]))) / U < 1.5)

raise SystemExit(summary("check_prime_sums.py"))
