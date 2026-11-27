"""Monte Carlo sanity check of Theorem 1.1 at computationally accessible sizes.

For x in {1e5, 1e6} and alpha in {0.2, 0.35, 0.8} we sample many independent sign sequences
and compute S = sum_{x<n<=x+y} f(n) with y = x/(log x)^alpha. We report:
  * the exact variance ratio  #{squarefree n in (x,x+y]} / (c y),
  * mean, variance, skewness and excess kurtosis of S/sqrt(cy),
  * the Kolmogorov-Smirnov distance to N(0,1).

This cannot test the asymptotic statement (log log x is tiny at these sizes, and the
proof's error terms are only small for astronomically large x). It checks that nothing is
grossly wrong with the normalization and that the distribution looks close to Gaussian.
Usage:  python3 monte_carlo_clt.py [n_samples]
"""
import sys
import numpy as np
from scipy import sparse, stats
from common import check, summary, primes_upto

C = 6 / np.pi ** 2
nsamp = int(sys.argv[1]) if len(sys.argv) > 1 else 4000
rng = np.random.default_rng(2026)
print(f"monte_carlo_clt.py  (samples per case = {nsamp})")


def interval_incidence(x, X):
    """Sparse 0/1 matrix rows = squarefree n in (x,X], cols = primes, entry 1 if p|n."""
    n0, n1 = int(np.floor(x)) + 1, int(np.floor(X))
    L = n1 - n0 + 1
    primes = primes_upto(n1)
    sqf = np.ones(L, dtype=bool)
    rows, cols = [], []
    for j, p in enumerate(primes):
        p = int(p)
        start = (-n0) % p
        idx = np.arange(start, L, p)
        if idx.size == 0:
            continue
        if p * p <= n1:
            s2 = (-n0) % (p * p)
            sqf[s2::p * p] = False
        rows.append(idx)
        cols.append(np.full(idx.size, j, dtype=np.int32))
    rows = np.concatenate(rows); cols = np.concatenate(cols)
    keep = sqf[rows]
    A = sparse.csr_matrix((np.ones(keep.sum(), dtype=np.float32), (rows[keep], cols[keep])), shape=(L, len(primes)))
    A = A[np.nonzero(sqf)[0]]
    return A, len(primes), int(sqf.sum())


for x in [1e5, 1e6]:
    for alpha in [0.2, 0.35, 0.8]:
        y = x / np.log(x) ** alpha
        A, npr, nsq = interval_incidence(x, x + y)
        S = np.empty(nsamp)
        batch = 250
        for b0 in range(0, nsamp, batch):
            k = min(batch, nsamp - b0)
            bits = rng.integers(0, 2, size=(npr, k)).astype(np.float32)  # eps_p = 1 - 2*bit
            parity = np.asarray(A @ bits) % 2                              # number of -1 factors mod 2
            S[b0:b0 + k] = (1 - 2 * parity).sum(axis=0)
        Z = S / np.sqrt(C * y)
        ks = stats.kstest(S / np.sqrt(nsq), "norm").statistic
        print(f"  x={x:.0e}, alpha={alpha}: y={y:.0f}, #sqfree/(cy)={nsq/(C*y):.4f}, mean={Z.mean():+.3f}, "
              f"var={Z.var():.3f}, skew={stats.skew(Z):+.3f}, exkurt={stats.kurtosis(Z):+.3f}, KS={ks:.3f}")
        check(f"exact variance ratio close to 1 (x={x:.0e}, alpha={alpha})", abs(nsq / (C * y) - 1) < 0.01)
        check(f"sample variance of S/sqrt(#sqfree) within 10% of 1", abs((S / np.sqrt(nsq)).var() - 1) < 0.1)
        # KS critical value at level 0.001 is about 1.95/sqrt(n)
        check(f"KS distance to N(0,1) below 0.001-level critical value", ks < 1.95 / np.sqrt(nsamp), f"crit {1.95/np.sqrt(nsamp):.3f}")

raise SystemExit(summary("monte_carlo_clt.py"))
