"""Checks for the interval kernel W_T (eq:WT, eq:Wbounds, eq:What).

1. W_T(t)/T^2 = |e^{s/T} - e^{-s/T}|^2 / |s|^2 with s = 1/2 + it.
2. J_T = int W_T = 4 pi T^2 sinh(1/T).
3. W_T(t) <= C min(1, T^2/(1+t^2)) with a modest absolute C.
4. The Fourier transform of W_T equals 4 pi T^2 sinh(1/T - |xi|/2) on |xi| <= 2/T
   and vanishes for |xi| > 2/T.
5. The cosine integral int cos(bt)/(1/4+t^2) dt = 2 pi e^{-|b|/2}.
"""
import numpy as np
from scipy.integrate import quad
from common import check, summary, W_T, What_closed

print("check_kernel.py")

# 1. modulus identity
rng = np.random.default_rng(1)
worst = 0.0
for T in [2.0, 5.0, 37.0, 1000.0]:
    t = rng.uniform(-200, 200, 2000)
    s = 0.5 + 1j * t
    lhs = W_T(t, T) / T**2
    rhs = np.abs(np.exp(s / T) - np.exp(-s / T)) ** 2 / np.abs(s) ** 2
    worst = max(worst, np.max(np.abs(lhs - rhs) / np.maximum(rhs, 1e-300)))
check("W_T/T^2 = |e^{s/T}-e^{-s/T}|^2/|s|^2", worst < 1e-9, f"max rel err {worst:.2e}")

# 5. cosine integral (Fourier transform of the Cauchy-type density)
worst = 0.0
for b in [0.0, 0.3, 1.0, 4.0]:
    val = 2 * quad(lambda t: 1.0 / (0.25 + t * t), 0, np.inf, weight="cos", wvar=b)[0] if b > 0 else \
        2 * quad(lambda t: 1.0 / (0.25 + t * t), 0, np.inf)[0]
    worst = max(worst, abs(val - 2 * np.pi * np.exp(-abs(b) / 2)))
check("int cos(bt)/(1/4+t^2) dt = 2 pi e^{-|b|/2}", worst < 1e-7, f"max abs err {worst:.2e}")


# 2. J_T. Write W_T = T^2 [2cosh(1/T) - 2cos(2t/T)]/(1/4+t^2) and integrate the
# two pieces with QUADPACK's Fourier-weighted routine (robust for oscillatory tails).
def J_numeric(T):
    a = 2 * np.cosh(1 / T) * 2 * quad(lambda t: 1 / (0.25 + t * t), 0, np.inf)[0]
    b = 2 * 2 * quad(lambda t: 1 / (0.25 + t * t), 0, np.inf, weight="cos", wvar=2 / T)[0]
    return T * T * (a - b)


worst = 0.0
for T in [2.0, 3.0, 10.0, 100.0, 1e4]:
    Jn = J_numeric(T)
    Jc = 4 * np.pi * T * T * np.sinh(1 / T)
    worst = max(worst, abs(Jn - Jc) / Jc)
check("J_T = 4 pi T^2 sinh(1/T)", worst < 1e-7, f"max rel err {worst:.2e}")
check("J_T ~ 4 pi T", abs(4 * np.pi * 1e6**2 * np.sinh(1e-6) / (4 * np.pi * 1e6) - 1) < 1e-10)

# 3. pointwise bound; report the best constant on a fine grid
Cmax = 0.0
for T in [2.0, 7.0, 50.0, 1e3, 1e5]:
    t = np.concatenate([np.linspace(0, 10 * T, 200001), np.geomspace(10 * T, 1e4 * T, 20001)])
    ratio = W_T(t, T) / np.minimum(1.0, T * T / (1 + t * t))
    Cmax = max(Cmax, ratio.max())
check("W_T(t) <= C min{1, T^2/(1+t^2)}", Cmax < 10, f"empirical C = {Cmax:.3f}")
check("W_T >= 0 and even", np.all(W_T(np.linspace(-50, 50, 10001), 3.0) >= 0)
      and np.allclose(W_T(np.linspace(0, 50, 101), 3.0), W_T(-np.linspace(0, 50, 101), 3.0)))


# 4. Fourier transform: hat W(xi) = int e^{i t xi} W_T(t) dt.  Using linearity,
# hat W = T^2 [2cosh(1/T) F(xi) - F(xi + 2/T) - F(xi - 2/T)], with F(b) the cosine
# integral of (1/4+t^2)^{-1}.  We evaluate F numerically (not by its closed form)
# so that this is an independent check of (eq:What).
def F_num(b):
    b = abs(b)
    if b < 1e-9:  # QAWF is unreliable for a vanishing frequency; F is continuous there
        return 2 * quad(lambda t: 1 / (0.25 + t * t), 0, np.inf)[0]
    return 2 * quad(lambda t: 1 / (0.25 + t * t), 0, np.inf, weight="cos", wvar=b)[0]


worst_in, worst_out = 0.0, 0.0
# Large T makes the three-term combination cancel catastrophically in floating point,
# so we test moderate T only (the identity is algebraic in T).
for T in [2.0, 3.0, 5.0, 8.0]:
    for xi in np.linspace(-3 / T, 3 / T, 61):
        num = T * T * (2 * np.cosh(1 / T) * F_num(xi) - F_num(xi + 2 / T) - F_num(xi - 2 / T))
        cl = What_closed(xi, T)
        if abs(xi) <= 2 / T:
            worst_in = max(worst_in, abs(num - cl) / (T))
        else:
            worst_out = max(worst_out, abs(num) / T)
check("hat W_T = 4 pi T^2 sinh(1/T-|xi|/2) on |xi|<=2/T", worst_in < 1e-6, f"err/T {worst_in:.2e}")
check("hat W_T = 0 for |xi| > 2/T (compact support)", worst_out < 1e-6, f"err/T {worst_out:.2e}")
check("0 <= hat W_T <= J_T", all(0 <= What_closed(x, 3.0) <= 4 * np.pi * 9 * np.sinh(1 / 3) + 1e-12
                                for x in np.linspace(-1, 1, 101)))

raise SystemExit(summary("check_kernel.py"))
