"""Shared helpers for the verification scripts (numpy/scipy only)."""
import numpy as np

PASS = 0
FAIL = 0


def check(name, ok, detail=""):
    global PASS, FAIL
    if ok:
        PASS += 1
        print(f"  [PASS] {name}" + (f"  ({detail})" if detail else ""))
    else:
        FAIL += 1
        print(f"  [FAIL] {name}" + (f"  ({detail})" if detail else ""))
    return ok


def summary(script):
    print(f"{script}: {PASS} passed, {FAIL} failed")
    return 1 if FAIL else 0


def primes_upto(n):
    s = np.ones(n + 1, dtype=bool)
    s[:2] = False
    for p in range(2, int(n ** 0.5) + 1):
        if s[p]:
            s[p * p :: p] = False
    return np.nonzero(s)[0]


def W_T(t, T):
    """Kernel (1.2): T^2 (2cosh(1/T) - 2cos(2t/T)) / (1/4 + t^2)."""
    return T * T * (2 * np.cosh(1 / T) - 2 * np.cos(2 * t / T)) / (0.25 + t * t)


def What_closed(xi, T):
    """Closed form (3.2) of the Fourier transform of W_T."""
    xi = np.abs(xi)
    out = 4 * np.pi * T * T * np.sinh(1 / T - xi / 2)
    return np.where(xi <= 2 / T, out, 0.0)
