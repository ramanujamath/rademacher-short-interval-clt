"""Combinatorial and orthogonality facts, checked by exhaustive enumeration.

1. Lemma 8.2: squarefree quadruples (m1,m2,m3,m4) over a prime set with square product
   correspond bijectively to 7-tuples (k,g,h,a,b,c,d) of pairwise coprime squarefree
   integers via m1=kgab, m2=khac, m3=kgcd, m4=khbd; the weight is k^-2(ghabcd)^-1;
   the phases match the six-factor formula (eq:sixfactors); and when k=1 the
   constraints (eq:star) are equivalent to (eq:constraints).
2. Lemma 7.1: the Bonami 2-to-4 inequality ||A||_4^2 <= sum 3^|S| c_S^2.
3. Lemma 3.1: the Walsh expansion of int W_T |R_w|^2, the formula Var = M_w^2 sum_{d>1} c_d^2,
   and the vanishing of c_d unless some divisor ratio a/b lies in [e^{-2/T}, e^{2/T}].
"""
import itertools
import numpy as np
from common import check, summary, What_closed, W_T
from scipy.integrate import quad

print("check_combinatorics.py")
rng = np.random.default_rng(3)

# ---------- 1. six-variable parametrization ----------
P = [11, 13, 17, 19, 23]
subsets = [frozenset(c) for r in range(len(P) + 1) for c in itertools.combinations(P, r)]
prod = lambda S: int(np.prod(sorted(S))) if S else 1
quads = []
for m in itertools.product(subsets, repeat=4):
    if all(sum(p in mi for mi in m) % 2 == 0 for p in P):
        quads.append(m)
check("number of square-product quadruples = 8^|P| (each prime: none, or one of 7 patterns)",
      len(quads) == 8 ** len(P), f"{len(quads)} vs {8**len(P)}")

patterns = {"k": (0, 1, 2, 3), "g": (0, 2), "h": (1, 3), "a": (0, 1), "b": (0, 3), "c": (1, 2), "d": (2, 3)}
ok_param, ok_weight, ok_phase, ok_constr = True, True, True, True
t = rng.uniform(-5, 5, 4)
for m in quads:
    var = {v: frozenset(p for p in P if tuple(i for i in range(4) if p in m[i]) == pat) for v, pat in patterns.items()}
    # every prime that occurs is in exactly one variable
    used = set().union(*var.values())
    ok_param &= used == set().union(*m)
    k, g, h, a, b, c, d = (prod(var[v]) for v in "kghabcd")
    M = [prod(mi) for mi in m]
    ok_param &= M == [k * g * a * b, k * h * a * c, k * g * c * d, k * h * b * d]
    w = 1 / np.sqrt(np.prod(np.array(M, dtype=float)))
    ok_weight &= abs(w - 1 / (float(k) ** 2 * g * h * a * b * c * d)) < 1e-12 * w
    if k == 1:
        lhs = M[0] ** (-0.5 - 1j * t[0]) * M[1] ** (-0.5 - 1j * t[1]) * M[2] ** (-0.5 + 1j * t[2]) * M[3] ** (-0.5 + 1j * t[3])
        rhs = 1 / (g ** (1 + 1j * (t[0] - t[2])) * h ** (1 + 1j * (t[1] - t[3])) * a ** (1 + 1j * (t[0] + t[1]))
                   * b ** (1 + 1j * (t[0] - t[3])) * c ** (1 + 1j * (t[1] - t[2])) * d ** (1 - 1j * (t[2] + t[3])))
        ok_phase &= abs(lhs - rhs) < 1e-12 * abs(rhs)
        Pm = lambda S: max(S) if S else 1
        star = (all(mi for mi in m) and Pm(m[0]) == Pm(m[2]) and Pm(m[1]) == Pm(m[3])
                and m[0] != m[2] and m[1] != m[3])
        abcd = var["a"] | var["b"] | var["c"] | var["d"]
        constr = (g > 1 and h > 1 and len(abcd) > 0 and Pm(abcd) < min(Pm(var["g"]), Pm(var["h"])))
        ok_constr &= star == constr
check("unique representation m1=kgab, m2=khac, m3=kgcd, m4=khbd", ok_param)
check("absolute weight |m1m2m3m4|^{-1/2} = k^-2 (ghabcd)^-1", ok_weight)
check("phases agree with the six-factor formula (eq:sixfactors)", ok_phase)
check("for k=1, constraints (eq:star) <=> (eq:constraints)", ok_constr)

# ---------- 2. Bonami inequality ----------
n = 8
signs = np.array(list(itertools.product([-1, 1], repeat=n)))
S_all = [S for r in range(n + 1) for S in itertools.combinations(range(n), r)]
chars = np.array([np.prod(signs[:, list(S)], axis=1) if S else np.ones(len(signs)) for S in S_all]).T
worst = 0.0
for trial in range(300):
    c = rng.standard_normal(len(S_all)) * (rng.random(len(S_all)) < rng.uniform(0.05, 1))
    A = chars @ c
    lhs = np.sqrt(np.mean(A ** 4))
    rhs = sum(3 ** len(S) * ci ** 2 for S, ci in zip(S_all, c))
    if rhs > 0:
        worst = max(worst, lhs / rhs)
check("||A||_4^2 <= sum 3^|S| c_S^2 (300 random polynomials, 8 signs)", worst <= 1 + 1e-12, f"max ratio {worst:.4f}")
# the d_3 identity mu^2 3^omega = d_3 on squarefree integers
def d3(m):
    return sum(1 for a in range(1, m + 1) if m % a == 0 for b in range(1, m // a + 1) if (m // a) % b == 0)
check("mu(m)^2 3^omega(m) = d_3(m) for squarefree m", all(d3(m) == 3 ** len(f) for m, f in
      [(1, []), (6, [2, 3]), (30, [2, 3, 5]), (210, [2, 3, 5, 7]), (11, [11])]))

# ---------- 3. Lemma 3.1 Walsh expansion and exact variance ----------
pw = [2, 3, 5, 7, 11, 13]
Mw = np.prod([1 + 1 / p for p in pw])
for T in [2.0, 6.0, 30.0]:
    cd = {}
    for r in range(len(pw) + 1):
        for D in itertools.combinations(pw, r):
            d = int(np.prod(D)) if D else 1
            s = 0.0
            for rr in range(len(D) + 1):
                for A in itertools.combinations(D, rr):
                    a_ = int(np.prod(A)) if A else 1
                    s += What_closed(np.log(a_ / (d // a_)), T)
            cd[D] = np.sqrt(d) / np.prod([p + 1 for p in D]) * s
    # exact integral for every sign pattern: int W |R_w|^2 = sum_{m,n} f f (mn)^-1/2 hatW(log(n/m))
    sm = [int(np.prod(C)) if C else 1 for r in range(len(pw) + 1) for C in itertools.combinations(pw, r)]
    arr = np.array(sm, float)
    H = What_closed(np.log(arr[None, :] / arr[:, None]), T) / np.sqrt(np.outer(arr, arr))
    vals, walsh = [], []
    for eps in itertools.product([-1, 1], repeat=len(pw)):
        e = dict(zip(pw, eps))
        fv = np.array([np.prod([e[p] for p in pw if m % p == 0]) for m in sm], float)
        vals.append(fv @ H @ fv)
        walsh.append(Mw * sum(c * np.prod([e[p] for p in D]) for D, c in cd.items()))
    vals, walsh = np.array(vals), np.array(walsh)
    check(f"Walsh expansion reproduces int W_T|R_w|^2 exactly (T={T:g})", np.max(abs(vals - walsh)) < 1e-9 * np.max(abs(vals)))
    check(f"mean = M_w J_T (T={T:g})", abs(vals.mean() - Mw * 4 * np.pi * T * T * np.sinh(1 / T)) < 1e-9 * vals.mean())
    var_formula = Mw ** 2 * sum(c * c for D, c in cd.items() if D)
    check(f"Var = M_w^2 sum_(d>1) c_d^2 (T={T:g})", abs(vals.var() - var_formula) < 1e-8 * max(1, var_formula),
          f"Var={vals.var():.4e}")
    # support: c_d = 0 unless some a|d has |log(a/(d/a))| <= 2/T
    ok = True
    for D, c in cd.items():
        d = int(np.prod(D)) if D else 1
        near = any(abs(np.log((int(np.prod(A)) if A else 1) ** 2 / d)) <= 2 / T
                   for rr in range(len(D) + 1) for A in itertools.combinations(D, rr))
        if not near:
            ok &= abs(c) < 1e-12
    check(f"c_d vanishes unless a divisor ratio lies within e^(+-2/T) (T={T:g})", ok)

# one direct numerical check of the integral for a single sign pattern
fv = np.array([np.prod([(-1 if p in (3, 7) else 1) for p in pw if m % p == 0]) for m in sm], float)
T = 2.0
direct = 2 * sum(quad(lambda t: W_T(t, T) * abs(np.sum(fv * arr ** (-0.5 - 1j * t))) ** 2, k, k + 2, limit=200)[0]
                 for k in range(0, 3000, 2))
H = What_closed(np.log(arr[None, :] / arr[:, None]), T) / np.sqrt(np.outer(arr, arr))
check("closed-form integral matches direct quadrature", abs(direct - fv @ H @ fv) < 0.01 * (fv @ H @ fv),
      f"{direct:.4f} vs {fv @ H @ fv:.4f}")

raise SystemExit(summary("check_combinatorics.py"))
