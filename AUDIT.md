# Audit of "Gaussian limits for Rademacher random multiplicative sums in logarithmic short intervals"

This audit covers the submitted version (`paper/original_v1.pdf`, 261102-JhaSumit v1) and the rewrite in `paper/rademacher_short_interval_clt.tex`. Every step of the proof was re-derived by hand. Where a step reduces to an identity or an explicit computation, it was also checked numerically (see `verification/`). Section and equation references below are to the **rewritten** manuscript unless marked "(v1)".

## Verdict

No mathematical gap was found. The argument works as written. Every inequality was checked, as was the exponent bookkeeping behind each `o(1)` claim. The external inputs were checked against their sources:

| Input | Source checked | Status |
|---|---|---|
| Ψ(a+h,z) − Ψ(a,z) ≤ Ψ(h,z) for a, h > z, z large | Granville, *Smooth numbers*, §1.4 (Hildebrand 1985) | matches |
| Ψ(h,z) = h u^{−u+o(u)} for z ≥ (log h)^{1+ε} | Granville (1.12) (Canfield–Erdős–Pomerance) | matches |
| Sieve estimate with relative error (e/(s log s))^s, s < log M/(log log 3M)^6 | Harper–Soundararajan–Xu, arXiv:2606.29040, Lemma 2.2 (after Iwaniec) | matches; HSX use ∏_{p<z}, v1 uses ∏_{p≤z} (relative difference O(1/z), absorbed; now noted in the text) |
| Soundararajan–Xu range y ≤ x/(log x)^{log 4−1+ε}; HSX §1.4 obstacle (translation invariance in Prop. 3.7) | HSX §1 and §1.4 | matches |
| Durkan–Pearce-Crump arXiv:2607.29429 | arXiv listing | exists |

## Corrections made in the rewrite

These are presentation-level fixes. None changes a statement or the structure of the proof.

1. **Typo in Lemma 9.2 (v1 p. 16).** The opposite-sign bound was written as O(δL³R√(log z)/log X). The factor R = L⁶⁰/log z should not be there; the correct bound is O(δL³√(log z)/log X). The conclusion (< δL^{−10}) is unaffected.
2. **Lemma 7.1 (v1 p. 10).** The d₃ bound was garbled in typesetting. It now reads: μ(m)²3^{ω(m)} = d₃(m) for squarefree m, and Σ_{m≤Y} d₃(m) ≪ Y(log 2Y)².
3. **Removing μ(m)² (Prop. 6.1).** "Orthogonality and unique decomposition" is now spelled out: E(V_f* − V_f) counts the integers in (x, X] divisible by p² for some p > z.
4. **Notation clashes**, all resolved by renaming:
   - `H` meant both T² (§4–5) and (log X)^{1000} (§8–10); the first is now `H_0`.
   - `R` meant both the Euler product R_z and L⁶⁰/log z; the latter is now `r_2` (and `r` is `r_1`).
   - `b` meant both the barrier parameter √ℓ and a summation variable; the barrier parameter is now `β` (and `b₀` is `β₀`).
   - `a` meant p^{−1/2}, a(t), a_m and a summation variable; a(t) is now `κ(t)`, and p^{−1/2} is the local `γ`.
   - `D` meant D_z, D_p and D_±; these are now `Θ_z`, `D_p` and `𝒟_±`.
   - `J` meant both J_T and the random integral J; the latter is now `𝒥`, `𝒥_{r_1}`, `𝒞_{r_1}`.
   - `K` meant both log log x and the kernel K(t); log log x is replaced by L = log log z. This is valid because K^{−A} ≤ L^{−A} and only some A > 0 is needed.
   - `w` was reused for |u| in Lemma 9.2; `I` for both the integral and the intervals I_r; `Q_t` (tilted measure) next to `Q_x`; `m` both a vector and M_zJ_T/2; `h` both a function and a length; `c_r` next to the constant c.
5. **Undefined objects are now defined:** d₃, ϑ, Euler's constant γ in the divisor-sum asymptotic, the Fourier/Plancherel convention, E_{>z}, "Walsh characters" (now described as an orthonormal family of sign products), "predictable quadratic variation", "frozen coefficient", 𝓔_Y with its logarithm, and the meaning of the barrier event in words.
6. **Short justifications added:**
   - the kernel bound W_T ≪ min{1, T²/(1+t²)}, via 2(cosh(1/T)−1) + 4sin²(t/T);
   - a ≥ T/2 in Lemma 3.1;
   - why |t₂−t₃|, |t₃+t₄| ≥ r₂/2 outside 𝒟_±;
   - how B̃_p decomposes over m₁ with P(m₁) = p.
7. **Prose.** An outline of the proof was added to the introduction. The abstract no longer uses jargon such as "barriers with negligible loss under a tilted measure" without explanation. Section and lemma statements are self-contained.

## Observations for the author (no change made)

- **The hypothesis α < log 4 − 1 does not appear to be used.** The only upper bound on α used in the argument is α < 1: in Prop. 6.1 it gives x/(y log z) → 0 and √(x/y)/log z → 0. The lower bound α > 0 is used in Lemma 9.2 (δ²L ≪ δL^{−10}) and §10 (log T ∼ αL). If this is right, the theorem holds for every fixed 0 < α < 1, which overlaps the Soundararajan–Xu range. Either state the wider range, or say that the restriction only marks the range not previously covered.
- **The proof is purely asymptotic.** `verification/check_asymptotic_bookkeeping.py` sets all implied constants to 1 and finds where each `o(1)` ratio drops below 1, in terms of L₀ = log log x. Several thresholds are astronomically large: about 2·10¹¹ for the bad-frequency term U⁶L^{−3/4}, and log L₀ ≈ 7·10⁵ for Ψ(h,z) ≪ hL^{−1} (Lemma 6.2(a)). This is normal for this kind of argument and is not an error. It does mean that numerical experiments cannot probe the regime of the proof.
- **Monte Carlo** at x ≤ 10⁶ gives normalized sums whose variance is correct to within 3% and whose distribution passes a Kolmogorov–Smirnov test against N(0,1). The excess kurtosis is visibly positive (0.5–4), which is consistent with slow, log-log-scale convergence driven by the random Euler product. This is a sanity check only.
- **Acknowledgements.** The rewrite adds one sentence saying that Claude (Anthropic) helped with the expository revision and the numerical checks. Edit or remove it as you see fit.

## Numerical verification summary

Run `python3 verification/run_all.py`. All 124 checks pass.

| Script | What it checks | Checks |
|---|---|---|
| `check_kernel.py` | W_T as squared modulus; J_T = 4πT² sinh(1/T); pointwise bound; Fourier transform of W_T in closed form with compact support [−2/T, 2/T] (computed independently by quadrature) | 9 |
| `check_local_moments.py` | Exact one-prime identities (eq:4mom), (eq:low4); tilted mean and variance (eq:tiltmean); two-point expansion in Lemma 9.1; Taylor bound in Lemma 7.2; relaxation factor in Lemma 8.2; Walsh expansion | 13 |
| `check_prime_sums.py` | Oscillatory prime sums (eq:pnt), Mertens (eq:mertens), Euler-product bound (eq:Ebound), on a finite range | 4 |
| `check_combinatorics.py` | Six-variable parametrization of square-product quadruples (bijection, weight, phases, equivalence of the constraints); Bonami inequality; Lemma 3.1 Walsh coefficients, exact variance, and support of c_d | 20 |
| `check_transfer_identities.py` | Parseval identity of Prop. 6.1 (exact left side vs. closed-form and quadrature right sides); normalization to cy; smooth/rough and martingale decompositions; tower identity E(Q_x \| F_z) = V_f, exact | 10 |
| `check_asymptotic_bookkeeping.py` | Every `o(1)` exponent comparison in §§4–10, for α = 0.05 and 0.38 | 56 |
| `monte_carlo_clt.py` | Distribution of the normalized sum at x = 10⁵, 10⁶ and α = 0.2, 0.35 | 12 |
