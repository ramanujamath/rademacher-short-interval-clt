# Gaussian limits for Rademacher random multiplicative sums in logarithmic short intervals

**Author:** Sumit Kumar Jha (Independent Researcher, Greater Noida, India). Contact: sumitkumarjha75@gmail.com

Let *f* be a Rademacher random multiplicative function. It is supported on the squarefree integers, and its values at the primes are independent fair ±1 signs. Fix 0 < α < log 4 − 1 and put y = x/(log x)^α. The main theorem says that

$$\frac{1}{\sqrt{6y/\pi^2}}\sum_{x<n\le x+y} f(n) \;\Rightarrow\; N(0,1)\qquad (x\to\infty).$$

The proof follows the martingale and conditional-variance method of Harper, Soundararajan and Xu. They use the fact that in the Steinhaus model the law of the random Euler product does not change when the frequency is shifted. That fails for Rademacher signs, and the paper replaces it with:

* a **weighted concentration theorem** for ∫ W_T(t)|R_z(t)|² dt, proved by a barrier argument under a change of measure, at each fixed frequency;
* a **conditional fourth-moment estimate** that controls correlations between both equal and opposite frequencies. A perfect-square quadruple of real characters is described by six free variables rather than four.

MSC 2020: 11K65 (primary); 11N37, 60F05.

## Repository layout

| Path | Contents |
|---|---|
| `paper/rademacher_short_interval_clt.tex` | Rewritten manuscript (LaTeX source) |
| `paper/rademacher_short_interval_clt.pdf` | Compiled manuscript |
| `paper/original_v1.pdf` | Version as submitted (v1), kept for reference |
| `AUDIT.md` | Step-by-step audit: verdict, corrections made, observations for the author |
| `verification/` | Numerical sanity checks (Python, numpy + scipy) |

## Building the paper

```bash
cd paper
latexmk -pdf rademacher_short_interval_clt.tex
```

## Running the checks

```bash
cd verification
python3 run_all.py           # about 3 minutes; all 124 checks should pass
python3 run_all.py --quick   # smaller Monte Carlo sample
```

The scripts check exact identities, both symbolically and by exhaustive enumeration over signs. They also check standard analytic estimates on finite ranges, the exponent bookkeeping behind every `o(1)` claim, and the distribution of the normalized sum by simulation at x ≤ 10⁶. They are sanity checks and do not replace the proof. In particular, the regime of the theorem is far beyond what can be simulated; see `AUDIT.md`.

## References

* S. Chatterjee, K. Soundararajan, *Random multiplicative functions in short intervals*, IMRN 2012.
* K. Soundararajan, M. W. Xu, *Central limit theorems for random multiplicative functions*, J. Anal. Math. 151 (2023).
* A. J. Harper, K. Soundararajan, M. W. Xu, *Distribution of random multiplicative functions in short intervals, with proper normalization*, arXiv:2606.29040.
