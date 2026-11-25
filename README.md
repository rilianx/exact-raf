# Exact multiplication of reduced affine forms

Code, data and scripts of the paper

> I. Araya, D. Miranda, V. Reyes, J. Geraldo.
> *Exact multiplication of reduced affine forms.* Submitted to Numerical Algorithms.

A reduced affine form `x0 + x1 e1 + ... + xn en + xr ex` keeps a fixed number of noise symbols
and accumulates all the approximation errors in a single error term. The product of two such forms
keeps the linear part `x0 y0 + e^T(y0 x + x0 y)` and must enclose the nonlinear part

    N = x0 yr ey + y0 xr ex + (e^T x + xr ex)(e^T y + yr ey)

by an interval. The multiplication of Skalna and Hladik (SH, Numer. Algorithms 76, 2017) computes
the exact range of the quadratic part `(e^T x + xr ex)(e^T y + yr ey)` but bounds
`x0 yr ey + y0 xr ex` separately. The paper gives a multiplication **E** that computes the exact
range of the whole of `N` in O(n log n) time, by a single walk along the boundary of a planar
zonotope, and shows that E coincides with SH when one operand has no error or a zero center, and
that `wid SH <= 3/2 wid E` otherwise (tight). S denotes the standard multiplication of the affine
arithmetic of Ibex (fAF2).

## Repository layout

```
raf/                    Python implementation (double precision, no outward rounding): reduced affine
                        forms written in the variables, set_product('faf2' | 'skalna' | 'exact')
                        selects S, SH or E; Chebyshev linearization of squares and unary functions
experiments/            range-enclosure experiments (Section 6.1) and the scripts of Tables 1 and 2
tests/test_soundness.py sampling checks of the three multiplications
results/                JSON results and the LaTeX tables included in the paper
ibex/                   rigorous implementation in Ibex and solver experiments (Sections 6.2-6.3,
                        Tables 3 and 4), see ibex/README.md
benchmarks/ibex/        synthetic IbexOpt instances (synq, synt, synp) and their generators
paper/numa/             the paper, its figures (fig/), the flat source folder for the journal
                        (make_submission.py) and a small-format version for e-readers (kindle/)
```

## Requirements

- Python 3.9 or later with `numpy` and `matplotlib` (`pip install -r requirements.txt`).
- For Table 2: [IbexOpt](https://github.com/ibex-team/ibex-lib), which computes the exact ranges
  (`export IBEXOPT=/path/to/ibexopt`).
- For the Ibex implementation and the solver experiments: a C++ compiler, CMake and an LP solver
  (CLP or SoPlex); `ibex/build.sh` fetches, patches and builds everything.

All the results used in the paper are stored in `results/` and `ibex/results/`, so that every table
can be regenerated without running the experiments again.

## Reproducing the results of the paper

From the root of the repository:

| Result | Experiment | Table from the stored results |
|---|---|---|
| Table 1 (random products) | `python experiments/random_products.py` | `python experiments/make_table_random.py` |
| Table 2 (test functions of Skalna and Hladik) | `IBEXOPT=... python experiments/sh_benchmark.py` | `python experiments/make_table_sh.py` |
| Corollary 1 (attainment of the bound) | `python experiments/check_gain_bound.py` | |
| Table 3 (time per product) | `bench_mult` (built by `ibex/build.sh`) | `python ibex/analysis/make_table_cost.py` |
| Table 4 (IbexOpt) | `ibex/run_interleaved.sh`, see `ibex/README.md` | `python ibex/analysis/make_table_solver_interleaved.py` |
| Figures 1 and 2 | `python paper/numa/fig/make_fig_zonotope.py`, `python paper/numa/fig/make_fig_square.py` | |

Checks: `python tests/test_soundness.py` (Python) and `affine_containment_check` (Ibex, see
`ibex/README.md`).

The three multiplications on the same representation:

```python
from raf import set_box, raf_variables, set_product

for rule in ["faf2", "skalna", "exact"]:        # S, SH, E
    set_box([(0.0, 2.0), (1.0, 3.0)])
    set_product(rule)
    x1, x2 = raf_variables()
    f = (x1 + 2 * x1 * x2) * (3 * x2 - x1 * x2)
    print(rule, f.get_range(), f.residual_width())
# faf2   (-62.0, 96.0)    102.0
# skalna (-46.0, 74.225)  64.225
# exact  (-36.0, 74.225)  54.225
```

## Using or porting the multiplications

The three multiplications only change the interval that encloses the nonlinear part of the product.
For `x = x0 + sum_i a_i e_i + xr ex` and `y = y0 + sum_i b_i e_i + yr ey`, every rule keeps the linear
part `x0 y0 + sum_i (x0 b_i + y0 a_i) e_i`, computes an interval `[lo, hi]` that encloses

    N = x0 yr ey + y0 xr ex + (a.e + xr ex)(b.e + yr ey),

and returns the form with center `x0 y0 + (lo + hi)/2` and error term `(hi - lo)/2`.

- **Reference (Python, self-contained):** `raf.nonlinear_enclosure(rule, x0, a, xr, y0, b, yr)` in
  `raf/forms.py` takes plain floats and lists and returns `(lo, hi)` for `rule` = `'faf2'` (S),
  `'skalna'` (SH) or `'exact'` (E). E is Algorithm 1 of the paper: the generators `(a_i, b_i)` are
  turned to the half-plane `u > 0` and sorted by angle once (`_sorted_generators`), and the boundary
  of the zonotope is walked for the two shifts `(xr, yr)` and `(xr, -yr)` (`_uv_range_sorted`), checking
  the vertices and the stationary point of each edge. About 70 lines without dependencies.
- **Rigorous (C++, floating point with outward rounding):** `ibex/fast/ibex_ExactMult_fast.h`,
  function `raf_mult_enclosure_fast(a, b, n, x0, ex, y0, ey, exact)`, which returns an
  `ibex::Interval` containing the range of `N` (`exact = false` gives SH). It depends only on
  `ibex::Interval` for the result and can be adapted to any interval type; the rounding errors of the
  walk in `double` are covered by an a priori bound (Section 5 of the paper, `ibex/fast/README.md`).
  Its use in Ibex is in `operator*=` of `ibex_Affine2_fAF2.cpp` (patches 0016 and 0017, `ibex/patches/`):
  the enclosure replaces the nonlinear error terms of fAF2, and its midpoint is added to the center.

To port E to another affine arithmetic library, compute the linear part as the library already does,
call the enclosure of `N` above with the center, noise coefficients and error term of each operand,
and add the midpoint and the radius of the interval to the center and to the new error term, with
the library's own handling of rounding errors.

## License

See `LICENSE`.
