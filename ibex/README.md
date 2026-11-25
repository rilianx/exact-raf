# Implementation in Ibex and solver experiments

Implementation of the multiplications SH (Skalna and Hladík, 2017) and E (this paper) in the affine
arithmetic `AF_fAF2` of [ibex-affine](https://github.com/ibex-team/ibex-affine), on top of the RAA-K
patch chain of `JCGeraldo/ibex-lib` (branch `raak-experimentos`, folder `raak/`), and the IbexOpt
experiments of Sections 6.2 and 6.3 of the paper.

## Build

```
ROOT=$HOME/ibex-local LP_LIB=clp bash ibex/build.sh
```

`build.sh` clones Ibex and ibex-affine, applies the RAA-K patches (0001–0015) and the patches below,
and builds the command-line optimizer `ibexopt-raak`, whose option `--arith=af2|sh|exact` selects
S, SH or E.

| Patch | Content |
|---|---|
| `0016-ibex-affine-multiplicacion-exacta.patch` | Modes `AFFINE_SH` and `AFFINE_EXACT` (`ibex_ExactMult.h`): rigorous enclosure of the whole nonlinear part `N`. |
| `0017-ibex-affine-multiplicacion-rapida.patch` | Fast version (`ibex_ExactMultFast.h`): walk in double precision with an a priori bound of the rounding errors and a filtered angular predicate (see `fast/`). |
| `0018-ibex-affine-cuadrado-chebyshev-af2.patch` | `IBEX_AF2_SQR=cheb`: S with the Chebyshev square (S⁺ in the paper). |
| `0019-ibex-affine-dos-esquinas-una-pasada.patch` | The two corners of Algorithm 1 in a single pass. |
| `0020-ibex-affine-cuadrado-exacto.patch` | `IBEX_SQR_EXACT=1`: exact range of the square (Proposition 4). |
| `0021-ibex-affine-diagnostico-s-mas-ancho.patch` | Diagnostic counters (no effect on the results): products in which S is wider than E. |

Validation: `tests/affine_containment_check.cpp` (containment of 20 000 random functions per family,
in interval arithmetic), `fast/verify.py` (fast walk against rational arithmetic),
`tests/test_sqexact.cpp` (exact square) and `tests/bench_mult.cpp` (time per product, Table 3).

## Solver experiments (Table 4)

All the variants are run in the same session, one after the other for each instance:

```
TASKS=tasks_all1.txt OUT=$HOME/ibex-local/results/interleaved_all P=4 bash ibex/run_interleaved.sh
```

Each line of a task file is `set file.bch time_limit variant`, with variants `af2` (S), `splus` (S⁺),
`sh` (SH) and `exact` (E). The synthetic instances are in `benchmarks/ibex/` (generators in
`benchmarks/ibex/generators/` and `analysis/gen_synp.py`); the COCONUT set is the optimization
benchmark distributed with Ibex (`benchs/optim`: easy, medium, hard, others, unsolved, blowup and
coconutbenchmark-library2; list in `results/coconut_manifest.csv`).

| File | Content |
|---|---|
| `results/table5_interleaved.csv` | Runs on synq, synt and synp |
| `results/table5_interleaved_all.csv` | Runs on the 298 COCONUT instances (60 s) |
| `results/sdiag/` | Counters of patch 0021 (5 s with E) used to select the 110 COCONUT instances of Table 4 |
| `analysis/make_table_solver_interleaved.py` | `results/table_solver.tex` (Table 4) |
| `analysis/summary_all_by_set.py`, `analysis/summary_all_subsets.py` | COCONUT results by subset |
| `results/bench_mult*.csv`, `results/containment_*` | Cost per product and containment tests |
