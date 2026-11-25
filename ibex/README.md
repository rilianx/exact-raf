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

Validation and cost (programs built by `build.sh` in `$ROOT`):

```
$ROOT/affine_containment_check 20000 11     # containment of 20 000 random functions per family, in intervals
$ROOT/bench_mult                            # time per product (Table 3): results/bench_mult_0019.csv
python ibex/fast/verify.py                  # fast walk of patch 0017 against rational arithmetic
```

## Solver experiments (Table 4)

The COCONUT set is the optimization benchmark distributed with Ibex (`ibex-lib/benchs/optim`: easy,
medium, hard, others, unsolved, blowup and coconutbenchmark-library2); copy its 298 distinct
instances, listed in `results/coconut_manifest.csv`, into one folder `$COCONUT`. The synthetic
instances are in `benchmarks/ibex/` (generators in `benchmarks/ibex/generators/`).

```
python ibex/make_tasks.py $COCONUT          # ibex/tasks_synthetic.txt, ibex/tasks_coconut.txt
TASKS=ibex/tasks_synthetic.txt OUT=$ROOT/results/interleaved     P=4 COCONUT=$COCONUT bash ibex/run_interleaved.sh
TASKS=ibex/tasks_coconut.txt   OUT=$ROOT/results/interleaved_all P=4 COCONUT=$COCONUT bash ibex/run_interleaved.sh
cat $ROOT/results/interleaved/rows/*.csv     > ibex/results/table5_interleaved.csv
cat $ROOT/results/interleaved_all/rows/*.csv > ibex/results/table5_interleaved_all.csv
python ibex/analysis/make_table_solver_interleaved.py     # results/table_solver.tex (Table 4)
```

All the variants of an instance (`af2` = S, `splus` = S⁺, `sh` = SH, `exact` = E) run one after the
other in the same session, four runs in parallel; the runs are resumable. The 110 COCONUT instances
of Table 4 are those in which S encloses some product more widely than E, counted by patch 0021 in
a run of 5 s with E:

```
for f in $COCONUT/*.bch; do
  $ROOT/ibexopt-raak/build/ibexopt-raak --lr=art --arith=exact --knomial_k=0 -t 5 $f | grep EXACT-MULT \
    > ibex/results/sdiag/$(basename $f .bch).txt
done
```

| File | Content |
|---|---|
| `results/table5_interleaved.csv` | Runs on synq, synt and synp |
| `results/table5_interleaved_all.csv` | Runs on the 298 COCONUT instances (60 s) |
| `results/sdiag/` | Counters of patch 0021 (5 s with E) used to select the 110 COCONUT instances of Table 4 |
| `analysis/make_table_solver_interleaved.py` | `results/table_solver.tex` (Table 4) |
| `analysis/summary_all_by_set.py`, `analysis/summary_all_subsets.py` | COCONUT results by subset |
| `analysis/make_table_cost.py` | `results/table_cost.tex` (Table 3) from `results/bench_mult_0019.csv` |
| `results/bench_mult*.csv`, `results/containment_*` | Cost per product and containment tests |
