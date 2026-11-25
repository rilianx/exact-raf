# Fast version of `raf_mult_enclosure` (patch 0017)

`ibex_ExactMult_fast.h` implements the same algorithm as `ibex_ExactMult.h` (patch 0016), but:

- the walk along the boundary is done in `double`, and an a priori bound of the rounding errors,
  `M = (8m + 64)·eps·(U·V + |x0|ey + |y0|ex)` with `eps = 2^-52` (valid for any rounding mode),
  is added at the end, followed by outward rounding with `nextafter`;
- the angular predicate uses a filter in `double` and falls back to `__float128` only when the sign is doubtful;
- the buffers are `thread_local` (no allocation per product);
- E walks its corners in a single pass.

Verification (`test_fast 0 <seed> <cases>` + `verify.py`, exact range with `fractions.Fraction`):
86 000 cases, including nearly parallel generators and scales from 1e-6 to 1e6: **0 violations**;
relative slack ≤ 1.1e-12. Times (`test_fast 1`): `../results/bench_fast.csv`, 5 to 10 times faster
than the version with `ibex::Interval`.

Build: `g++ -std=gnu++11 -O2 -msse3 -frounding-math -I. test_fast.cpp $(pkg-config --cflags ibex) -o test_fast $(pkg-config --libs ibex) -lquadmath`
(with `ibex_ExactMult_interval.h` = a copy of `ibex-affine/src/arithmetic/ibex_ExactMult.h`).
