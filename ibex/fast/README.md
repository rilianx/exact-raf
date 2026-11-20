# Versión rápida de `raf_mult_enclosure`

`ibex_ExactMult_fast.h`: mismo algoritmo que `ibex_ExactMult.h` (parche 0016), pero:

- el recorrido del borde se hace en `double` y se agrega al final una cota a priori del error de redondeo
  `M = (8m + 64)·eps·(U·V + |x0|ey + |y0|ex)`, con `eps = 2^-52` (vale para cualquier modo de redondeo),
  y se redondea hacia afuera con `nextafter`;
- el predicado angular usa un filtro en `double` y recurre a `__float128` solo si el signo es dudoso;
- los buffers son `thread_local` (sin asignaciones por producto);
- `exact` recorre las cuatro esquinas en una sola pasada.

Verificación (`test_fast 0 <semilla> <casos>` + `verify.py`, rango exacto con `fractions.Fraction`):
86 000 casos, incluidos generadores casi paralelos y escalas de 1e-6 a 1e6: **0 violaciones**;
holgura relativa ≤ 1.1e-12. Tiempos (`test_fast 1`): `results/bench_fast.csv`, de 5 a 10 veces más rápida
que la versión con `ibex::Interval`.

Compilar: `g++ -std=gnu++11 -O2 -msse3 -frounding-math -I. test_fast.cpp $(pkg-config --cflags ibex) -o test_fast $(pkg-config --libs ibex) -lquadmath`
(con `ibex_ExactMult_interval.h` = copia de `ibex-affine/src/arithmetic/ibex_ExactMult.h`).
