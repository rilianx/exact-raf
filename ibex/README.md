# Multiplicación exacta en Ibex (parche 0016)

Implementación en `ibex-affine` (modelo `AF_fAF2`) de la multiplicación exacta del paper
*Exact multiplication of reduced affine forms* (`paper/numa/`) y de la de Skalna y Hladík (2017, Ec. 10),
sobre la cadena de parches RAA-K de JCGeraldo (`JCGeraldo/ibex-lib`, rama `raak-experimentos`,
carpeta `raak/`). Especificación: `paper/numa/IBEX_EXACT_SPEC.md`.

| Archivo | Contenido |
|---|---|
| `patches/0016-ibex-affine-multiplicacion-exacta.patch` | Sobre `ibex-affine` `c219009c` + `raak-ibex-affine-consolidado.patch` (0001–0015). Agrega `ibex_ExactMult.h`, los modos `AFFINE_SH` y `AFFINE_EXACT` y la rama correspondiente en `operator*=`; `Asqr` usa en ambos modos el cuadrado exacto de `AFFINE_RAA`. |
| `patches/0017-ibex-affine-multiplicacion-rapida.patch` | Sobre 0016. `ibex_ExactMultFast.h`: recorrido en `double` con cota a priori del redondeo, predicado angular filtrado, sin asignaciones y las cuatro esquinas en una pasada (ver `fast/`). Con `IBEX_RAF_SQR=af2`, `sh`/`exact` usan el cuadrado de fAF2 (ablación). |
| `cli/` | `ibexopt-raak` e `ibexsolve-raak` de JCGeraldo con `--arith=sh` y `--arith=exact` y la línea `EXACT-MULT stats`. |
| `tests/affine_containment_check.cpp` | Prueba de contención de JCGeraldo (sin tolerancia, en intervalos) con las variantes `sh_k0` y `exact_k0` y la comparación de anchos. |
| `tests/bench_mult.cpp` | Tiempo por producto según `n`, con y sin error en los factores (resultados en `results/bench_mult.csv`). |
| `build.sh` | Clona, parcha y compila todo (Ibex con `LP_LIB=clp` o `soplex`). |
| `run_experiments.sh` | Barridos de `ibexopt` (reales, synq, synl, synt) y de `ibexsolve` (17 instancias con material) con `af2`, `raa`, `sh`, `exact`. |
| `results/` | Resultados (CSV) de las pruebas y los barridos. |

## Rigor

`raf_mult_enclosure` devuelve un intervalo que contiene el rango de la parte no lineal
`N = x0 ey η + y0 ex ξ + (aᵀε + ex ξ)(bᵀε + ey η)`:
orden angular con predicado exacto (`__float128`), vértices y valores en las aristas con `ibex::Interval`,
y el centro y el radio se suman con los mismos `twoSum` y factores `(1+2·AF_EM)` de fAF2.

## Validación

- `affine_containment_check 20000 11` (9 familias × 20 000 funciones × 7 variantes): **0 violaciones seguras**;
  con el error achicado al 90 % detecta 1696 (`sh`) y 1719 (`exact`) violaciones, así que la prueba es sensible.
- En productos aislados (200 000 casos), el encierro de `exact` está siempre contenido en el de `sh`.
- `ex5_4_3` (`ibexopt`, `art`): mismo óptimo con las cuatro aritméticas; 42 celdas con `af2`/`raa`, 40 con `sh`/`exact`.

## Parche 0017 (versión rápida)

- `affine_containment_check 20000 11`: **0 violaciones seguras** (`results/containment_20000_0017_stderr.txt`).
- Verificación contra el rango exacto en aritmética racional: 86 000 casos, 0 violaciones (`fast/`).
- Costo por producto (`results/bench_mult_0017.csv`): `exact` cuesta 2.4 a 5.4 veces `af2` (antes 7 a 45)
  y es de 1.5 a 18 veces más barato que `raa`.
