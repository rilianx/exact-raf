# RAA-K: lo que cambia con Skalna–Hladík y hacia dónde reenfocar

Nota para la línea RAA-K (`JCGeraldo/ibex-lib`, rama `raak-experimentos`, carpeta `raak/`). Resume lo que
aprendimos al implementar en Ibex la multiplicación de Skalna y Hladík (SH) y la multiplicación exacta (E) del paper
*Exact multiplication of reduced affine forms* (`paper/numa/` de este repo), y propone un cambio de foco.

## En una frase

Para **un producto** de dos formas afines, la parte cuadrática `(aᵀε)(bᵀε)` tiene matriz de rango 2 y SH ya calcula
su rango **exacto** en O(n log n). Ninguna cota k-nomial puede ganarle ahí. El k-nomial debería competir donde no
hay un método exacto barato: **formas cuadráticas generales `εᵀQε` con Q de rango > 2**, es decir, funciones
cuadráticas completas (objetivos y restricciones) o sumas de varios productos, no la multiplicación afín operación
por operación.

## 1. Qué hace hoy RAA-K

En `operator*=` de `AF_fAF2` (ver la tabla de `ibex_KNomial.h`):

- **RAA, k = 0:** acota la parte cuadrática del producto término a término (forma simétrica), con un corner test
  para los cruces con los errores de los factores. Para el cuadrado usa Chebyshev, `[0, (Sx+e)²]`.
- **RAA, k ≥ 2:** descompone `εᵀQε` del producto en bloques SDSOS de hasta k variables
  (`knomial_quad_bound`), cada bloque con su rango exacto y el residuo término a término.

Todo es **por producto**: Q es la matriz de un solo producto, `Q = ½(abᵀ + baᵀ)`.

## 2. Qué cambia SH (y E)

- **SH** (Skalna y Hladík, *Numer. Algorithms* 76, 2017, Ec. 10) calcula el rango exacto de
  `R = (aᵀε + ex ξ)(bᵀε + ey η)` como el rango de `u·v` sobre un zonotopo plano. Lo hace en O(n) con selección
  de medianas, u O(n log n) ordenando por ángulo.
- **E** (nuestro paper) va un paso más allá: también incluye exactamente los términos `x0·ey·η + y0·ex·ξ`. Calcula
  el rango exacto de **toda** la parte no lineal del producto.
- **Consecuencia:** como Q tiene rango 2, `εᵀQε = (aᵀε)(bᵀε)` y su rango exacto es barato. Cualquier cota por
  bloques (k-nomial, SDSOS, término a término) es mayor o igual que el rango exacto. **Dentro de un producto, el
  k-nomial a lo más se acerca a SH, nunca lo supera, y cuesta más.**

### Evidencia (IbexOpt, relajación `art`, misma máquina)

| | celdas / af2 | CPU / af2 | fuente |
|---|---|---|---|
| synq, RAA k = 0 | 0,669 | 0,739 | 45 instancias resueltas por todas, 300 s |
| synq, SH | **0,259** | **0,304** | ídem |
| synq, E | 0,259 | 0,315 | ídem |
| reales, RAA k = 0 | 1,000 | 0,963 | 6 instancias |
| reales, SH | 0,916 | 0,907 | ídem |

- **k = 6 frente a k = 0 (estudio de semillas de RAA-K):** en synq con n = 8, k = 6 usa 0,37–0,38 de las celdas de
  RAA k = 0 (`raak/results/20-benchs-solver/README.md`, resultado 4). Combinado con RAA k = 0 / af2 ≈ 0,67, da
  **k = 6 / af2 ≈ 0,25**, prácticamente lo mismo que SH (0,26). Son sesiones y subconjuntos distintos, así que es
  solo una indicación, pero calza con la teoría: subir k converge a SH, que ya es exacto.
- **Costo por producto** (µs, C++, factores con error):

  | n | af2 | RAA k = 0 | SH | E |
  |---|---|---|---|---|
  | 2 | 0,08 | 0,39 | 0,26 | 0,26 |
  | 20 | 0,60 | 10,3 | 1,51 | 1,68 |
  | 100 | 2,1 | 204 | 10,6 | 11,4 |

  RAA es O(n²); SH y E son O(n log n).

Datos: `ibex/results/` de este repo (`opt_reales_t300.csv`, `opt_sinteticas_t300.csv`, `bench_mult*.csv`).

## 3. Hacia dónde reenfocar: `xᵀQx`, no la multiplicación

### Dónde aparece una Q de rango > 2

1. **Funciones cuadráticas completas.** Si un objetivo o una restricción es un polinomio de grado 2,
   `f(v) = c + bᵀv + vᵀQv`, sobre la caja queda `f = f0 + gᵀε + εᵀQ̃ε`, con Q̃ de rango hasta n. synq es
   exactamente esto: una suma de K productos de formas lineales da una Q de rango hasta 2K.
2. **Sumas de productos dentro de una expresión.** La aritmética afín encierra cada producto por separado
   (con SH o E) y suma los errores, así que pierde las cancelaciones entre productos. Acotar la suma de las partes
   cuadráticas en conjunto puede ser bastante más ajustado.

### Por qué ahí sí tiene sentido el k-nomial

El rango exacto de `εᵀQε` sobre `[-1,1]ⁿ` con Q indefinida es NP-difícil (es un QP no convexo en caja). Entre
"producto por producto" (barato, pierde correlaciones) y "exacto" (caro), el SDSOS por bloques de k variables es
una cota intermedia razonable, y el parámetro k regula costo y precisión.

### Propuesta concreta

**Relajación lineal de funciones cuadráticas por bloques k-nomiales.**

- **Entrada:** una restricción u objetivo cuadrático (detectado simbólicamente) y una caja.
- **Parte lineal:** `f0 + gᵀε`, la misma que da la aritmética afín, para comparar en igualdad de condiciones.
- **Residuo:** `εᵀQ̃ε` acotado por k-nomial SDSOS sobre la **Q completa** de la función, con el rigor que ya
  tienen los parches 0012–0013 (residuo en intervalos).
- **Integración en Ibex:** como un linearizador propio para restricciones cuadráticas (en la línea de
  `LinearizerAffine2`), no dentro de `operator*=`. Las restricciones no cuadráticas siguen con la aritmética afín
  usando SH o E.

Otra alternativa, más general y más cara de implementar, es una aritmética que guarde los términos de segundo orden
hasta el final de la expresión. Antecedentes a revisar: las formas cuadráticas confiables de Messine y Touhami
(*Reliable Computing*, 2006; verificar la referencia exacta) y los modelos de Taylor de orden 2.

### Líneas base obligatorias

Cualquier experimento nuevo debe compararse con:

1. **S (fAF2)**, la multiplicación por omisión de Ibex;
2. **SH o E producto por producto**, sumadas: es la línea base correcta, **no af2**;
3. **el rango exacto** de `εᵀQ̃ε` en la raíz (IbexOpt, o enumeración para n chico), para medir cuánto margen queda;
4. si es posible, **otra cota conocida para formas cuadráticas en caja**, como cotas espectrales o SDP.

## 4. Primer paso sugerido (barato, antes de implementar nada)

En la caja raíz de las 78 instancias synq, medir el ancho del residuo cuadrático del objetivo con:

- la suma de SH por producto;
- el k-nomial sobre la Q completa, con k = 2, 4, 6;
- el rango exacto.

Si el k-nomial sobre la Q completa no reduce claramente el ancho frente a la suma de SH (digamos menos de 10–20 %),
el enfoque no vale la pena en el solver y conviene saberlo antes de integrarlo.

## 5. Qué se puede reutilizar de este repo

- `ibex/patches/0016-…` y `0017-…`: modos `AFFINE_SH` y `AFFINE_EXACT` en fAF2.
  - `--arith=sh|exact` en `ibexopt-raak` e `ibexsolve-raak`.
  - `IBEX_RAF_SQR=af2` para usar el cuadrado de fAF2 en esos modos (ablación).
- `ibex/fast/`: implementación en double con cota a priori del redondeo y predicado angular exacto, verificada con
  aritmética racional. Es un patrón de rigor reutilizable para el k-nomial.
- `ibex/tests/affine_containment_check.cpp`: incluye las variantes `sh_k0` y `exact_k0`.
- `ibex/cli/cli_ibexsolve_raak.cpp`: `--ctc=default-xn+affine`, el contractor por omisión con la relajación afín
  en lugar de X-Newton.
- `ibex/analysis/`: scripts de resumen (media geométrica, resueltas, consistencia de f*) y la tabla del solver.

## 6. Lecciones de los experimentos

- **La línea base debe ser SH, no af2.** Contra af2 casi todo parece mejorar, porque fAF2 usa la fórmula S, que
  pierde hasta un factor 2 en la parte cuadrática.
- **Los cuadrados de expresiones compuestas pesan mucho.** En los benchmarks clásicos (medium+hard de Ibex),
  toda la mejora de SH frente a fAF2 viene del cuadrado de Chebyshev de expresiones como `(xi − xj)²`. Con el
  cuadrado de fAF2, SH no mejora nada (1,00). Hay que separar siempre los efectos de producto y de cuadrado.
- **S ya es exacto para productos y cuadrados de variables sueltas** (`xi·xj`, `xi²`). Las diferencias solo
  aparecen con factores compuestos, así que conviene elegir instancias que los tengan.
- **Hay mucho ruido por instancia.** Cambios del orden de 10⁻¹³ en los encierros cambiaron las celdas de algunas
  instancias en un factor mayor que 5 (por ejemplo `himmel16`). Solo hay que reportar medias geométricas sobre
  muchas instancias y desconfiar de los casos sueltos.
- **NCSP (17 instancias con material, `default-xn+affine`, 120 s):**
  - SH/E gana claro en Pramanik (2 609 celdas frente a 5 335 de af2 y 21 251 del contractor por omisión de Ibex),
    en Directkin (1 161 frente a 1 993) y en Dietmaier (48–49 de 54 soluciones validadas en 120 s, frente a 33);
  - en Trigexp2 la multiplicación no cambia nada.

  Datos: `ibex/results/csp17_*.csv`.
