# Formas afines de grado 2 en Ibex (`LinearizerQuad2`)

Integración en Ibex de la propuesta de `../../notas/raak-proto-cpp/`: formas afines de grado 2 con un encierro
final estructurado, usadas como relajación lineal (`art`) en `ibexopt-raak` e `ibexsolve-raak`.

| archivo | contenido |
|---|---|
| [`ibex_Quad2.h`](ibex_Quad2.h) | aritmética `c + aᵀz + zᵀQz` (z = variables y símbolos de error), cinta desde `ibex::Function`, encierros finales |
| [`ibex_LinearizerQuad2.h`](ibex_LinearizerQuad2.h) | linearizador, con la misma estructura que `LinearizerAffineT::linearize` |
| [`prodk.h`](prodk.h) | bloques producto (copia de `../../notas/raak-prodk-cpp/prodk.h`) |
| [`test_q2_validez.cpp`](test_q2_validez.cpp) | prueba de validez con soluciones conocidas |
| [`barrido_opt.py`](barrido_opt.py), [`barrido_solve.py`](barrido_solve.py), [`resumir_barrido.py`](resumir_barrido.py), [`resumir_csp.py`](resumir_csp.py) | barridos y resúmenes |

`../build.sh` compila todo. Las opciones nuevas, en los dos ejecutables, son `--arith=q2|q2kn|q2pc|q2full`:

| `--arith` | encierro final de la parte cuadrática de cada restricción | análogo a |
|---|---|---|
| `q2` | barato: cada símbolo exacto con su parte lineal y su diagonal, más el resto término a término | RAA k = 0, sobre la expresión completa |
| `q2kn` | k-nomial k = 6 (`ibex_KNomial.h`) | RAA-K, sobre la expresión completa |
| `q2pc` | bloques producto, pares + soporte completo (m ≤ 12) | SH, aplicado a productos de dos sumas extraídos de la Q completa |
| `q2full` | bloques producto completos (agrega triples) | lo mismo; demasiado lento en el solver |

## Cómo funciona

1. **Cinta:** al construirse, el linearizador traduce cada restricción a una cinta, recorriendo el DAG de la
   función.
   - Operaciones soportadas: símbolos (también índices de vectores), constantes (también intervalos), +, −, ×, /,
     `^n`, `sqr`, exp, log, sqrt, sin, cos y abs.
   - Las demás funciones unarias (tan, atan, cosh, …) se evalúan por intervalos.
   - Lo que no se puede traducir (vectores, max, min, apply) deja esa restricción sin fila.
2. **Evaluación en cada caja:** se calcula la forma de grado 2 y, como hace `AffineEval`, un intervalo por nodo:
   la operación por intervalos intersectada con el rango de la forma. Ese intervalo es el dominio de las
   linealizaciones y el `ev` del linearizador. Las linealizaciones son de Chebyshev, igual que el modo por
   omisión de ibex-affine.
3. **Fila del LP:** la parte lineal en las variables da la fila. El resto (parte cuadrática y símbolos de error)
   se acota con el encierro final y hace de error. Las comprobaciones con `ev` son las mismas de
   `LinearizerAffine2`.

## Rigor (pendiente)

La aritmética es en dobles. Se agrega un margen de seguridad de 10⁻¹³ relativo a la **magnitud acumulada** del
cálculo (la suma del tamaño de los términos de todos los nodos), en la intersección de cada nodo y en las filas.
No es una cota rigurosa del redondeo. Lo riguroso sería el esquema de `knomial_quad_bound`: el residuo en
intervalos.

**Validación** ([`test_q2_validez.cpp`](test_q2_validez.cpp)). Toma las soluciones que encuentra af2 y, en cajas
aleatorias que las contienen, comprueba que el intervalo `ev` y la relajación contengan el valor de cada
restricción:
- ponts-geo, Pramanik, Directkin, Neveu1 y Dietmaier: 0 violaciones en más de 2,2 millones de comprobaciones;
- bearing, sobre el sistema extendido del optimizador: 0 violaciones.

**Dos errores corregidos en el camino:**
1. Un margen relativo al rango del resultado perdía las 128 soluciones de ponts-geo, por cancelación cerca de 0.
2. Un margen de 10⁻⁹ era tan grande que impedía alcanzar la precisión pedida en f*.

**Coherencia en los barridos:** ninguna cota inferior supera a la mejor cota superior, en las instancias que
resuelven todas las variantes coincide f*, y en NCSP coinciden las soluciones.

## Resultados

### Optimización

Datos: [`opt_q2_t30.csv`](../results/opt_q2_t30.csv) y [`opt_q2_t30_resumen.txt`](../results/opt_q2_t30_resumen.txt).
Corridas de `ibexopt-raak --lr=art`, 30 s, 4 en paralelo, todas las variantes en la misma máquina. Las razones son
medias geométricas sobre las instancias que resuelven todas las variantes, relativas a af2.

| conjunto | variante | resueltas | celdas / af2 | CPU / af2 |
|---|---|---|---|---|
| Ibex medium+hard y COCONUT (41; 25 resueltas por todas) | af2 | 28 | 1 | 1 |
| | sh | 29 | 0,977 | 0,956 |
| | q2 | 27 | 0,802 | 0,801 |
| | **q2kn** | **28** | **0,776** | **0,787** |
| | q2pc | 27 | 0,765 | 0,903 |
| synq (78; 27 resueltas por todas) | af2 | 41 | 1 | 1 |
| | sh | 48 | 0,208 | 0,228 |
| | q2 | 78 | 0,034 | 0,043 |
| | **q2kn** | **78** | **0,023** | **0,035** |
| | q2pc | 77 | 0,022 | 0,040 |
| generador fqf, n = 10 y 20 (54; 23 resueltas por todas) | af2 | 33 | 1 | 1 |
| | sh | 33 | 0,998 | 0,949 |
| | q2 | 33 | 0,747 | 0,751 |
| | q2kn | 34 | 0,541 | 0,551 |
| | **q2pc** | **34** | **0,469** | **0,498** |

Casos que vale la pena destacar en Ibex medium+hard (celdas, af2 → q2kn): ex6_2_8 17 886 → 5830, ex6_2_12
3968 → 2228, ex7_3_4bis 1084 → 636, ex14_2_7 5432 → 4478. La pérdida es bearing: af2 la resuelve en 15 s y q2kn
necesita 56 s.

### NCSP

Datos: [`csp17_q2_t120.csv`](../results/csp17_q2_t120.csv) y [`csp17_q2_t120_resumen.txt`](../results/csp17_q2_t120_resumen.txt).
`ibexsolve-raak --ctc=default-xn+affine`, 120 s.

| instancia | af2 | sh | q2 | q2kn | q2pc |
|---|---|---|---|---|---|
| Pramanik | 5335 celdas, 1,7 s | 2609, 1,0 s | 1721, 1,2 s | **1087**, 3,0 s | 1215, 3,4 s |
| Directkin | 1993, 6,8 s | 1161, 4,0 s | 1319, 4,4 s | **1051**, 6,1 s | 1157, 5,7 s |
| Neveu1, ponts-geo | sin cambios (las mismas soluciones) | | | | |
| Dietmaier (120 s) | 30 soluciones | 46 | 30 | 42 | 34 |
| Trigexp2-\* | casi sin cambios; q2pc entre 20 y 40 % más lenta | | | | |

En las 4 instancias que resuelven todas, la media de celdas relativa a af2 es: sh 0,72, q2 0,68, q2kn 0,57 y
q2pc 0,60. En tiempo, q2kn y q2pc son entre 1,1 y 1,2 veces af2 (sh 0,78).


## Resultados con SoPlex 4.0.2 (resultados de referencia)

Compilación con `LP_LIB=soplex` (SoPlex 4.0.2, incluido en ibex-lib) y el parche 0018, que porta las
optimizaciones de SoPlex de `ibex-team/ibex-lib` develop del 4-11-2026:
- reloj apagado (`2db5d751`);
- Neumaier-Shcherbina sin matrices densas (`b86f8667`);
- cota de los LPs detenidos antes del óptimo (`9ea5ef6c`).

La evaluación de grado 2 es la optimizada: sin `std::map`, con memoria reutilizada, sin la evaluación por
intervalos redundante y con k = 3 en q2kn. Las cinco variantes corrieron en la misma compilación y la misma
máquina.

Datos:
- [`opt_q2_soplex_t30.csv`](../results/opt_q2_soplex_t30.csv) y [resumen](../results/opt_q2_soplex_t30_resumen.txt)
- [`csp17_q2_soplex_t120.csv`](../results/csp17_q2_soplex_t120.csv) y [resumen](../results/csp17_q2_soplex_t120_resumen.txt)
- costo por componente: [`soplex_costo_por_componente.txt`](../results/soplex_costo_por_componente.txt)

**Las filas son lineales y su costo de LP es el mismo en todas las aritméticas.** Las llamadas al LP por celda
(1,6–2,5), las filas por llamada y el costo del LP por celda son prácticamente iguales con af2, sh y q2\*. El
sobrecosto de LP que se veía con Clp no aparece con SoPlex.

### Optimización (30 s; medias geométricas relativas a af2 en las instancias que resuelven todas)

| conjunto | variante | resueltas | celdas / af2 | CPU / af2 |
|---|---|---|---|---|
| Ibex medium+hard y COCONUT (41; 25 resueltas por todas) | af2 | 29 | 1 | 1 |
| | sh | 29 | 0,919 | 0,966 |
| | **q2** | **30** | 0,732 | **0,691** |
| | **q2kn** | **30** | **0,716** | 0,728 |
| | q2pc | 29 | 0,870 | 1,128 |
| synq (78; 27 resueltas por todas) | af2 | 42 | 1 | 1 |
| | sh | 52 | 0,206 | 0,253 |
| | q2 | 78 | 0,028 | 0,038 |
| | **q2kn** | **78** | 0,019 | **0,028** |
| | q2pc | 77 | 0,018 | 0,033 |
| generador fqf (54; 22 resueltas por todas) | af2 | 33 | 1 | 1 |
| | sh | 33 | 1,001 | 0,983 |
| | q2 | 34 | 0,749 | 0,693 |
| | q2kn | 36 | 0,533 | 0,486 |
| | **q2pc** | **36** | **0,465** | **0,437** |

### NCSP (120 s)

| instancia | af2 | sh | q2 | q2kn | q2pc |
|---|---|---|---|---|---|
| Pramanik | 5335 celdas, 1,2 s | 2609, 0,9 s | 1721, **0,6 s** | **1087**, 0,8 s | 1215, 2,8 s |
| Directkin | 1991, 4,2 s | 1153, **2,7 s** | 1323, 2,8 s | **1051**, 2,9 s | 1155, 4,6 s |
| Dietmaier | sin terminar en 120 s, 50 soluciones | 54 soluciones en 84 s | sin terminar, 49 | 54 en 89 s | 54 en 111 s |
| Trigexp2-\*, Neveu1, ponts-geo | prácticamente iguales; q2pc más lenta | | | | |

En las 4 instancias que resuelven todas, relativo a af2:

| variante | celdas | tiempo |
|---|---|---|
| sh | 0,72 | 0,88 |
| q2 | 0,68 | **0,74** |
| q2kn | **0,57** | 0,83 |
| q2pc | 0,60 | 1,32 |

### Lectura con SoPlex

1. **q2 y q2kn son las mejores variantes en general.** En Ibex medium+hard usan 0,72–0,73 veces las celdas de af2
   y 0,69–0,73 veces el tiempo, y resuelven 30 instancias frente a 29. En synq resuelven las 78 (af2, 42). En NCSP
   bajan celdas y tiempo frente a af2 y a SH.
2. **q2pc es la mejor en formas cuadráticas densas** (el generador: 0,47 de las celdas y 0,44 del tiempo), pero
   cuesta más por celda en el resto, porque recalcula los bloques producto en cada evaluación.
3. **SH por producto,** en cambio, casi no gana en Ibex medium+hard (0,92 de las celdas) y nada en el generador.

## Lectura (barrido anterior, con Clp)

1. **En optimización, q2kn es la mejor variante:** con unas 0,78 veces las celdas y el tiempo de af2 en Ibex
   medium+hard, unas 0,55 veces en el generador y unas 0,03 veces en synq, resolviendo las 78 de synq frente a 41
   con af2. Le gana a SH por producto en los tres conjuntos.
2. **q2pc reduce más las celdas** (en el generador, 0,47), pero cuesta más por celda. Conviene cuando la forma es
   chica o densa.
3. **En NCSP baja las celdas a la mitad o menos donde hay productos de factores compuestos**, pero el costo por
   celda todavía no compensa frente a SH en tiempo. Ahí falta optimizar la evaluación de grado 2.
4. **synq exagera la ventaja** por cómo se genera, con productos que se cancelan a propósito. Los números
   representativos son los de Ibex medium+hard y el generador.

## Pendiente

- Rigor en punto flotante: el residuo y las filas con aritmética de intervalos.
- Bajar el costo por celda: vectores dispersos, menos `std::map` y reutilizar memoria entre restricciones.
- Elegir el encierro final según el tamaño de la forma (k-nomial o bloques producto).
- Barridos más largos y con más instancias, incluidas COCONUT completa y QPLIB.

## Experimento: aproximación cuadrática de las funciones (resultado negativo)

**Qué se probó.** En q2, cada función se linealiza con Chebyshev sobre el intervalo del nodo:
`f(X) ≈ α·X + β ± δ`. Eso conserva la estructura cuadrática del argumento (escalada por α), pero la curvatura
de la propia función va al error.

Se probó una aproximación cuadrática `f(X) ≈ f(m) + α(X−m) + γ(X−m)² ± δ`, interpolando en los 3 nodos de
Chebyshev del dominio:
- `(X−m)²` se calcula como un producto q2 (grado 2 conservado, resto de grado 3–4 a error);
- solo se usa cuando el error total es menor que el de la linealización.

Se activa con la variable de entorno `Q2_UNARY` (por omisión, `lin`):

| `Q2_UNARY` | qué hace |
|---|---|
| `lin` | linealización de Chebyshev (por omisión) |
| `quad` | aproximación cuadrática con cota rigurosa `max|f'''|·r³/24` (`f'''` acotada con intervalos) |
| `quadmuestreo` | aproximación cuadrática con el error medido por muestreo denso; no es rigurosa, sirve como techo de la idea |

**Medición** ([`test_q2_ancho.cpp`](test_q2_ancho.cpp)): ancho del encierro de cada restricción en cajas
aleatorias de escala 1, 10⁻¹, 10⁻² y 10⁻³, comparado con `lin` en las mismas cajas.

| instancias | resultado |
|---|---|
| ex6_2_6, ex6_2_8, chem (`ln`), Pramanik | **nunca se usa**, ni siquiera con el error ideal: `ln` es singular en 0 y el dominio llega a 10⁻⁶, así que ninguna parábola le gana a la recta de Chebyshev |
| ex14_2_7, bearing, Trigexp2-9 | se usa miles de veces, pero el ancho no cambia (1,000) |
| Directkin | entre 0,98 y 1,02 según la escala |

**Conclusión:** no vale la pena. El error de `ln` en ex6_2_\* está dominado por la singularidad cerca de 0 (por
ejemplo, en los términos `x·ln(x)`), y eso no lo arregla una aproximación cuadrática. Se queda la
linealización de Chebyshev.
