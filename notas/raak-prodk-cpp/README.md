# Bloques producto en C++

Port a C++ de la reducción ávida por bloques producto (PRODK), que sirve de encierro final de las formas
cuadráticas que deja M3 centrado (`../raak-errores-simbolos/`). Se mide su precisión y su costo frente al
k-nomial en C++ y frente a Ibex por omisión.

| archivo | contenido |
|---|---|
| [`prodk.h`](prodk.h) | la reducción, solo cabecera y en dobles (ver Rigor) |
| [`prodk_driver.cpp`](prodk_driver.cpp) | una forma por línea: `m h11 … hmm`. Salida: `L U bloques tiempo_us ancho_tt`. Argumentos: `[REPS] [LAZY] [MAX_TRIPLES_M] [FULL_MAX]` |
| [`volcar_formas.py`](volcar_formas.py) | vuelca las 5037 formas homogeneizadas H de M3 centrado (mismas cajas que las mediciones; `formas_meta.csv`) |
| [`comparar.py`](comparar.py), [`validar.py`](validar.py) | fidelidad frente a Python y validez frente al rango exacto |
| [`tabla_cpp.py`](tabla_cpp.py) | precisión final de las expresiones y tiempos, en [`tabla_cpp.txt`](tabla_cpp.txt) |

## Qué hace

Para acotar `zᵀHz` en [−1, 1]ᵐ:

**1. Soportes.** Se recorren:
- pares de símbolos con entrada no nula;
- triples con 2 o más aristas;
- el soporte completo, si m ≤ 12.

**2. Candidatos por soporte.** Con los autovalores de la submatriz, calculados por Jacobi:
- el par espectral `(a+b)(a−b)`;
- el autovector dominante como cuadrado;
- las filas de factor común `zᵢ·(…)`.

Cada uno se prueba escalado por 1, ¾, ½ y ¼.

**3. Ganancia.** Para cada candidato se calcula la reducción del encierro término a término de la submatriz,
menos el rango SH del bloque. El rango SH usa orden angular exacto por producto cruz.

**4. Bucle ávido.** Se extrae el mejor bloque, se actualiza el residuo y se repite.

**5. Encierro.** La suma de los rangos SH de los bloques más el término a término del residuo.

## Optimizaciones

Todas dan exactamente el mismo resultado que la versión anterior (razón de anchos 1,000000), salvo la primera:

| optimización | efecto en el tiempo |
|---|---|
| **Búsqueda perezosa:** los soportes afectados se recalculan solo cuando llegan al tope del montículo | 2–8 veces más rápida; ancho casi igual (razón 0,9999) |
| **Sin reservas de memoria en el bucle:** arreglos fijos de hasta 12 símbolos | 2,5–3 veces |
| **SH una sola vez por candidato:** `SH(t·u, v) = t·SH(u, v)` | ~1,5 veces |

## Validez y fidelidad

- **Validez:** en las 4361 formas con m ≤ 12, ningún encierro queda por debajo del rango exacto (enumeración de
  caras KKT). La mediana del ancho frente al exacto es **1,000** (p90 1,12): en la mayoría de esas formas, los
  bloques producto alcanzan el rango exacto.
- **Fidelidad frente a Python** (`prodk_fast.py`): los anchos difieren en el 31 % de las formas. La media
  geométrica es 0,989, es decir, C++ es algo más angosto y siempre válido. La diferencia viene de los
  autovectores: Jacobi y LAPACK no los devuelven iguales, eso cambia qué candidato gana y la búsqueda ávida sigue
  otro camino.

## Resultados

**Precisión final de las expresiones** (7532 combinaciones de expresión y caja, 58 instancias). El encierro final
de M3 centrado es la menor entre el barato y la variante indicada. Referencia: M00, Ibex por omisión (fAF2 con la
fórmula S).

| encierro final de M3 centrado | ancho/M00, optim. | ancho/M00, NCSP | instancias que mejoran > 2 % | ancho/INNER, optim. | ancho/INNER, NCSP | tiempo C++ del encierro (mediana / p90, µs) |
|---|---|---|---|---|---|---|
| (fAF2 + SH, M0, referencia) | 0,987 | 0,947 | 12 | 1,74 | 1,50 | — |
| encierro barato | 0,950 | 0,955 | 17 | 1,67 | 1,52 | < 0,1 |
| k-nomial k = 6 | 0,946 | 0,922 | 17 | 1,66 | 1,46 | 0,6 / 2,4 |
| bloques producto: solo pares | 0,947 | 0,923 | 17 | 1,67 | 1,47 | 3,8 / 19 |
| bloques producto: pares + completo (m ≤ 12) | 0,934 | 0,914 | 27 | 1,64 | 1,45 | 12 / 54 |
| bloques producto: pares + triples | 0,926 | 0,905 | 31 | 1,63 | 1,44 | 16 / 238 |
| **bloques producto: pares + triples + completo** | **0,925** | **0,905** | **32** | 1,63 | 1,44 | 24 / 327 |
| SDP (Python, referencia) | 0,910 | 0,903 | 30 | 1,60 | 1,43 | — |

**Tiempo por forma, según su tamaño** (mediana / p90, µs):

| m | formas | k-nomial k = 6 | pares + completo | pares + triples | pares + triples + completo |
|---|---|---|---|---|---|
| 2–5 | 2691 | 0,4 / 1 | 4,8 / 14 | 7,1 / 17 | 9,3 / 25 |
| 6–10 | 825 | 0,7 / 1 | 24 / 44 | 39 / 84 | 54 / 128 |
| 11–20 | 1303 | 1,7 / 3 | 37 / 146 | 156 / 403 | 215 / 425 |
| 21–40 | 118 | 5,8 / 359 | 55 / 286 | 725 / 14 154 | 723 / 14 122 |
| 41–200 | 91 | 16 / 85 | 54 / 209 | 53 / 187 | 53 / 184 |

Con m > 40 no se usan triples; con m > 12, tampoco el soporte completo.

**Escala de referencia.** En C++, un producto de fAF2 cuesta entre 0,08 y 0,6 µs, y uno SH entre 0,26 y 1,5 µs,
para n de 2 a 20 (`../../ibex/results/bench_mult*.csv`). Una expresión de algunas decenas de operaciones se
evalúa en pocos µs.

## Lectura

1. **En C++, los bloques producto completos casi alcanzan a la SDP:** 0,925 frente a 0,910 en optimización y
   0,905 frente a 0,903 en NCSP, con 32 instancias que mejoran.
2. **Pero cuestan bastante más que el k-nomial:** 24 µs de mediana frente a 0,6 µs, y colas de milisegundos con
   m entre 21 y 40. Eso equivale a varias evaluaciones completas de la expresión.
3. **Variante intermedia: pares + soporte completo.** Con 12 µs de mediana y sin colas largas (p90 54 µs) da
   0,934 y 0,914, y 27 instancias que mejoran. Es la mejor relación entre precisión y costo de las medidas.
4. **Solo pares** no aporta sobre el k-nomial: la misma precisión, 6 veces más cara.
5. **Para Ibex, lo razonable sería:**
   - usar el k-nomial siempre, porque es casi gratis;
   - usar bloques producto (pares + completo) solo cuando la forma es chica o en nodos seleccionados;
   - de forma más general, decidir según el tamaño m de la forma.

## Rigor (pendiente para Ibex)

La versión es para medir: no usa redondeo dirigido. La cota vale para bloques cualesquiera, porque se acota el
residuo completo. Por eso en Ibex basta el mismo esquema de `knomial_quad_bound` (parches 0012 y 0013): calcular
los bloques en dobles y, al final, recalcular el residuo y los rangos SH con aritmética de intervalos. Para los
rangos SH de los bloques sirve la implementación rigurosa de `../../ibex/fast/`.
