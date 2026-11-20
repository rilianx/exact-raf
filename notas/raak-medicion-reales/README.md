# Medición en la raíz de instancias reales: sumas de productos

Continuación de [`../raak-medicion-synq/`](../raak-medicion-synq). La pregunta es si vale la pena acotar en
conjunto las sumas de productos que aparecen dentro de una expresión, ya sea expandiéndolas o manteniendo los
términos bilineales en la forma afín, en lugar de encerrar cada producto con SH y sumar. Es una medición
estática en la caja raíz declarada en cada `.bch`: no se corre el solver.

**Instancias (41):** las 38 medium+hard de `benchs/optim` de Ibex y las 3 de COCONUT de
`raak/instancias/reales` (`3pk-borne`, `mistake`, `s365mod`), todas de `JCGeraldo/ibex-lib`, rama
`raak-experimentos`.

## Método

1. En cada objetivo y restricción se buscan los **nodos suma**: combinaciones lineales maximales, con +, −, y
   producto o división por constante.
2. Se miden los que tienen **2 o más términos cuadráticos**, es decir, `A*B` con A y B no constantes, o `A^2`.
   A y B pueden ser subexpresiones compuestas.
3. Cada factor se evalúa en aritmética afín: `A = a₀ + aᵀε + e_A ξ_A`, con min-range para exp, ln, sqrt y 1/x.
   El error `ξ_A` es un símbolo propio de la subexpresión A, y es el mismo si A se repite.
4. El residuo que se compara es `N = Σₖ cₖ (uₖᵀζ)(vₖᵀζ)`, con `ζ = (ε, ξ)`. El centro y los términos lineales
   son iguales para todos los métodos.

Métodos, con las mismas columnas que en synq:

| columna | qué es |
|---|---|
| `S` | suma de la fórmula de fAF2 por producto |
| `SH` | suma del rango exacto de cada producto (**línea base**) |
| `Qk0`, `Qk2`, `Qk6` | Q agregada: término a término, y k-nomial SDSOS con k = 2 y k = 6 (`ibex_KNomial.h`) |
| `EIG` | cota espectral de la Q agregada |
| `BEST` | la menor entre SH, Qk6 y EIG |
| `EXACT` | rango exacto si hay a lo más 12 símbolos. Con más, una cota interior obtenida por búsqueda local, marcada (≥). |

`like` y `like-1` dividen por expresiones cuyo intervalo contiene 0 en la raíz, así que no tienen sumas medibles.

## Resultados

Se encontraron **93 nodos suma con 2 o más términos cuadráticos, en 36 de las 41 instancias**. En 37 de ellos
los términos comparten algún símbolo.

Ancho relativo a SH, primero media geométrica sobre los nodos de cada instancia y después sobre las instancias
(cada instancia pesa 1):

| | S | **SH** | Qk0 | Qk2 | Qk6 | EIG | BEST | EXACT |
|---|---|---|---|---|---|---|---|---|
| 36 instancias | 1,05 | 1 | 0,99 | 0,97 | 0,95 | 1,74 | 0,93 | 0,92 |
| 25 con símbolos compartidos | 1,06 | 1 | 0,98 | 0,95 | 0,93 | 1,77 | 0,90 | 0,89 |

- **SH ya es exacto en 87 de los 93 nodos.** El margen que queda es, en media geométrica, de un 8 % del ancho, y
  se concentra en pocas instancias:

  | instancia | Qk0 | Qk6 | EIG | EXACT | de dónde sale la ganancia |
  |---|---|---|---|---|---|
  | ex6_2_6 | 0,35 | 0,35 | 0,77 | 0,35 | productos idénticos con coeficientes que se cancelan (`x2*ln(x2)` aparece 4 veces y sus coeficientes suman ≈ 0) y el mismo `ln(…)` en varios productos |
  | ex6_2_8 | 0,40 | 0,40 | 0,97 | 0,40 | ídem |
  | ex2_1_9, ex2_1_9fact | 1 | 1 | **0,76** | 0,73 | suma de 22 bilineales `xᵢxⱼ` con diagonal cero: el k-nomial SDSOS no forma bloques, la cota espectral sí aprovecha |
  | ex8_5_1bis, ex8_5_2_1 | 1 | 0,95 | 1,07 | 0,88 | k-nomial y espectral; la menor de las dos da 0,90 |

- **En ninguna de estas instancias el k-nomial es lo decisivo.**
  - En ex6_2_6 y ex6_2_8 basta con juntar los términos iguales (Qk0), lo que hace también una simplificación
    simbólica.
  - En ex2_1_9 el k-nomial no sirve, porque necesita diagonal.
  - Solo en ex8_5_* aporta algo (5 %).
- **El k-nomial agregado es peor que SH en 15 de los 93 nodos.** Ocurre en ex6_1_3bis, ex8_4_4bis y 3pk-borne,
  con hasta 1,09× SH en 3pk-borne y 1,20× en ex6_1_3bis. Por eso, si se implementa, hay que quedarse con la menor
  de las cotas, y SH tiene que estar entre ellas.
- **La cota espectral sola es mala en instancias reales** (1,74× SH en media): las cajas no son cubos y Q no está
  normalizada. Solo sirve como una cota más dentro del mínimo.

## Comparación con synq

| | synq (78) | reales (36) |
|---|---|---|
| SH / exacto | 6,0 | 1,08 |
| mejor agregado / SH | 0,21 | 0,93 |

La gran ganancia de synq viene de su generador, que factoriza `Q = AᵀM` con A aleatoria y crea a propósito
mucha cancelación entre productos. En los benchmarks reales esa estructura casi no aparece. Las sumas de
productos son en su mayoría de dos tipos:
- sumas de términos con poca interacción, donde SH ya es exacto;
- términos repetidos, que se resuelven juntando términos iguales y no requieren el k-nomial.

## Salvedades

- **Solo la caja raíz.** En cajas internas los radios cambian, pero Q en ε se escala por los radios, así que la
  estructura (qué interactúa con qué) es la misma.
- **Errores compartidos.** La ganancia de ex6_2_6 y ex6_2_8 supone que el error de `ln(…)` es un símbolo compartido
  entre productos. fAF2 junta todos los errores en un solo término por forma, así que hoy esa correlación se
  pierde igual. Un acumulador cuadrático en fAF2 tendría que mantener esos símbolos.
- **El exacto es una cota interior en 16 nodos** (más de 12 símbolos). En todos ellos coincide con SH o con BEST,
  así que no cambia las conclusiones.

## Reproducir

```bash
# knomial_driver: ver ../raak-medicion-synq/README.md
B=jcgeraldo-ibex-lib
python3 medir_reales.py ./knomial_driver nodos.csv $B/benchs/optim/medium/*.bch $B/benchs/optim/hard/*.bch \
    $B/raak/instancias/reales/{3pk-borne,mistake,s365mod}.bch
python3 resumir.py nodos.csv > resumen.txt
```

Datos por nodo en [`nodos.csv`](nodos.csv) y tabla por instancia en [`resumen.txt`](resumen.txt).
