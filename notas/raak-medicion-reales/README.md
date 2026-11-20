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
2. Se miden los que tienen **2 o más términos cuadráticos**. Un término es cuadrático si es:
   - un producto `A*B` con A y B no constantes;
   - un cuadrado `A^2`;
   - una división `A/B` con B no constante, porque fAF2 la calcula como `A * (1/B)`.

   A y B pueden ser subexpresiones compuestas, por ejemplo `sin(x)*y`, `x*ln(x)` o `x/(x+y)`: en aritmética afín
   cada factor tiene parte lineal en ε, así que el producto genera términos bilineales igual que `x*y`.
3. Cada factor se evalúa en aritmética afín: `A = a₀ + aᵀε + e_A ξ_A`.
   - exp, ln, sqrt y 1/x usan min-range.
   - sin, cos, `x^k` con k ≥ 3 y abs (cuando su argumento cruza el 0) usan la secante, con el residuo acotado
     por muestreo. Así conservan su parte lineal en ε.

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

`like` y `like-1`, y un nodo de ex8_5_2_1, dividen por expresiones cuyo intervalo contiene 0 en la raíz, así que
no se pueden medir.

## Resultados

Se encontraron **113 nodos suma con 2 o más términos cuadráticos, en 38 de las 41 instancias**. En 59 de ellos
los términos comparten algún símbolo. Las divisiones aportan 20 de esos nodos.

Ancho relativo a SH, primero media geométrica sobre los nodos de cada instancia y después sobre las instancias
(cada instancia pesa 1):

| | S | **SH** | Qk0 | Qk2 | Qk6 | EIG | BEST | EXACT |
|---|---|---|---|---|---|---|---|---|
| 38 instancias | 1,05 | 1 | 1,00 | 0,97 | 0,95 | 1,69 | 0,93 | 0,92 |
| 27 con símbolos compartidos | 1,06 | 1 | 0,99 | 0,95 | 0,93 | 1,69 | 0,90 | 0,89 |

- **SH ya es exacto en 86 de los 113 nodos.** El margen que queda es, en media geométrica, de un 8 % del ancho, y
  se concentra en pocas instancias:

  | instancia | Qk0 | Qk6 | EIG | EXACT | de dónde sale la ganancia |
  |---|---|---|---|---|---|
  | ex6_2_6 | 0,35 | 0,35 | 0,77 | 0,35 | productos idénticos con coeficientes que se cancelan (`x2*ln(x2)` aparece 4 veces y sus coeficientes suman ≈ 0) y el mismo `ln(…)` en varios productos |
  | ex6_2_8 | 0,40 | 0,40 | 0,97 | 0,40 | ídem |
  | ex2_1_9, ex2_1_9fact | 1 | 1 | **0,76** | 0,73 | suma de 22 bilineales `xᵢxⱼ` con diagonal cero: el k-nomial SDSOS no forma bloques, la cota espectral sí aprovecha |
  | ex8_5_1bis, ex8_5_2_1 | 1 | 0,92–0,95 | 1,1–1,5 | 0,82–0,88 | k-nomial; con la espectral, 0,86–0,90 |
  | ex14_2_3 | 1,10 | **0,94** | 1,57 | 0,89 | sumas de cocientes `xᵢ/(Σ cⱼxⱼ)` con el mismo denominador: el k-nomial sí gana a SH y a Qk0 |
  | ex14_2_7 | 1,15 | 0,99 | 1,63 | 0,90 | ídem; el k-nomial casi empata con SH y queda un 10 % de margen |

- **El k-nomial es lo decisivo solo en ex8_5_\* y ex14_2_\*, y gana poco: entre 1 y 8 %.**
  - En ex6_2_6 y ex6_2_8 basta con juntar los términos iguales (Qk0), lo que hace también una simplificación
    simbólica.
  - En ex2_1_9 el k-nomial no sirve, porque necesita diagonal.
  - En ex14_2_\* las sumas son de cocientes con el mismo denominador, que es un caso con estructura de verdad: Qk0
    pierde frente a SH, el k-nomial lo recupera y gana algo, y el exacto está un 10 % más abajo.
- **El k-nomial agregado es peor que SH en 21 de los 113 nodos.** Ocurre sobre todo en ex6_1_3bis, ex8_4_4bis y
  3pk-borne, con hasta 1,20× SH. Por eso, si se implementa, hay que quedarse con la menor
  de las cotas, y SH tiene que estar entre ellas.
- **La cota espectral sola es mala en instancias reales** (1,69× SH en media): las cajas no son cubos y Q no está
  normalizada. Solo sirve como una cota más dentro del mínimo.

## Comparación con synq

| | synq (78) | reales (36) |
|---|---|---|
| SH / exacto | 6,0 | 1,09 |
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
