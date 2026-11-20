# Medición en la raíz de synq: residuo cuadrático del objetivo

Primer paso de la sección 4 de [`../raak-cambio-de-enfoque.md`](../raak-cambio-de-enfoque.md). En la caja raíz
`[-1,1]ⁿ` de las 78 instancias synq (`JCGeraldo/ibex-lib`, rama `raak-experimentos`, `raak/instancias/sinteticas`),
se mide el ancho del encierro de la parte cuadrática del objetivo.

## Qué se mide

El objetivo de synq es `f(x) = Σₖ (aₖᵀx)(bₖᵀx) + gᵀx`, con K = n productos de formas lineales sin constante. En la
raíz, ε = x, así que la parte lineal es `gᵀx` para todos los métodos y el residuo es `xᵀQx`, con
`Q = Σₖ ½(aₖbₖᵀ + bₖaₖᵀ)` de rango n.

| columna | qué es |
|---|---|
| `S` | suma por producto de la fórmula de fAF2 (por omisión en Ibex) |
| `SH` | suma por producto del rango exacto de `(aᵀx)(bᵀx)`. Aquí SH = E, porque los centros son 0 y no hay errores. **Línea base.** |
| `RAAp0`, `RAAp6` | suma por producto de `knomial_quad_bound` con k = 0 y k = 6 (RAA-K actual) |
| `Qk0` | Q completa, término a término. Es lo que daría la aritmética afín sobre el objetivo **expandido en monomios**. |
| `Qk2` … `Qk8` | Q completa con la reescritura k-nomial SDSOS de `ibex_KNomial.h` (sin modificar) y el residuo término a término |
| `EIG` | cota espectral `[n·min(λmin,0), n·max(λmax,0)]` |
| `EXACT` | rango exacto de `xᵀQx` en la caja, por enumeración de las caras KKT (3ⁿ) |

Las 78 instancias pasan dos comprobaciones: la Q armada desde el `.bch` coincide con la `Q_efectiva` del
`manifest.csv`, y ninguna cota es más angosta que el rango exacto.

## Resultados (media geométrica del ancho)

Ancho relativo a **SH** (suma por producto):

| grupo | inst | S | SH | RAAp6 | EIG | Qk0 | Qk2 | Qk6 | EXACT |
|---|---|---|---|---|---|---|---|---|---|
| todas | 78 | 1,90 | 1 | 1,21 | 0,23 | 0,25 | 0,22 | **0,21** | 0,17 |
| n = 3 | 18 | 1,80 | 1 | 1,16 | 0,38 | 0,40 | 0,34 | 0,33 | 0,30 |
| n = 5 | 24 | 1,87 | 1 | 1,21 | 0,24 | 0,25 | 0,22 | 0,21 | 0,17 |
| n = 8 | 18 | 1,96 | 1 | 1,22 | 0,18 | 0,20 | 0,18 | 0,17 | 0,13 |
| n = 10 | 18 | 1,96 | 1 | 1,23 | 0,16 | 0,19 | 0,17 | 0,17 | 0,12 |

Ancho relativo al **exacto**:

| grupo | SH | EIG | Qk0 | Qk2 | Qk6 | min(EIG, Qk6) |
|---|---|---|---|---|---|---|
| todas | 6,0 | 1,36 | 1,49 | 1,31 | 1,27 | 1,19 |
| n = 3 | 3,4 | 1,28 | 1,34 | 1,13 | 1,11 | 1,10 |
| n = 10 | 8,7 | 1,39 | 1,67 | 1,51 | 1,46 | 1,25 |
| d = 0,2 | 7,2 | 1,51 | 1,29 | 1,11 | 1,11 | — |
| d = 0,9 | 5,7 | 1,27 | 1,78 | 1,60 | 1,51 | — |

Las tablas completas, por densidad y por convexidad, están en [`resumen.txt`](resumen.txt), y los datos por
instancia en [`synq_raiz.csv`](synq_raiz.csv).

## Lectura

1. **El criterio de la nota (k-nomial sobre Q completa al menos 10–20 % más angosto que SH) se cumple de sobra:**
   Qk6 / SH = 0,21, es decir, casi 5 veces más angosto, y en las 78 instancias.
2. **Pero casi toda esa ganancia viene de juntar Q, no del k-nomial.** synq se genera como `Q = AᵀM` con A
   aleatoria (`generar_sinteticas.py`, `factor`), así que los productos se cancelan mucho entre sí. SH, al
   encerrar cada producto por separado, pierde esa cancelación y queda a 6× del exacto. Basta con expandir y
   acotar término a término (Qk0) para bajar a 0,25 de SH. La línea base justa para el k-nomial sobre Q completa
   es **Qk0**, no SH.
3. **Ganancia propia del k-nomial sobre Qk0:** Qk6 / Qk0 = 0,85 (entre 0,63 y 1,00). Cierra en promedio el 37 % de
   la brecha entre Qk0 y el exacto. Aun así queda a 1,27× del exacto, y el margen restante crece con n y con la
   densidad.
4. **k ≥ 4 no aporta nada en synq.** Qk4 = Qk6 = Qk8 en todas las instancias, y los bloques extraídos tienen a lo
   más 4 variables (2 en 52 instancias, 3 en 23 y 4 en 3). La búsqueda solo arma cliques balanceados en signo con
   diagonales del mismo signo, y eso escasea.
5. **La cota espectral, mucho más simple, compite con el k-nomial.** EIG / SH = 0,23 frente a 0,21. EIG es mejor
   en 30 de 78 instancias (en 10 de 18 con n = 10) y en las densas (d = 0,9). Tomar el mínimo de las dos baja de
   1,27 a 1,19 veces el exacto.

## Qué sugiere para el siguiente paso

- El argumento a favor de reenfocar está en **acotar la Q agregada**. Que el k-nomial sea la mejor manera de
  hacerlo no está demostrado: la línea base 4 de la nota (espectral o SDP) ya empata con él. Antes de integrar el
  k-nomial en un linearizador, conviene comparar con una cota SDP por diagonal,
  `max xᵀQx ≤ min { Σ dᵢ : diag(d) − Q ⪰ 0 }`, que suele ser muy ajustada para n ≤ 10.
- synq favorece mucho a "juntar Q" por construcción. Para no sobrevender el resultado, hay que medir también en
  instancias con productos que **no** vengan de una factorización artificial, y en cajas internas del
  branch-and-bound, no solo en la raíz.

## Reproducir

```bash
# instancias y cabecera del k-nomial
git clone --depth 1 -b raak-experimentos https://github.com/JCGeraldo/ibex-lib jcgeraldo-ibex-lib
g++ -O2 -std=c++11 -Istub -Ijcgeraldo-ibex-lib/raak/affine knomial_driver.cpp -o knomial_driver
python3 medir_synq.py jcgeraldo-ibex-lib/raak/instancias/sinteticas ./knomial_driver synq_raiz.csv
python3 resumir.py synq_raiz.csv > resumen.txt
```

- `knomial_driver.cpp` llama a `knomial_rewrite` (Q completa) y a `knomial_quad_bound` (por producto) de
  `ibex_KNomial.h` sin modificarlos. La cota de la Q completa replica la de `knomial_quad_bound` en modo 1 (dobles
  con redondeo hacia +∞).
- `stub/ibex_Interval.h` existe solo para compilar fuera de Ibex; con el modo 1 no se ejecuta.
- Requiere numpy. Tarda unos 3 s.
