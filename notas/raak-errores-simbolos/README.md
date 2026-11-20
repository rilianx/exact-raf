# Errores como símbolos y sumas de productos agregadas, en cajas internas

Pregunta: en sistemas no lineales genéricos, ¿cuánto se gana si los errores de las subexpresiones no lineales se
guardan como **símbolos propios** en vez del escalar de fAF2, y si además las **sumas de productos** se acotan en
conjunto?

Es una medición estática sobre la raíz y cajas internas, sin correr el solver.

**Instancias (58):**
- **Optimización (41):** las 38 medium+hard de Ibex y 3 de COCONUT (`3pk-borne`, `mistake`, `s365mod`).
- **NCSP (17):** las mismas del barrido 20 de JCGeraldo (`benchs/solver`): Directkin, 12 Trigexp2, Dietmaier,
  Neveu1, Pramanik y ponts-geo.

## Evaluadores ([`aa_sym.py`](aa_sym.py))

| modo | qué hace |
|---|---|
| `M0` | lo de hoy: fAF2 con el error como un escalar por forma, y SH en cada producto |
| `M1` | cada operación no lineal (producto, función, potencia) crea un símbolo de error nuevo, que se propaga linealmente; SH en cada producto sobre todos los símbolos |
| `M2` | M1, y cada suma de 2 o más productos se acota en conjunto con la intersección de tres cotas: suma de SH, término a término de la Q agregada y bloques producto |

Detalles comunes a los tres modos:
- Las subexpresiones repetidas se evalúan una sola vez por expresión, como en el DAG de Ibex, así que en M1 y M2
  comparten sus símbolos.
- Las linealizaciones unarias son las mismas en los tres: min-range para exp, ln, sqrt y 1/x, y secante para sin,
  cos, `x^k` y abs.
- La referencia `INNER` es el rango interior de la expresión, obtenido por muestreo y búsqueda local. Es una cota
  inferior del rango verdadero.

**Cajas ([`medir_cajas.py`](medir_cajas.py)):** la raíz, más 4 caminatas aleatorias de bisección por instancia
(se biseca la variable de mayor ancho relativo y se elige un hijo al azar), con las cajas a profundidad 5, 10 y 20.
En cada caja se mide el ancho del encierro de cada objetivo y restricción.

**Validez:** ningún encierro de M0, M1 o M2 queda por debajo del rango interior en ninguna caja.

## Resultados

Media geométrica, primero sobre las expresiones y cajas de cada instancia y después sobre las instancias:

| conjunto | prof | M1 / M0 | M2 / M0 | INNER / M0 | M0 / INNER | M2 / INNER |
|---|---|---|---|---|---|---|
| optimización (41) | 0 | 0,99 | 0,975 | 0,40 | 2,50 | 2,44 |
| | 5 | 0,99 | 0,954 | 0,51 | 1,96 | 1,87 |
| | 10 | 0,99 | 0,952 | 0,61 | 1,65 | 1,57 |
| | 20 | 0,99 | 0,966 | 0,74 | 1,35 | 1,31 |
| NCSP (17) | 0 | 0,98 | 0,965 | 0,62 | 1,62 | 1,56 |
| | 5 | 0,99 | 0,970 | 0,62 | 1,61 | 1,57 |
| | 10 | 0,99 | 0,978 | 0,68 | 1,48 | 1,44 |
| | 20 | 1,00 | 0,994 | 0,74 | 1,34 | 1,33 |

**Instancias que cambian** (todas las profundidades juntas):

| instancia | M1 / M0 | M2 / M0 | INNER / M0 |
|---|---|---|---|
| ex6_2_6 | 0,77 | **0,46** | 0,13 |
| ex6_2_8 | 0,78 | **0,47** | 0,07 |
| Pramanik | 0,93 | **0,77** | 0,39 |
| Directkin | 0,94 | 0,91 | 0,59 |
| ex6_2_11 | 1,00 | 0,94 | 0,10 |
| ex6_2_9 | 1,00 | 0,97 | 0,07 |
| ex6_2_12 | 1,00 | 0,98 | 0,21 |

Solo 4 de las 58 instancias mejoran más de un 2 % con M1, y 7 con M2. Tabla completa en [`resumen.txt`](resumen.txt)
y datos en [`cajas.csv`](cajas.csv).

## Lectura

1. **En promedio la ganancia es chica:** con M2 se quita entre un 2 y un 5 % del ancho, y con M1 alrededor de un
   1 %.
2. **Es grande en pocas instancias, y son las de siempre.** ex6_2_\* es la familia donde RAA-K tenía material.
   Pramanik y Directkin son donde SH y E ya habían ganado en el solver. Ahí M2 quita entre un 23 y un 54 %, y la
   ganancia se mantiene en las cajas internas.
3. **M2 aporta más que M1.** Guardar los errores como símbolos casi no sirve por sí solo. Sirve cuando hay una
   operación que aprovecha las correlaciones, como la suma de productos agregada.
4. **La mayor parte de la sobreestimación no está en lo bilineal.** M0 está entre 1,3 y 2,5 veces por encima del
   rango interior, y M2 baja eso muy poco: de 2,50 a 2,44 en la raíz de optimización. El resto viene de las
   linealizaciones de funciones (exp, ln, sin sobre intervalos anchos), de potencias y cocientes, y de la
   dependencia entre subexpresiones que la aritmética afín no captura. La parte bilineal es una fracción menor del
   problema en estos benchmarks.

## Reproducir

```bash
B=jcgeraldo-ibex-lib   # rama raak-experimentos
python3 medir_cajas.py cajas.csv $B/benchs/optim/medium/*.bch $B/benchs/optim/hard/*.bch \
    $B/raak/instancias/reales/{3pk-borne,mistake,s365mod}.bch $(sed "s#^#$B/#" $B/raak/results/20-benchs-solver/csp-900s/instancias.txt)
python3 resumir.py cajas.csv > resumen.txt
```

Usa `../raak-qplib/prodk_fast.py` para los bloques producto. Tarda unos 5 minutos.
