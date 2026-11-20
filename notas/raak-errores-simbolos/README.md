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

## Presupuesto de error: de dónde viene el ancho

[`atribuir.py`](atribuir.py) descompone el radio del encierro de M1 en dos partes: la parte lineal en las
variables y la contribución de cada símbolo de error, ya propagado hasta el final. Cada símbolo lleva la etiqueta
de la operación que lo creó:

| etiqueta | operación |
|---|---|
| `mul_lin` | producto de factores afines sin error |
| `mul_comp` | producto con al menos un factor que ya trae error: el factor es una expresión no lineal |
| `sqr_lin`, `sqr_comp` | lo mismo para cuadrados |
| `recip`, `pow`, `exp`, `log`, `sqrt`, `sincos`, `abs` | funciones |

Hay dos salvedades:
- La parte lineal es en buena medida el rango verdadero de la expresión, no sobreestimación.
- Lo relevante es el reparto en las expresiones que sí están sobreestimadas.

Tablas completas en [`resumen_atribucion.txt`](resumen_atribucion.txt).

**Fracción del radio (%) en expresiones con M0 / INNER ≥ 1,5** (media sobre instancias):

| conjunto | prof | inst | lineal | mul_lin | **mul_comp** | sqr_comp | pow | exp | log | sincos |
|---|---|---|---|---|---|---|---|---|---|---|
| optimización | 0 | 21 | 31 | 5 | **45** | 0 | 2 | 1 | 15 | 0 |
| | 10 | 22 | 26 | 2 | **37** | 0 | 1 | 1 | 28 | 4 |
| | 20 | 20 | 39 | 1 | **26** | 0 | 1 | 2 | 26 | 3 |
| NCSP | 0 | 13 | 44 | 0 | **46** | 9 | 0 | 0 | 0 | 1 |
| | 10 | 12 | 47 | 1 | **41** | 6 | 0 | 5 | 0 | 1 |
| | 20 | 11 | 48 | 0 | **37** | 2 | 0 | 12 | 0 | 0 |

Sin contar la parte lineal, el origen dominante de cada instancia se reparte así:

| origen dominante | instancias |
|---|---|
| `mul_comp` | 23 |
| `sqr_lin` | 10 |
| `mul_lin` | 10 |
| `log` | 10, casi todas de ex6_2_\* |
| `sqr_comp` | 2 |
| `pow` | 2 |
| `exp` | 1 |

### Lectura

1. **La fuente principal son los productos de factores que ya traen error** (`mul_comp`): entre el 26 y el 46 % del
   radio en las expresiones sobreestimadas, en los dos conjuntos. Son productos de dos expresiones no lineales, es
   decir, polinomios de grado mayor que 2 o productos de funciones.
2. **Cada producto ya se encierra exacto** (SH sobre todos los símbolos, que equivale a E). Lo que se pierde es la
   dependencia: el error de un factor se trata como un símbolo independiente de las variables, aunque en realidad
   es una función de ellas, típicamente cuadrática. Ninguna mejora de la multiplicación de dos formas afines puede
   recuperar eso.
3. **Los productos y cuadrados de factores lineales** (`mul_lin`, `sqr_lin`) dominan en 20 instancias, pero esas
   instancias casi no están sobreestimadas (M0 / INNER entre 1,03 y 1,15). SH ya resuelve ese caso.
4. **ln aparece solo en ex6_2_\***, donde los términos `x·ln(x)` sobre [10⁻⁶, 1] tienen linealizaciones muy malas.
   Es un caso particular.
5. **exp crece con la profundidad en NCSP** (12 % a profundidad 20, por Trigexp2), pero es secundario.

### Qué sugiere

Atacar la pérdida de dependencia en los productos de factores no lineales. La forma general de hacerlo es una
aritmética que **guarde los términos de segundo orden** `εᵀQε` en la forma (formas afines cuadráticas, como
Messine y Touhami, o modelos de Taylor de orden 2), en vez de convertir cada parte no lineal en un símbolo nuevo.
Con eso, el error de un factor sigue siendo una función cuadrática de las variables, y el producto siguiente
puede aprovecharlo.

El costo es O(n²) por operación. El k-nomial, los bloques producto o la SDP volverían a tener sentido, porque
acotarían esa Q al final, una vez por expresión. Antes de implementarlo se puede medir con el mismo esquema: un
modo M3 con formas cuadráticas, truncando lo de grado mayor que 2 a un símbolo, en las mismas cajas.
