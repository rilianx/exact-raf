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

## M3: formas cuadráticas ([`aa_quad.py`](aa_quad.py), [`medir_m3.py`](medir_m3.py))

**Qué cambia respecto de M1:** M3 guarda los términos de segundo orden. Cada forma es `c + aᵀz + zᵀQz` (más
símbolos de error):
- **Producto:** se conserva todo el grado 2. Lo de grado 3 y 4 se acota con intervalos y pasa a un símbolo nuevo.
- **Funciones:** usan la misma linealización que M1, sobre el rango de la forma completa.

**Dos variantes de encierro final:**
- **M3a (barato):** cada símbolo con su parte lineal y su diagonal, exacto en una dimensión, más el resto término
  a término.
- **M3b:** la menor entre M3a, los bloques producto y la SDP de Shor sobre la forma homogeneizada.

**Medición:** todos los modos se midieron en las **mismas cajas**, con la semilla determinista (crc32). La
atribución de arriba se rehízo con esas cajas y no cambia ([`resumen_atribucion.txt`](resumen_atribucion.txt)).
Ningún encierro quedó por debajo del rango interior en ninguna de las 7534 combinaciones de expresión y caja.

**Error encontrado y corregido.** Con esta medición apareció un error en el rango SH: el orden angular con
`atan2` empataba cuando u y v tenían escalas muy distintas, y el rango salía demasiado angosto. Se corrigió en
`../raak-qplib/prodk_fast.py`, con producto cruz y normalización, y se recalcularon todas las mediciones que lo
usaban. Solo cambiaron un poco los pares espectrales de `../raak-pares`, y las conclusiones se mantienen.

**Ancho / M0** (media geométrica; `M2|M3b` es la intersección de M2 y M3b, ambos válidos):

| conjunto | prof | M1 | M2 | M3a | M3b | M2\|M3b | M0 / INNER | (M2\|M3b) / INNER |
|---|---|---|---|---|---|---|---|---|
| optimización (41) | 0 | 0,99 | 0,975 | 0,99 | 0,92 | **0,92** | 2,50 | 2,31 |
| | 5 | 0,99 | 0,96 | 0,97 | 0,91 | **0,91** | 2,01 | 1,83 |
| | 10 | 0,99 | 0,96 | 0,96 | 0,91 | **0,91** | 1,74 | 1,58 |
| | 20 | 0,99 | 0,97 | 0,96 | 0,93 | **0,93** | 1,40 | 1,31 |
| NCSP (17) | 0 | 0,98 | 0,965 | 1,04 | 0,98 | **0,96** | 1,62 | 1,55 |
| | 5 | 0,99 | 0,97 | 1,01 | 0,95 | **0,94** | 1,61 | 1,51 |
| | 10 | 1,00 | 0,98 | 1,01 | 0,94 | **0,94** | 1,50 | 1,41 |
| | 20 | 1,00 | 0,99 | 1,01 | 0,95 | **0,95** | 1,38 | 1,31 |

Por instancia, en [`resumen_m3.txt`](resumen_m3.txt). Las que más bajan con `M2|M3b`:

| instancias | ancho / M0 |
|---|---|
| ex6_2_6, ex6_2_8 | 0,45–0,55 |
| Pramanik | 0,66 |
| s365mod | 0,77 |
| Directkin | 0,83 |
| 3pk-borne, ex6_1_3, ex6_2_11, ex6_1_4, Dietmaier, ex2_1_9 | 0,84–0,89 |

### Lectura de M3

1. **M3 mejora más instancias, pero recupera poco de la brecha.** Con `M2|M3b`, 37 de las 58 instancias mejoran
   más de un 2 %, contra 6 con M2. Aun así, la brecha frente al rango interior solo se cierra entre un 13 y un
   22 %: por ejemplo, en optimización a profundidad 10 pasa de 1,74 a 1,58.
2. **La ganancia viene sobre todo del encierro final fuerte (M3b), no de la aritmética.**
   - M3a, con el encierro barato, queda en 0,96–1,04 y es **peor que M1 en el 10 %** de los casos. El resto de
     grado 3 y 4 acotado con intervalos pierde más de lo que gana el grado 2.
   - En las instancias cuadráticas (ex2_1_9, ex5_3_2, ex6_1_4) M3b es justamente la cota SDP de la Q completa, que
     ya habíamos visto en QPLIB.
3. **Costo.** M3 es O(n²) por operación, y el encierro final de M3b (SDP y bloques en Python) tiene una mediana
   de 0,03 s por expresión y caja. Es demasiado para hacerlo en cada nodo.
4. **Queda mucha sobreestimación sin explicar por lo bilineal:** `M2|M3b` sigue entre 1,3 y 2,3 veces el rango
   interior.

### Qué haría falta para que M3 valga la pena

- **Acotar mejor lo de grado 3 y 4 en el producto de formas cuadráticas**, en vez de intervalos sueltos. Por
  ejemplo, acotar `(aᵀz)(zᵀQz)` con la forma conjunta. Sin eso, M3 barato no le gana a M1.
- **Un encierro final barato pero fuerte:** un certificado tipo SDP reutilizable, o el k-nomial y los bloques
  producto en C++, en lugar de la SDP en cada evaluación. Falta medir M3 solo con bloques producto, sin la SDP,
  para separar esos dos efectos.

## M3 con resto centrado y encierros finales separados ([`medir_m3c.py`](medir_m3c.py), [`resumir_m3c.py`](resumir_m3c.py))

**Resto centrado (`M3c`).** En el producto, cada parte cuadrática se separa como `q = m + q̂`: el centro m de su
rango pasa a la constante, que es exacta, y el resto de grado 3 y 4 se acota con q̂ ∈ [−r, r]. En los cuadrados,
el término `q̂²` queda en [0, r²]. `M3o` es la versión original, que acotaba el resto con el rango completo de q.

**Encierros finales**, calculados por separado sobre la misma forma:

| sufijo | encierro |
|---|---|
| `_tt` | barato: cada símbolo con su parte lineal y su diagonal, exacto, más el resto término a término |
| `_prodk` | la menor entre `_tt` y los bloques producto |
| `_kn6` | la menor entre `_tt` y el k-nomial con k = 6 de `ibex_KNomial.h`, en C++ |
| `_sdp` | la menor entre `_tt` y la SDP de Shor |

Medido en las mismas 7534 combinaciones de expresión y caja. Ninguna queda por debajo del rango interior. Tabla
completa en [`resumen_m3c.txt`](resumen_m3c.txt).

**Ancho / M0** (media geométrica):

| conjunto | prof | M2 | M3o_tt | M3c_tt | M3c_prodk | M3c_kn6 | M3c_sdp | M2\|M3c_sdp | M0 / INNER → (M2\|M3c_sdp) / INNER |
|---|---|---|---|---|---|---|---|---|---|
| optimización | 0 | 0,975 | 0,991 | 0,990 | 0,965 | 0,977 | 0,927 | 0,923 | 2,50 → 2,31 |
| | 10 | 0,959 | 0,958 | 0,957 | 0,941 | 0,953 | 0,914 | 0,912 | 1,74 → 1,58 |
| | 20 | 0,968 | 0,962 | 0,962 | 0,950 | 0,961 | 0,937 | 0,936 | 1,40 → 1,31 |
| NCSP | 0 | 0,965 | 1,043 | 1,032 | 0,981 | 0,978 | 0,967 | 0,958 | 1,62 → 1,55 |
| | 10 | 0,979 | 1,005 | 1,004 | 0,957 | 0,970 | 0,946 | 0,945 | 1,50 → 1,42 |
| | 20 | 0,991 | 1,010 | 1,010 | 0,965 | 0,983 | 0,965 | 0,965 | 1,38 → 1,33 |

**Instancias que bajan más de un 2 %:**

| variante | instancias |
|---|---|
| M2 | 6 |
| M3c_tt | 10 |
| M3c_kn6 | 12 |
| **M3c_prodk** | **27** |
| M3c_sdp | 28 |
| M2\|M3c_sdp | 30 |

### Lectura

1. **Centrar el resto ayuda poco en promedio**, aunque sí en algunas instancias: en Pramanik, M3_tt pasa de 0,77 a
   0,74. M3c_tt sigue siendo peor que M1 en el 10 % de los casos (773 frente a 783). Así que el resto de grado 3 y 4
   no era la causa principal.
2. **La pérdida está en el encierro final barato.** M1 encierra cada producto exacto con SH. M3_tt, en cambio,
   acota la Q agregada término a término y pierde esa exactitud. Con los bloques producto, que recuperan productos
   de dos sumas exactos, M3c baja a 0,94–0,965 en optimización y a 0,957–0,98 en NCSP, y mejora 27 instancias:
   casi lo mismo que la SDP (28), sin resolver ninguna SDP.
3. **El k-nomial (k = 6) aporta menos que los bloques producto** sobre estas formas: 12 instancias frente a 27. Los
   productos de factores no lineales generan formas con poca diagonal útil para cuadrados.
4. **La SDP sigue siendo la mejor**, con un 3 a 4 % más que los bloques producto, pero no se puede pagar por nodo.
5. **En conjunto, la vía bilineal y cuadrática reduce el ancho entre un 5 y un 9 % en promedio**, mejora unas 30 de
   las 58 instancias y llega a un 35–55 % en pocas (ex6_2_6, ex6_2_8, Pramanik, s365mod). La brecha frente al rango
   interior se cierra solo en una fracción: a profundidad 10, de 1,74 a 1,58 en optimización y de 1,50 a 1,42 en
   NCSP.

## Tiempos ([`tiempos.py`](tiempos.py), [`resumen_tiempos.txt`](resumen_tiempos.txt))

Se midieron 1158 combinaciones de expresión y caja (caminata 1, profundidades 0 y 10); 772 de ellas tienen parte
cuadrática. Todo está en Python salvo el k-nomial, que es C++. Los tiempos sirven para comparar entre sí los
componentes, no como valores absolutos.

| componente | mediana (µs) | p90 (µs) | mediana / evaluación M0 |
|---|---|---|---|
| evaluación M0 (fAF2 + SH) | 375 | 2404 | 1 |
| evaluación M1 | 257 | 2168 | 0,82 |
| evaluación M3 centrado | 225 | 2389 | 0,75 |
| encierro barato (`_tt`) | 6 | 16 | 0,02 |
| bloques producto (Python) | 14 503 | 57 566 | 18 |
| k-nomial k = 6 (C++, en lote) | 70 | — | 0,16 |
| SDP (SCS) | 62 674 | 255 116 | 32 |

**Tamaño de la forma homogeneizada** (número de símbolos + 1): mediana 5, p90 13 y máximo 102. Las formas son
chicas.

**Lectura:**
- Evaluar con formas cuadráticas (M3) no cuesta más que con M1 en estos tamaños.
- El encierro barato es despreciable.
- Los bloques producto, en este prototipo en Python, cuestan unas 18 evaluaciones. Una versión en C++ debería
  quedar en el orden del k-nomial, que es ~70 µs.
- La SDP cuesta unas 30 evaluaciones y solo es viable fuera del bucle de cada nodo.
