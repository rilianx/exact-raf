# Reescribir sumas de productos como productos de pares de sumas

Idea: antes de evaluar, reescribir una suma de productos como **pocos productos de dos sumas largas**
(`x*y + x*z → x*(y+z)`), para que SH o E, que son exactos para un producto de dos formas afines, hagan el trabajo.
No hace falta una aritmética nueva: es una reescritura de la expresión.

Se usa la misma Q de las mediciones anteriores, en la caja raíz y sin correr el solver:
- en synq, la Q del objetivo de las 78 instancias;
- en las reales, los 113 nodos suma con 2 o más términos cuadráticos de 38 instancias.

Se prueban dos reescrituras. En las dos, cada producto se encierra con su rango exacto (SH) y se suman los anchos.

| columna | reescritura |
|---|---|
| `SH` | los productos tal como vienen en la expresión (**línea base**) |
| `STAR` | factor común: `zᵢ·(Qᵢᵢzᵢ + 2Σⱼ Qᵢⱼzⱼ)`, eligiendo con avidez la variable de mayor peso restante. Es lo que hace una factorización simbólica. |
| `EIGP` | pares espectrales: cada autovalor positivo se empareja con uno negativo, `λ₊(pᵀz)² − |λ₋|(qᵀz)² = (a+b)(a−b)` |
| `EIGSQ` | todos los autovectores como cuadrados, sin emparejar |
| `PARES` | la menor entre SH, STAR y EIGP |
| `TODO` | la menor entre PARES, el k-nomial con k = 6 (`Qk6`) y la cota espectral (`EIG`) |

`Qk0`, `Qk6`, `EIG` y `EXACT` son los de [`../raak-medicion-synq`](../raak-medicion-synq) y
[`../raak-medicion-reales`](../raak-medicion-reales). Todas las reescrituras dan cotas válidas: ninguna queda por
debajo del exacto.

## Resultados

**synq** (media geométrica):

| | SH | STAR | EIGP | Qk0 | Qk6 | EIG | PARES | TODO | EXACT |
|---|---|---|---|---|---|---|---|---|---|
| ancho / SH | 1 | **0,22** | 0,27 | 0,25 | 0,21 | 0,23 | 0,22 | 0,19 | 0,17 |
| ancho / exacto, todas | 6,0 | **1,34** | 1,62 | 1,49 | 1,27 | 1,36 | 1,30 | 1,16 | 1 |
| ancho / exacto, n = 3 | 3,4 | 1,16 | **1,08** | 1,34 | 1,11 | 1,28 | 1,06 | 1,03 | 1 |
| ancho / exacto, n = 10 | 8,7 | 1,54 | 2,31 | 1,67 | 1,46 | 1,39 | 1,54 | 1,25 | 1 |

**Reales** (38 instancias; media geométrica por instancia, ancho / SH):

| SH | STAR | EIGP | Qk0 | Qk6 | EIG | PARES | TODO | EXACT |
|---|---|---|---|---|---|---|---|---|
| 1 | 0,97 | 1,18 | 1,00 | 0,95 | 1,69 | 0,945 | 0,928 | 0,920 |

La tabla por instancia está en [`resumen.txt`](resumen.txt).

## Lectura

1. **Factorizar y usar SH funciona casi tan bien como el k-nomial, y sin aritmética nueva.**
   - En synq, STAR baja a 0,22 del ancho de SH, mejor que acotar la Q completa término a término (0,25) y cerca
     del k-nomial con k = 6 (0,21).
   - En las reales, STAR reproduce las ganancias de ex6_2_6 y ex6_2_8 (0,35 y 0,40) y de ex8_5_\* (0,90–0,93).
   - Es solo SH sobre productos de dos sumas: O(n log n) por producto, y la implementación de SH y E en fAF2 ya
     existe.
2. **La reescritura no depende de la caja.** Si `xᵢ = mᵢ + rᵢεᵢ`, la parte cuadrática de `xᵢ·(Σ cⱼxⱼ)` es
   `rᵢεᵢ·Σ cⱼrⱼεⱼ`, sea cual sea el centro. Por eso la factorización simbólica `x*y + x*z → x*(y+z)` da en cada
   caja la misma descomposición que STAR en la raíz, y se puede hacer una sola vez, antes del solver. (La elección
   ávida usa los pesos de la raíz.) EIGP, en cambio, depende de los radios.
3. **Hay que elegir, no reescribir siempre.**
   - STAR es peor que la expresión original en 27 de los 113 nodos reales (himmel16 1,13; 3pk-borne 1,42).
   - EIGP lo es en 36 de 113.

   Como en la relajación lineal cada forma da una relajación válida, se pueden **agregar las dos**, la original y
   la reescrita, y quedarse en la práctica con la mejor sin tener que decidir.
4. **Combinado con el k-nomial y la espectral (TODO)** se llega a 1,16 veces el exacto en synq y a 0,928 de SH en
   las reales, contra 0,920 del exacto.
5. **En los benchmarks reales el margen sigue siendo chico.** SH sobre la expresión original ya es exacto en 87 de
   los 113 nodos, así que la reescritura en pares no puede mejorar mucho ahí. Su valor está en problemas con
   formas cuadráticas densas que no vienen factorizadas.

## Siguiente paso posible

STAR es una heurística ávida muy simple, y en synq queda a 1,34 veces el exacto. Cabe buscar mejores
descomposiciones en pares, por ejemplo:
- elegir el factor común y el agrupamiento minimizando la suma de los anchos de SH, con búsqueda local;
- combinar factores comunes con pares espectrales;
- permitir que el factor sea una suma, `(xᵢ + xⱼ)·(…)`.

Como la reescritura se puede hacer una vez, antes del solver, se puede gastar más tiempo en buscarla.

## Reproducir

```bash
B=jcgeraldo-ibex-lib
python3 medir_pares.py synq $B/raak/instancias/sinteticas ../raak-medicion-synq/synq_raiz.csv synq_pares.csv
python3 medir_pares.py reales ../raak-medicion-reales/nodos.csv reales_pares.csv \
    $B/benchs/optim/medium/*.bch $B/benchs/optim/hard/*.bch $B/raak/instancias/reales/{3pk-borne,mistake,s365mod}.bch
python3 resumir.py > resumen.txt
```

## Generalización del k-nomial: reducción por bloques producto (`PRODK`)

[`medir_prodk.py`](medir_prodk.py) extrae de Q, de forma ávida, bloques que son productos de dos sumas
`(uᵀz)(vᵀz)`, cada uno encerrado exacto con SH, y acota el residuo término a término. Es lo mismo que hace el
k-nomial, que solo extrae cuadrados (u = v).

- **Si u y v no comparten variables, no se gana nada:** el encierro de SH coincide con el término a término. La
  ganancia viene de que las sumas compartan variables.
- **No necesita diagonal en todas las variables del bloque**, a diferencia de un cuadrado.

**Candidatos por soporte.** Para pares, triples conectados y el soporte completo si hay a lo más 12 símbolos:
- el par espectral de la submatriz, `(a+b)(a−b)`;
- el autovector dominante, como cuadrado;
- filas de factor común;
- cada uno con fracciones 1, ¾, ½ y ¼.

Tabla completa en [`resumen_prodk.txt`](resumen_prodk.txt).

**synq** (ancho / exacto; `MIN2` es la menor entre PRODK y el k-nomial con k = 6):

| | SH | Qk6 | PRODK | MIN2 | la menor de todas | PRODK < Qk6 |
|---|---|---|---|---|---|---|
| todas | 6,01 | 1,27 | **1,22** | 1,21 | 1,14 | 58 / 78 |
| n = 3 | 3,37 | 1,11 | **1,01** | 1,01 | 1,01 | 14 / 18 |
| d = 0,2 | 7,16 | 1,11 | **1,04** | 1,03 | 1,03 | 21 / 24 |
| d = 0,9 | 5,74 | 1,51 | 1,46 | 1,44 | 1,22 | 15 / 24 |

**Reales** (ancho / SH, media por instancia): Qk6 0,952, PRODK 0,959, MIN2 0,944, la menor de todas 0,927,
exacto 0,920.

| instancia | Qk6 | PRODK | exacto |
|---|---|---|---|
| ex8_5_2_1 | 0,92 | **0,82** | 0,82 |
| ex8_5_1bis | 0,95 | **0,88** | 0,88 |
| ex14_2_3 | **0,94** | 1,00 | 0,89 |
| ex14_2_7 | **0,99** | 1,06 | 0,90 |
| ex6_2_6, ex6_2_8 | 0,35–0,40 | 0,35–0,40 | 0,35–0,40 |

**Lectura:**
- Los bloques producto le ganan al k-nomial en synq (58 de 78 instancias, sobre todo cuando Q es dispersa).
- En ex8_5_\* llegan al exacto.
- En ex14_2_\* pierden, porque los candidatos ávidos no encuentran la estructura que sí aprovecha el k-nomial.
- Los dos métodos se complementan: lo natural es extraer bloques de ambos tipos, cuadrados y productos, en la
  misma reducción.
- Con diagonal cero (ex2_1_9) ninguno de los dos sirve: cualquier producto con variables compartidas crea
  diagonal, y el residuo la paga. Ahí solo ayuda la cota espectral.

**Costo:** el prototipo en Python tarda alrededor de 1 s por Q con n = 10. Hay que reimplementarlo antes de pensar
en el solver.
