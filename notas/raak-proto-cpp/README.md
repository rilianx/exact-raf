# Prototipo completo en C++: formas afines de grado 2 con encierro final estructurado

Compara, **todo en C++** y en las mismas expresiones y cajas, la aritmética de Ibex por omisión con la propuesta:
formas afines de grado 2 y un encierro final que aprovecha la estructura de la forma. Está fuera de Ibex y en
dobles, sin redondeo dirigido.

| archivo | contenido |
|---|---|
| [`emitir_cintas.py`](emitir_cintas.py) | convierte cada expresión Minibex en una cinta de operaciones (DAG con subexpresiones compartidas y constantes resueltas) y escribe las mismas cajas de las mediciones (semilla crc32) |
| [`aa2.cpp`](aa2.cpp) | evaluador; usa `ibex_KNomial.h` (k-nomial) y `../raak-prodk-cpp/prodk.h` (bloques producto) |
| [`resumir.py`](resumir.py) | tabla final; comprueba que ningún encierro quede por debajo del rango interior |
| [`proto_cpp.csv`](proto_cpp.csv), [`resumen.txt`](resumen.txt) | resultados |

**Modos:**

| modo | qué es |
|---|---|
| `M00` | fAF2 de Ibex por omisión, con la fórmula (5) de Skalna-Hladik en el producto |
| `M0` | fAF2 con el producto SH (parches 0016/0017) |
| `M3c` | formas `c + aᵀz + zᵀQz`, con z = variables más símbolos de error, y resto de grado 3 y 4 centrado |

**Encierros finales de M3c:**
- **barato:** cada símbolo exacto con su parte lineal y su diagonal, más el resto término a término;
- **k-nomial** k = 6;
- **bloques producto, pares + completo:** pares y soporte completo si m ≤ 12;
- **bloques producto completos:** pares, triples y soporte completo.

**Funciones** (iguales en todos los modos): min-range para exp, ln, sqrt y 1/x; secante con residuo exacto por
puntos críticos para sin, cos, abs y `x^k`. Contra la versión en Python (que muestreaba el residuo), los anchos
coinciden en las instancias polinomiales y difieren menos de un 0,5 % con sin y cos.

```bash
B=jcgeraldo-ibex-lib   # rama raak-experimentos
python3 emitir_cintas.py cintas.txt <lista de .bch>    # la misma lista de ../raak-errores-simbolos
g++ -O2 -std=c++14 -I../raak-medicion-synq/stub -I$B/raak/affine -I../raak-prodk-cpp aa2.cpp -o aa2
./aa2 cintas.txt 5 > proto_cpp.csv                     # 5 repeticiones por medición; unos 6 s en total
python3 resumir.py > resumen.txt
```

## Resultados

7532 combinaciones de expresión y caja en 58 instancias: 41 de optimización y 17 de NCSP, raíz más cajas
internas a profundidad 5, 10 y 20. Ningún encierro queda por debajo del rango interior. Las precisiones son
medias geométricas, por instancia y después entre instancias. Los tiempos son por expresión y caja: evaluación
más encierro final.

| método | ancho/M00, optim. | ancho/M00, NCSP | instancias que mejoran > 2 % (de 58) | ancho/INNER, optim. | ancho/INNER, NCSP | tiempo, mediana / p90 (µs) | tiempo / M00 (mediana) |
|---|---|---|---|---|---|---|---|
| fAF2 de Ibex (fórmula S) | 1 | 1 | — | 1,76 | 1,59 | 0,43 / 8,2 | 1 |
| fAF2 + SH | 0,987 | 0,947 | 12 | 1,74 | 1,50 | 0,57 / 8,4 | 1,1 |
| grado 2 + encierro barato | 0,950 | 0,955 | 17 | 1,67 | 1,51 | 0,79 / 14,8 | 1,5 |
| **grado 2 + k-nomial k = 6** | **0,946** | **0,922** | 17 | 1,66 | 1,46 | **1,35 / 18,2** | **2,2** |
| grado 2 + bloques producto (pares + completo) | 0,934 | 0,914 | 27 | 1,64 | 1,45 | 5,0 / 58 | 8,0 |
| **grado 2 + bloques producto (completo)** | **0,925** | **0,905** | **32** | 1,63 | 1,43 | 8,2 / 222 | 12 |

**Por profundidad** (ancho / M00, las 58 instancias juntas):

| prof | fAF2 + SH | grado 2 + barato | + k-nomial | + bloques (pares + completo) | + bloques (completo) |
|---|---|---|---|---|---|
| 0 | 0,964 | 0,966 | 0,942 | 0,934 | 0,926 |
| 5 | 0,971 | 0,949 | 0,934 | 0,923 | 0,913 |
| 10 | 0,974 | 0,945 | 0,933 | 0,922 | 0,911 |
| 20 | 0,982 | 0,958 | 0,950 | 0,938 | 0,931 |

## Lectura

1. **Grado 2 con k-nomial cuesta unas 2 veces fAF2** (1,35 µs frente a 0,43 µs de mediana por expresión) y quita
   un 5 % del ancho en optimización y un 8 % en NCSP.
2. **Grado 2 con bloques producto completos cuesta unas 12 veces fAF2** (8 µs de mediana) y quita un 7,5 % en
   optimización y un 9,5 % en NCSP. Mejora 32 de las 58 instancias.
3. **Pares + completo es el punto intermedio:** cuesta 8 veces fAF2, con colas mucho más cortas que la versión
   completa (p90 58 µs frente a 222 µs).
4. **La ganancia se mantiene en las cajas internas:** entre un 6 y un 9 % a profundidad 5–10, y menos a
   profundidad 20, donde la caja ya es chica y todo se acerca a lineal.
5. **Los tiempos son de un prototipo sin optimizar para Ibex.** Los vectores dispersos de M3c y los `std::map` de
   la homogeneización se pueden mejorar. Aun así, el orden de magnitud está claro: el k-nomial es casi gratis y
   los bloques producto cuestan unas pocas evaluaciones de la expresión.

## Qué falta para Ibex

- **Rigor en punto flotante.** Hay que acumular los errores de redondeo de M3c y recalcular con intervalos el
  residuo del encierro final, con el mismo esquema que `knomial_quad_bound` (parches 0012 y 0013).
- **Integración como aritmética** (o como modo de fAF2) y evaluación en el solver: número de celdas y tiempo total.
  Que un encierro sea más angosto no garantiza menos celdas.
- **Regla para elegir el encierro final** según el tamaño m de la forma: k-nomial siempre, y bloques producto
  cuando m es chico.
