# Formas cuadráticas generales en la raíz: QPLIB y un generador ad hoc

synq factoriza `Q = AᵀM` y crea a propósito cancelación entre productos, así que no sirve como banco general. Esta
medición usa formas cuadráticas `xᵀQx` generales **escritas en monomios**, que es como aparecen en QPLIB y en los
modelos reales. Con la Q en monomios, la línea base de Ibex (fAF2) es el encierro término a término (`Qk0`).

Es una medición estática en la caja raíz, sin correr el solver.

## Instancias

- **QPLIB (58):** las mismas que evaluó JCGeraldo en `raak/results/03-qplib/` (rama `raak-experimentos`).
  - Se mide la Q del objetivo en la caja de las variables. Las binarias se relajan a [0, 1].
  - n va de 30 a 175, y 21 de las 58 tienen diagonal cero.
  - Los `.qplib` se obtuvieron de dos espejos públicos en GitHub (`ZhenglingYangli/NLIPSat-Toolkit` y
    `Esnilg/QPLIB`), porque qplib.zib.de no es accesible desde el entorno. No están en el repo.
  - El encierro término a término coincide con el `af2_mono` de JCGeraldo, salvo un factor 2 de convención.
- **Generador fqf (81):** [`gen_fqf.py`](gen_fqf.py), en [`instancias_fqf/`](instancias_fqf), con la Q y un
  `.bch` en monomios de cada instancia.
  - n = 10, 20 y 40.
  - Densidad fuera de la diagonal de 0,1, 0,3 y 1,0, con valores U[−1, 1].
  - Diagonal cero, positiva o con ambos signos.
  - 3 semillas por combinación.

## Métodos

Todos se aplican a `εᵀ(DQD)ε` en [−1,1]ⁿ, con D la matriz diagonal de los radios de la caja.

| columna | qué es |
|---|---|
| `Qk0` | término a término, lo que da Ibex hoy (**línea base**) |
| `Qk2`, `Qk3`, `Qk6` | k-nomial SDSOS de `ibex_KNomial.h` |
| `PRODK` | reducción por bloques producto ([`prodk_fast.py`](prodk_fast.py), versión incremental de `../raak-pares/medir_prodk.py`) |
| `BLOQ` | la menor entre `Qk6` y `PRODK` |
| `EIG` | cota espectral |
| `SDP` | relajación de Shor: `max xᵀQx ≤ min Σdᵢ` con `d ≥ 0` y `diag(d) − Q ⪰ 0`. Se resuelve con SCS y se repara para que d sea factible, así que la cota es válida. |
| `INNER` | rango interior por búsqueda local, que es una cota inferior del exacto. Con n = 10 es el exacto (enumeración de caras). |

## Resultados (media geométrica)

**Ancho / término a término (Qk0):**

| | Qk3 | Qk6 | PRODK | BLOQ | EIG | SDP | INNER |
|---|---|---|---|---|---|---|---|
| QPLIB (58) | 0,97 | 0,97 | 0,98 | 0,97 | 0,41 | **0,34** | 0,30 |
| QPLIB, diagonal cero (21) | 1,00 | 1,00 | 1,00 | 1,00 | 0,52 | 0,47 | 0,43 |
| QPLIB, diagonal no nula (37) | 0,96 | 0,96 | 0,97 | 0,96 | 0,35 | 0,28 | 0,25 |
| fqf (81) | 0,85 | 0,85 | 0,89 | 0,84 | 0,96 | **0,56** | 0,53 |
| fqf, diagonal cero | 1,00 | 1,00 | 1,00 | 1,00 | 0,85 | 0,63 | 0,58 |
| fqf, diagonal positiva | 0,73 | 0,73 | 0,84 | 0,73 | 0,85 | 0,53 | 0,51 |
| fqf, diagonal con ambos signos | 0,84 | 0,84 | 0,85 | 0,81 | 1,21 | 0,52 | 0,51 |

**Ancho / INNER** (distancia al exacto, sobreestimada un poco porque INNER es una cota inferior):

| | Qk0 | BLOQ | EIG | SDP |
|---|---|---|---|---|
| QPLIB | 3,31 | 3,22 | 1,34 | **1,11** |
| fqf | 1,87 | 1,57 | 1,79 | **1,04** |
| fqf, n = 40 | 2,28 | 1,90 | 1,83 | 1,07 |
| fqf, densidad 1,0 | 2,69 | 2,16 | 1,55 | 1,05 |

Las tablas completas, por tipo de QPLIB, n, densidad y diagonal, están en [`resumen.txt`](resumen.txt), y los
datos en [`qplib_raiz.csv`](qplib_raiz.csv) y [`fqf_raiz.csv`](fqf_raiz.csv).

**Tiempos medianos en QPLIB**, en Python, salvo el k-nomial que es C++:

| k-nomial | PRODK | SDP (SCS) |
|---|---|---|
| 0,02 s | 2,8 s | 0,5 s (máximo 79 s con n = 175) |

## Lectura

1. **El término a término (Ibex) está muy lejos del exacto en formas cuadráticas generales:** 3,3 veces en QPLIB y
   1,9 en fqf. Hay mucho margen.
2. **Los métodos por bloques recuperan poco de ese margen.**
   - En QPLIB, el k-nomial y los bloques producto quitan alrededor de un 3 % del ancho (un 4 % con diagonal no
     nula, y nada con diagonal cero). Es lo mismo que había encontrado JCGeraldo, con 0,97 para k = 3.
   - En fqf quitan un 16 %, y hasta un 27 % con diagonal positiva, que es el caso favorable del k-nomial.
   - Los bloques producto le ganan al k-nomial solo con diagonal de ambos signos y Q dispersa. Su combinación
     (`BLOQ`) mejora poco al k-nomial solo.
3. **La SDP casi cierra la brecha:** queda a 1,04–1,11 veces el rango interior. Le gana a los bloques en las 58
   instancias de QPLIB y en todas las de fqf. Incluso la cota espectral, mucho más simple, es mucho mejor que los
   bloques en QPLIB (0,41 frente a 0,97).
4. **Diagonal cero:** es la situación de 21 de las 58 instancias de QPLIB. Ningún método por bloques hace nada, y
   solo ayudan la SDP y la espectral.
5. **Conclusión para RAA-K:** para formas `xᵀQx` generales, la herramienta adecuada es una cota tipo SDP de la Q
   completa, no la extracción de bloques de pocas variables. Los bloques serían, a lo más, un complemento barato
   cuando no se puede pagar una SDP.

## Idea para el solver (por medir)

La solución d de Shor en la raíz da `xᵀQx ≤ Σ dᵢxᵢ²` y `xᵀQx ≥ −Σ d'ᵢxᵢ²` en toda la caja raíz. Como esas cotas
son separables, se acotan exactas en cualquier subcaja, sin volver a resolver la SDP. Hay que medir cuánto se
pierde en cajas internas al reutilizar la d de la raíz, frente a recalcular la SDP, frente al término a término.

## Reproducir

```bash
# driver del k-nomial: ver ../raak-medicion-synq/README.md;  pip install cvxpy (usa SCS)
python3 gen_fqf.py instancias_fqf
python3 medir_q.py gen   ./knomial_driver fqf_raiz.csv   instancias_fqf/manifest.csv
python3 medir_q.py qplib ./knomial_driver qplib_raiz.csv <carpeta>/QPLIB_*.qplib
python3 resumir.py > resumen.txt
```
