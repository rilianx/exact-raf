# Multiplicación exacta en Ibex (`AF_fAF2`): especificación

Documento para implementar en `ibex-affine` (cadena de parches `raak/patches`, sobre `c219009c`) la
multiplicación exacta del paper *Exact multiplication of reduced affine forms* (`paper/numa/exact_raf.tex`,
Sección 3) y la de Skalna y Hladík (2017, Ec. 10) como línea base. Con esto se cubren tres cosas que le
faltan al paper: implementación rigurosa, tiempos en código compilado y experimentos en el solver.

## 1. Qué se calcula

Formas de `AF_fAF2`: `x = x0 + Σ a_i ε_i + ex·ξ` (`_val[0]`, `_val[i]`, `_err`), igual para `y` con `b_i`, `ey`, `η`.

```
x·y = x0·y0 + Σ (x0 b_i + y0 a_i) ε_i + N
N(ε, ξ, η) = x0·ey·η + y0·ex·ξ + (aᵀε + ex ξ)(bᵀε + ey η)
```

La parte lineal es la que ya calcula `operator*=` (con su contabilidad de redondeo `ttt`/`sss`). Lo único
que cambia es el intervalo `[nlo, nhi]` que encierra `N`: el centro nuevo es `x0·y0 + (nlo+nhi)/2` y el
error nuevo `(nhi−nlo)/2` más los errores de redondeo que ya se acumulan.

Reglas a comparar (todas con la misma parte lineal):

| `--arith` | Encierro de `N` |
|---|---|
| `af2` | fAF2 actual = fórmula (5) de Skalna-Hladík |
| `raa` | RAA actual: cuadrática término a término + corner test de los cruces `e·L`; `\|y0\| ex + \|x0\| ey` aparte |
| `sh` (nuevo) | rango exacto de `R = (aᵀε + ex ξ)(bᵀε + ey η)` + `[−1,1]·(\|y0\| ex + \|x0\| ey)` |
| `exact` (nuevo) | rango exacto de `N` |

`raa_corner_cross` y `knomial_quad_bound` (k=0) quedan como están para `raa`; `sh` y `exact` no los usan.

## 2. Algoritmo (Algoritmo 1 del paper)

1. Generadores `g_i = (a_i, b_i)` no nulos; voltearlos (`g ← −g`) para que queden en `{u>0} ∪ {u=0, v>0}`;
   ordenarlos por ángulo creciente. **Una sola vez** por producto.
2. Polígono del zonotopo `Z`: `p_0 = −Σ g_k`, `p_k = p_{k−1} + 2 g_k` (k = 1..m), y luego `p_{m+k} = p_{m+k−1} − 2 g_k`.
3. Para cada esquina `(s, t) ∈ {−1,1}²`: desplazamiento `c = (ex·s, ey·t)` y constante `κ = x0·ey·t + y0·ex·s`.
   Recorrer las `2m` aristas `[p, q]` del polígono desplazado por `c`; en cada una evaluar
   `φ(u,v) = u·v` en `p`, en `q` y, si existe, en el punto estacionario de `τ ↦ φ(p + τ(q−p))`, `τ ∈ (0,1)`.
   Actualizar `nlo`, `nhi` con `κ + valor`.
4. `sh`: lo mismo, pero una sola vuelta sobre el zonotopo de `{g_i} ∪ {(ex,0), (0,ey)}` sin desplazamiento,
   y sumar `±(|y0| ex + |x0| ey)`.

Costo: `O(n log n)` (orden) + 4 recorridos `O(n)`. Si `m = 0` (sin generadores), el polígono es un punto.

Atajo (Corolario 11 del paper): `ancho(SH) − ancho(E) ≤ 2·min(|x0|·ey, |y0|·ex)`. Si ese mínimo es 0
(un factor sin error o con centro 0), `exact` = `sh` y basta una vuelta.

## 3. Rigor (lo delicado)

El encierro debe contener el rango exacto aunque el cálculo sea en punto flotante:

- **Orden por ángulo con predicado exacto.** Un orden equivocado entre generadores casi paralelos produce
  un camino que entra al interior de `Z` y puede **subestimar** el máximo. Comparar con el signo exacto de
  `u1·v2 − u2·v1` (producto sin error con `fma`/`twoProd`, como en los predicados de Shewchuk), no con `atan2`.
  El volteo `g ← −g` es exacto.
- **Vértices como intervalos.** `p_k` se acumula con aritmética de intervalos (`ibex::Interval`), igual que
  el desplazamiento `c`.
- **Valores en aristas.** Con `p`, `q` intervalos: `φ(p)`, `φ(q)` por producto de intervalos. Para el punto
  estacionario, escribir `φ(p + τ d) = P + B τ + A τ²` con `P = p_u p_v`, `B = p_u d_v + p_v d_u`,
  `A = d_u d_v` (intervalos). Si `A` no contiene 0 y `τ* = −B/(2A)` corta `(0,1)`, incluir
  `P − B²/(4A)` (en intervalos). Si `A` contiene 0, incluir la extensión de Horner `P + τ(B + Aτ)` con `τ = [0,1]`.
- **Resultado.** `nlo` = mínimo de los extremos inferiores, `nhi` = máximo de los superiores; centro y radio
  con redondeo hacia afuera y sumados a `_val[0]` y `_err` con los mismos `twoSum`/`(1+2·AF_EM)` que hoy.

Con esto el sobrecosto de rigor es del orden del redondeo, como en el parche 0012.

## 4. Validación

- `affine_containment_check.cpp`: agregar `sh` y `exact` (cero violaciones).
- Contenciones en un producto (deben valer siempre): `exact ⊆ sh ⊆ af2` y `exact ⊆ raa`.
- Coincidencias: `exact == sh` cuando `ex = 0`, `ey = 0`, `x0 = 0` o `y0 = 0`; `exact` exacto (igual al óptimo)
  cuando `ex = ey = 0`.
- Ejemplo del paper: `x = −1 + ξ`, `y = 1 + η` (sin `ε`): `exact` da `[−4, 0]`, `sh` da `[−4, 2]`.
- Comparar contra la referencia en Python de abajo en productos aleatorios (tolerancia relativa ~1e-12).

## 5. Experimentos

1. **Tiempo por producto** (C++), `n = 2, 5, 10, 20, 50, 100`, con y sin error en los factores:
   `af2`, `raa`, `sh`, `exact`. Es lo que el paper reporta hoy solo en Python.
2. **Solver**, repitiendo los barridos 22 (`ibexopt`, `art`, synq/synl/synt y benchmarks, 1800 s) y 20
   (`ibexsolve`, las 17 instancias con material y las 264 a 60 s), con `--arith=af2|raa|sh|exact` y k = 0.
   Métricas: celdas y CPU (media geométrica), resueltas, cota inferior en las no resueltas.
3. Opcional: contar en cuántos productos `min(|x0|ey, |y0|ex) > 0` (los únicos donde `exact` puede ganar a `sh`).

## 6. Referencia en Python (coordenadas de ruido, sin redondeo)

```python
import math

def uv_range(gens, c0=0.0, c1=0.0):
    """Rango exacto de (u+c0)(v+c1) sobre el zonotopo sum_k e_k g_k, e en [-1,1]^K."""
    g = [(-u, -v) if (u < 0 or (u == 0 and v < 0)) else (u, v) for u, v in gens if u or v]
    if not g:
        return c0 * c1, c0 * c1
    g.sort(key=lambda t: math.atan2(t[1], t[0]))
    pu, pv = c0 - sum(t[0] for t in g), c1 - sum(t[1] for t in g)
    lo = hi = pu * pv
    for sgn in (2.0, -2.0):
        for gu, gv in g:
            du, dv = sgn * gu, sgn * gv
            qu, qv = pu + du, pv + dv
            lo, hi = min(lo, qu * qv), max(hi, qu * qv)
            A, B = du * dv, pu * dv + pv * du
            if A != 0 and 0 < -B / (2 * A) < 1:
                t = -B / (2 * A)
                val = (pu + t * du) * (pv + t * dv)
                lo, hi = min(lo, val), max(hi, val)
            pu, pv = qu, qv
    return lo, hi

def mul_exact_noise(x0, a, ex, y0, b, ey):
    """Encierro [nlo, nhi] de N; el producto es x0*y0 + sum (x0 b_i + y0 a_i) e_i + [nlo, nhi]."""
    gens = list(zip(a, b))
    nlo, nhi = math.inf, -math.inf
    for s in (-1, 1):
        for t in (-1, 1):
            lo, hi = uv_range(gens, ex * s, ey * t)
            k = x0 * ey * t + y0 * ex * s
            nlo, nhi = min(nlo, lo + k), max(nhi, hi + k)
    return nlo, nhi

def mul_sh_noise(x0, a, ex, y0, b, ey):
    lo, hi = uv_range(list(zip(a, b)) + [(ex, 0.0), (0.0, ey)])
    e = abs(y0) * ex + abs(x0) * ey
    return lo - e, hi + e
```

La implementación completa en coordenadas de variables (con pruebas) está en `zcaa/zcaa.py`
(`mul_exact`, `mul_skalna`) del repo `rilianx/raa-affine-arithmetic`.
