"""Reescritura de una suma de productos como pocos productos de PARES DE SUMAS, para que SH/E los encierren.

Para la misma Q de las mediciones anteriores (synq: objetivo; reales: cada nodo suma con >= 2 terminos
cuadraticos) se comparan descomposiciones Q = sum_p (u_p v_p^T + v_p u_p^T)/2, cada producto (u_p^T z)(v_p^T z)
encerrado con su rango exacto (SH) y los anchos sumados:

  SH     los productos originales de la expresion (linea base)
  STAR   factor comun: z_i * (Q_ii z_i + 2 sum_j Q_ij z_j), eligiendo con avidez la variable de mayor peso
         restante (lo que haria una factorizacion simbolica x*y + x*z -> x*(y+z))
  EIGP   pares espectrales: Q = sum l_k q_k q_k^T; el i-esimo autovalor positivo con el i-esimo negativo
         como diferencia de cuadrados (a+b)(a-b); los que sobran, como cuadrados a*a
  EIGSQ  todos los autovectores como cuadrados, sin emparejar

Uso:
  python3 medir_pares.py synq   <dir_synq>  <synq_raiz.csv>  <salida.csv>
  python3 medir_pares.py reales <nodos.csv> <salida.csv> <bch>...
"""
import csv, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "raak-medicion-reales"))
import medir_reales as MR  # noqa: E402  (product_range por zonotopo, lector de Minibex, evaluador afin)


def width(pairs):
    return sum(hi - lo for lo, hi in (MR.product_range(u, v) for u, v in pairs))


def star_pairs(Q):
    R = Q.copy()
    m = R.shape[0]
    pairs = []
    off = lambda: np.abs(R - np.diag(np.diag(R)))
    while off().max() > 0:
        i = int(np.argmax(off().sum(axis=1) + np.abs(np.diag(R))))
        v = 2 * R[i].copy(); v[i] = R[i, i]
        u = np.zeros(m); u[i] = 1.0
        pairs.append((u, v))
        R[i, :] = 0; R[:, i] = 0
    for j in np.nonzero(np.diag(R))[0]:
        u = np.zeros(m); u[j] = 1.0
        pairs.append((u, R[j, j] * u))
    return pairs


def eig_pairs(Q, pair=True, tol=1e-12):
    lam, P = np.linalg.eigh(Q)
    sc = max(1.0, np.abs(lam).max()) * tol
    pos = [(l, P[:, k]) for k, l in sorted(enumerate(lam), key=lambda t: -t[1]) if l > sc]
    neg = [(l, P[:, k]) for k, l in sorted(enumerate(lam), key=lambda t: t[1]) if l < -sc]
    pairs = []
    if pair:
        for (lp, p), (ln, q) in zip(pos, neg):
            a, b = np.sqrt(lp) * p, np.sqrt(-ln) * q
            pairs.append((a + b, a - b))          # (a+b)(a-b) = lp (p^T z)^2 - |ln| (q^T z)^2
        k = min(len(pos), len(neg))
        pos, neg = pos[k:], neg[k:]
    for l, p in pos + neg:
        a = np.sqrt(abs(l)) * p
        pairs.append((a, np.sign(l) * a))
    return pairs


def measure(Q, orig_pairs):
    r = {"SH": width(orig_pairs), "STAR": width(star_pairs(Q)),
         "EIGP": width(eig_pairs(Q)), "EIGSQ": width(eig_pairs(Q, pair=False))}
    r["n_STAR"], r["n_EIGP"] = len(star_pairs(Q)), len(eig_pairs(Q))
    return r


def main_synq(inst_dir, prev_csv, out_csv):
    sys.path.insert(0, str(HERE.parent / "raak-medicion-synq"))
    import medir_synq as MS
    prev = {r["instancia"]: r for r in csv.DictReader(open(prev_csv))}
    rows = []
    for f in sorted(Path(inst_dir).glob("synq_*.bch")):
        n, prods, g = MS.parse_objective(f)
        Q = sum((np.outer(a, b) + np.outer(b, a)) / 2 for a, b in prods)
        r = {"instancia": f.stem, "n": n, **measure(Q, prods)}
        p = prev[f.stem]
        assert abs(r["SH"] - float(p["SH"])) <= 1e-9 * max(1, r["SH"]), f.stem
        for k in ("Qk0", "Qk6", "EIG", "EXACT"): r[k] = float(p[k])
        rows.append(r)
    write(rows, out_csv)


def main_reales(prev_csv, out_csv, files):
    prev = list(csv.DictReader(open(prev_csv)))
    rows, idx = [], 0
    for f in files:
        box, consts, trees = MR.read_bch(f)
        ev = MR.Ev(box, consts)
        for name, tree in trees:
            for terms in ev.sum_nodes(tree, []):
                try:
                    U, V = node_uv(ev, terms)
                except (ValueError, OverflowError, ZeroDivisionError):
                    continue
                if U is None: continue
                p = prev[idx]; idx += 1
                assert p["instancia"] == Path(f).stem and p["expr"] == name and int(p["m"]) == U.shape[1], (f, name)
                Q = sum((np.outer(u, v) + np.outer(v, u)) / 2 for u, v in zip(U, V))
                r = {"instancia": p["instancia"], "expr": name, "K": p["K"], "m": p["m"], **measure(Q, list(zip(U, V)))}
                assert abs(r["SH"] - float(p["SH"])) <= 1e-9 * max(1, r["SH"]), (f, name)
                for k in ("Qk0", "Qk6", "EIG", "EXACT", "exacto"): r[k] = p[k]
                rows.append(r)
    assert idx == len(prev)
    write(rows, out_csv)


def node_uv(ev, terms):
    """Los mismos vectores (u_k, v_k) que measure_node de medir_reales.py."""
    quads = [(c, *ev.quad_factors(t)) for c, t in terms if ev.quad_factors(t)]
    if len(quads) < 2: return None, None
    syms, vecs = {}, []
    for c, A, B in quads:
        uv = []
        for X in (A, B):
            fx = ev.af(X)
            d = dict(fx.lin)
            if fx.err > 0: d["E:" + MR.ast.unparse(X)] = fx.err
            uv.append(d)
            for k in d: syms.setdefault(k, len(syms))
        vecs.append((c, uv[0], uv[1]))
    U = np.zeros((len(vecs), len(syms))); V = np.zeros_like(U)
    for i, (c, du, dv) in enumerate(vecs):
        for k, val in du.items(): U[i, syms[k]] = c * val
        for k, val in dv.items(): V[i, syms[k]] = val
    return U, V


def write(rows, out_csv):
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    if sys.argv[1] == "synq": main_synq(*sys.argv[2:5])
    else: main_reales(sys.argv[2], sys.argv[3], sys.argv[4:])
