"""Reduccion avida de Q por bloques PRODUCTO (u^T z)(v^T z): la generalizacion del k-nomial (que solo extrae
cuadrados c (u^T z)^2) a productos de dos sumas, encerrados exactos con SH.

    z^T Q z = sum_p (u_p^T z)(v_p^T z) + z^T R z,    ancho <= sum_p SH(u_p, v_p) + TT(R)

con TT el encierro termino a termino del residuo (el mismo de knomial_quad_bound: R_ii [0,1] + 2 R_ij [-1,1]).
La cota es valida para cualquier eleccion de bloques; la eleccion solo decide el ancho.

En cada paso, para cada soporte S (pares y triples conectados de Q; el soporte completo si m <= 12) se generan
candidatos P sobre R_S:
  - par espectral de R_S: (a+b)(a-b) con el mayor autovalor positivo y el menor negativo;
  - el autovector de mayor |lambda| como cuadrado;
  - filas de factor comun z_i (R_ii z_i + 2 sum_j R_ij z_j) restringidas a S;
cada uno escalado por t en {1, 3/4, 1/2, 1/4} (el analogo del lambda del k-nomial). Se extrae el de mayor
reduccion de ancho y se repite hasta que ninguno reduzca.

Uso: python3 medir_prodk.py synq   <dir_synq>  <synq_pares.csv> <salida.csv>
     python3 medir_prodk.py reales <reales_pares.csv> <salida.csv> <bch>...
"""
import csv, itertools, sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent / "raak-medicion-reales"))
import medir_reales as MR  # noqa: E402
import medir_pares as MP   # noqa: E402

TS = (1.0, 0.75, 0.5, 0.25)


def tt(R):
    return np.abs(np.diag(R)).sum() + 2 * (np.abs(R).sum() - np.abs(np.diag(R)).sum())


def sh(pairs):
    return sum(hi - lo for lo, hi in (MR.product_range(u, v) for u, v in pairs))


def as_matrix(pairs, m):
    P = np.zeros((m, m))
    for u, v in pairs: P += (np.outer(u, v) + np.outer(v, u)) / 2
    return P


def candidates(R, S, m):
    S = list(S)
    RS = R[np.ix_(S, S)]
    out = []
    lam, V = np.linalg.eigh(RS)
    sc = max(1e-300, np.abs(lam).max()) * 1e-12
    def emb(x):
        y = np.zeros(m); y[S] = x; return y
    ip, ineg = int(np.argmax(lam)), int(np.argmin(lam))
    if lam[ip] > sc and lam[ineg] < -sc:
        a, b = np.sqrt(lam[ip]) * V[:, ip], np.sqrt(-lam[ineg]) * V[:, ineg]
        out.append([(emb(a + b), emb(a - b))])
    k = int(np.argmax(np.abs(lam)))
    if abs(lam[k]) > sc:
        a = np.sqrt(abs(lam[k])) * V[:, k]
        out.append([(emb(a), emb(np.sign(lam[k]) * a))])
    for t, i in enumerate(S):
        v = 2 * RS[t].copy(); v[t] = RS[t, t]
        if np.abs(v).sum() > 0:
            u = np.zeros(m); u[i] = 1.0
            out.append([(u, emb(v))])
    return out


def supports(Q, full_max=12):
    m = Q.shape[0]
    nz = np.abs(Q) > 0
    S = set()
    for i, j in itertools.combinations(range(m), 2):
        if nz[i, j] or (nz[i, i] and nz[j, j]): S.add((i, j))
    for i, j, k in itertools.combinations(range(m), 3):
        e = nz[i, j] + nz[i, k] + nz[j, k]
        if e >= 2: S.add((i, j, k))
    if m <= full_max: S.add(tuple(range(m)))
    return sorted(S)


def prodk(Q, max_iter=200):
    m = Q.shape[0]
    R = Q.copy()
    sup = supports(Q)
    blocks = []
    cur = tt(R)
    for _ in range(max_iter):
        best, best_val = None, cur * (1 - 1e-9)
        for S in sup:
            for pairs in candidates(R, S, m):
                for t in TS:
                    ps = [(t * u, v) for u, v in pairs]
                    val = sh(blocks + ps) - sh(blocks) + tt(R - as_matrix(ps, m)) + sh(blocks)
                    if val < best_val: best, best_val = ps, val
        if best is None: break
        blocks += best
        R = R - as_matrix(best, m)
        cur = sh(blocks) + tt(R)
    return cur, len(blocks)


def main_synq(inst_dir, prev_csv, out_csv):
    sys.path.insert(0, str(HERE.parent / "raak-medicion-synq"))
    import medir_synq as MS
    prev = {r["instancia"]: r for r in csv.DictReader(open(prev_csv))}
    rows = []
    for f in sorted(Path(inst_dir).glob("synq_*.bch")):
        n, prods, g = MS.parse_objective(f)
        Q = sum((np.outer(a, b) + np.outer(b, a)) / 2 for a, b in prods)
        w, nb = prodk(Q)
        p = prev[f.stem]
        assert w >= float(p["EXACT"]) * (1 - 1e-9), f.stem
        rows.append({**p, "PRODK": w, "n_PRODK": nb})
        print(f"{f.stem:34s} PRODK/EXACT={w / float(p['EXACT']):.3f} Qk6/EXACT={float(p['Qk6']) / float(p['EXACT']):.3f}", flush=True)
    MP.write(rows, out_csv)


def main_reales(prev_csv, out_csv, files):
    prev = list(csv.DictReader(open(prev_csv)))
    rows, idx = [], 0
    for f in files:
        box, consts, trees = MR.read_bch(f)
        ev = MR.Ev(box, consts)
        for name, tree in trees:
            for terms in ev.sum_nodes(tree, []):
                try:
                    U, V = MP.node_uv(ev, terms)
                except (ValueError, OverflowError, ZeroDivisionError):
                    continue
                if U is None: continue
                p = prev[idx]; idx += 1
                assert p["instancia"] == Path(f).stem and p["expr"] == name
                Q = sum((np.outer(u, v) + np.outer(v, u)) / 2 for u, v in zip(U, V))
                w, nb = prodk(Q)
                if p["exacto"] == "True": assert w >= float(p["EXACT"]) * (1 - 1e-9), (f, name)
                rows.append({**p, "PRODK": w, "n_PRODK": nb})
        print(f"{Path(f).stem:20s} listo", flush=True)
    assert idx == len(prev)
    MP.write(rows, out_csv)


if __name__ == "__main__":
    if sys.argv[1] == "synq": main_synq(*sys.argv[2:5])
    else: main_reales(sys.argv[2], sys.argv[3], sys.argv[4:])
