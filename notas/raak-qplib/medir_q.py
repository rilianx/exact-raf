"""Medicion en la raiz de formas cuadraticas generales x^T Q x (QPLIB y generador ad hoc).

En la caja [l, u], x = m + D eps con D = diag(radios), asi que la parte cuadratica en ruido es eps^T (D Q D) eps.
Se mide el ancho del encierro de eps^T Qe eps en [-1,1]^n con:

  Qk0         termino a termino = lo que da Ibex (fAF2) con la Q escrita en monomios  (LINEA BASE)
  Qk2, Qk3, Qk6  k-nomial SDSOS de ibex_KNomial.h (raak-experimentos)
  PRODK       reduccion avida por bloques producto (prodk_fast.py)
  EIG         cota espectral n [min(lmin,0), max(lmax,0)]
  SDP         relajacion de Shor: max <= min sum d s.t. diag(d) - Q >= 0 (y lo mismo para el min)
  INNER       rango interior (busqueda local desde vertices), o EXACT por enumeracion de caras si n <= 12

Uso:
  python3 medir_q.py qplib <driver> <salida.csv> <archivos .qplib>...
  python3 medir_q.py gen   <driver> <salida.csv> <manifest.csv del generador>
"""
import csv, itertools, math, os, subprocess, sys, time
from pathlib import Path
import numpy as np
import cvxpy as cp

sys.path.insert(0, str(Path(__file__).resolve().parent))
import prodk_fast as PF  # noqa: E402

KS = (0, 2, 3, 6)


# ------------------------------------------------------------------ QPLIB
def read_qplib(path):
    """Devuelve (nombre, tipo, Q simetrica con f = x^T Q x + ..., lb, ub). Solo objetivo y cotas."""
    toks = []
    for line in Path(path).read_text().splitlines():
        line = line.split("#")[0].strip()
        if line: toks.append(line.split())
    it = iter(toks)
    name = next(it)[0]; typ = next(it)[0]; next(it)  # sentido
    n = int(next(it)[0])
    cons = typ[2] not in ("B", "N")
    m = int(next(it)[0]) if cons else 0
    Q = np.zeros((n, n))
    for _ in range(int(next(it)[0])):
        i, j, v = next(it); i, j, v = int(i) - 1, int(j) - 1, float(v)
        # QPLIB: f = 1/2 x^T Q0 x con el triangulo inferior de Q0
        Q[i, j] += 0.5 * v
        if i != j: Q[j, i] += 0.5 * v
    next(it)                                   # b por omision
    for _ in range(int(next(it)[0])): next(it)
    next(it)                                   # constante
    if cons:
        if typ[2] in ("Q", "C"):               # terminos cuadraticos de las restricciones
            for _ in range(int(next(it)[0])): next(it)
        for _ in range(int(next(it)[0])): next(it)   # terminos lineales
    inf = float(next(it)[0])
    if cons:
        next(it)
        for _ in range(int(next(it)[0])): next(it)
        next(it)
        for _ in range(int(next(it)[0])): next(it)
    lb, ub = np.zeros(n), np.ones(n)
    if typ[1] != "B":                          # variables no todas binarias: cotas explicitas
        d = float(next(it)[0]); lb[:] = d
        for _ in range(int(next(it)[0])):
            i, v = next(it); lb[int(i) - 1] = float(v)
        d = float(next(it)[0]); ub[:] = d
        for _ in range(int(next(it)[0])):
            i, v = next(it); ub[int(i) - 1] = float(v)
        if typ[1] in ("M", "I", "G"): raise ValueError("tipos de variable mixtos no soportados")
    lb[lb <= -inf] = -np.inf; ub[ub >= inf] = np.inf
    return name, typ, Q, lb, ub


# ------------------------------------------------------------------ cotas
def run_driver(driver, Qe):
    n = Qe.shape[0]
    lines = [f"Q {n} {k} " + " ".join(f"{x:.17g}" for x in Qe.ravel()) for k in KS]
    env = dict(os.environ, RAAK_NODE_BUDGET=os.environ.get("RAAK_NODE_BUDGET", "200000"))
    out = subprocess.run([driver], input="\n".join(lines) + "\n", capture_output=True, text=True, check=True, env=env)
    return [tuple(map(float, l.split())) for l in out.stdout.strip().splitlines()]


def sdp_width(Q, eps=1e-5):
    """Shor: max x^T Q x en [-1,1]^n <= sum d si d >= 0 y diag(d) - Q >= 0 (x^T Q x <= sum d_i x_i^2 <= sum d_i).
    Se resuelve con SCS (aproximado) y se repara: d = max(d, 0) y d += max(0, -lambda_min(diag(d) - Q)), lo que
    da una d factible y por tanto una cota valida."""
    n = Q.shape[0]
    res = []
    for s in (1, -1):                          # max de s * x^T Q x
        d = cp.Variable(n)
        prob = cp.Problem(cp.Minimize(cp.sum(d)), [d >= 0, cp.diag(d) - s * Q >> 0])
        prob.solve(solver=cp.SCS, eps=eps, max_iters=20000)
        dv = np.asarray(d.value, dtype=float) if d.value is not None else np.abs(s * Q).sum(axis=1)
        dv = np.maximum(dv, 0.0)
        lmin = np.linalg.eigvalsh(np.diag(dv) - s * Q)[0]
        dv = dv + max(0.0, -lmin) * (1 + 1e-12) + 1e-14 * np.abs(dv).max()
        res.append(dv.sum())
    return res[0] + res[1]                     # max - min <= sdp(Q) + sdp(-Q)


def quad_range_exact(Q):
    n = Q.shape[0]
    lo, hi = math.inf, -math.inf
    for mask in range(1 << n):
        F = [i for i in range(n) if mask >> i & 1]
        B = [i for i in range(n) if not mask >> i & 1]
        XB = np.array(list(itertools.product((-1.0, 1.0), repeat=len(B)))).T if B else np.zeros((0, 1))
        X = np.zeros((n, XB.shape[1])); X[B] = XB
        if F:
            QFF, QFB = Q[np.ix_(F, F)], Q[np.ix_(F, B)]
            rhs = -QFB @ XB if B else np.zeros((len(F), 1))
            XF = np.linalg.lstsq(QFF, rhs, rcond=None)[0]
            sc = max(1.0, np.abs(Q).max())
            ok = (np.abs(QFF @ XF - rhs).max(axis=0) <= 1e-9 * sc) & (np.abs(XF).max(axis=0) <= 1 + 1e-12)
            X[F] = XF; X = X[:, ok]
            if X.shape[1] == 0: continue
        vals = np.einsum("ij,ik,kj->j", X, Q, X)
        lo, hi = min(lo, vals.min()), max(hi, vals.max())
    return lo, hi


def quad_range_inner(Q, starts=40, seed=0):
    """Rango interior: ascenso por coordenadas (optimo exacto en cada coordenada) desde vertices aleatorios."""
    rng = np.random.default_rng(seed)
    n = Q.shape[0]
    dq = np.diag(Q).copy()
    out = []
    for sgn in (1.0, -1.0):
        best = -math.inf
        for _ in range(starts):
            x = rng.choice([-1.0, 1.0], size=n)
            g = Q @ x
            for _ in range(100):
                changed = False
                for i in range(n):
                    gi = g[i] - dq[i] * x[i]   # f(t) = dq t^2 + 2 gi t + const
                    cand = [-1.0, 1.0]
                    if sgn * dq[i] < 0: cand.append(min(1.0, max(-1.0, -gi / dq[i])))
                    t = max(cand, key=lambda t: sgn * (dq[i] * t * t + 2 * gi * t))
                    if t != x[i]:
                        g += Q[:, i] * (t - x[i]); x[i] = t; changed = True
                if not changed: break
            best = max(best, sgn * (x @ Q @ x))
        out.append(sgn * best)
    return out[1], out[0]


def measure(Qe, driver, exact_max=12):
    n = Qe.shape[0]
    r = {"n": n, "nnz_diag": int((np.diag(Qe) != 0).sum()),
         "dens_off": float((np.abs(Qe) > 0).sum() - (np.diag(Qe) != 0).sum()) / max(1, n * (n - 1))}
    lam = np.linalg.eigvalsh(Qe)
    sc = max(1e-300, np.abs(lam).max()) * 1e-9
    r["n_pos"], r["n_neg"] = int((lam > sc).sum()), int((lam < -sc).sum())
    t0 = time.time()
    for k, (L, U, nb, kmax) in zip(KS, run_driver(driver, Qe)): r[f"Qk{k}"] = U - L
    r["t_knomial"] = time.time() - t0
    t0 = time.time(); r["PRODK"], r["n_PRODK"] = PF.prodk(Qe); r["t_prodk"] = time.time() - t0
    r["EIG"] = n * (max(lam[-1], 0) - min(lam[0], 0))
    t0 = time.time(); r["SDP"] = sdp_width(Qe); r["t_sdp"] = time.time() - t0
    if n <= exact_max:
        lo, hi = quad_range_exact(Qe); r["exacto"] = True
    else:
        lo, hi = quad_range_inner(Qe); r["exacto"] = False
    r["INNER"] = hi - lo
    for k in ("Qk0", "Qk2", "Qk3", "Qk6", "PRODK", "EIG", "SDP"):
        assert r[k] >= r["INNER"] * (1 - 1e-7), (k, r[k], r["INNER"])
    return r


def write(rows, out_csv):
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys); w.writeheader(); w.writerows(rows)


def main():
    mode, driver, out_csv, rest = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
    rows = []
    if mode == "qplib":
        for f in rest:
            name, typ, Q, lb, ub = read_qplib(f)
            if not (np.isfinite(lb).all() and np.isfinite(ub).all()):
                print(f"{name}: omitida (cotas infinitas)", flush=True); continue
            D = (ub - lb) / 2
            Qe = D[:, None] * Q * D[None, :]
            if not np.any(Qe):
                print(f"{name}: omitida (Q nula)", flush=True); continue
            r = {"instancia": name, "tipo": typ, **measure(Qe, driver)}
            rows.append(r)
            print(f"{name} {typ} n={r['n']:3d} Qk6/Qk0={r['Qk6']/r['Qk0']:.3f} PRODK/Qk0={r['PRODK']/r['Qk0']:.3f} "
                  f"SDP/Qk0={r['SDP']/r['Qk0']:.3f} INNER/Qk0={r['INNER']/r['Qk0']:.3f}  t_prodk={r['t_prodk']:.1f}s", flush=True)
            write(rows, out_csv)
    else:
        import gen_fqf
        for g in csv.DictReader(open(rest[0])):
            Q = gen_fqf.load(Path(rest[0]).parent / g["archivo"])
            r = {**g, **measure(Q, driver)}
            rows.append(r)
            print(f"{g['instancia']:28s} Qk6/Qk0={r['Qk6']/r['Qk0']:.3f} PRODK/Qk0={r['PRODK']/r['Qk0']:.3f} "
                  f"SDP/Qk0={r['SDP']/r['Qk0']:.3f} INNER/Qk0={r['INNER']/r['Qk0']:.3f}", flush=True)
            write(rows, out_csv)


if __name__ == "__main__":
    main()
