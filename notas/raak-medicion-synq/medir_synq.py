"""Medicion de la seccion 4 de notas/raak-cambio-de-enfoque.md.

En la caja raiz [-1,1]^n de cada instancia synq, el objetivo es
    f(x) = sum_k (a_k^T x)(b_k^T x) + g^T x,
asi que en coordenadas de ruido (eps = x) la parte lineal es g^T x y el residuo cuadratico es x^T Q x con
Q = sum_k (a_k b_k^T + b_k a_k^T)/2. Se mide el ancho del encierro de x^T Q x con:

  S      suma por producto de la formula de fAF2 (S): [sum a_i b_i /2 -/+ (Sa Sb - sum|a_i b_i|/2)]
  SH     suma por producto del rango exacto de (a^T x)(b^T x) (= Skalna-Hladik, y = E: centros 0, sin errores)
  RAAp0  suma por producto de knomial_quad_bound con k = 0 (RAA termino a termino)
  RAAp6  suma por producto de knomial_quad_bound con k = 6
  Qk0    Q completa, termino a termino (k = 0)
  Qk2, Qk4, Qk6, Qk8  Q completa con la reescritura k-nomial SDSOS de ibex_KNomial.h
  EIG    cota espectral: x^T Q x en [n min(lmin,0), n max(lmax,0)]
  EXACT  rango exacto de x^T Q x en la caja (enumeracion de caras KKT)

Uso: python3 medir_synq.py <dir_instancias> <driver> <salida.csv>
"""
import csv, itertools, re, subprocess, sys
from pathlib import Path
import numpy as np

KS_FULL = (0, 2, 4, 6, 8)


def parse_lin(s, n):
    v = np.zeros(n)
    s = s.replace(" ", "")
    for coef, i in re.findall(r"([+-]?[0-9.]*)\*?x(\d+)", s):
        c = 1.0 if coef in ("", "+") else (-1.0 if coef == "-" else float(coef))
        v[int(i) - 1] += c
    return v


def parse_objective(path):
    txt = Path(path).read_text()
    n = len(re.findall(r"^x\d+ in", txt, re.M))
    obj = re.search(r"minimize\s*(.*?);", txt, re.S).group(1).strip()
    terms, depth, cur = [], 0, ""
    for ch in obj:  # separa por '+' de nivel 0
        depth += ch == "("
        depth -= ch == ")"
        if ch == "+" and depth == 0:
            terms.append(cur.strip()); cur = ""
        else:
            cur += ch
    terms.append(cur.strip())
    prods, g = [], np.zeros(n)
    for t in terms:
        parts = re.findall(r"\(([^()]*)\)", t)
        if len(parts) == 2:
            prods.append((parse_lin(parts[0], n), parse_lin(parts[1], n)))
        elif len(parts) == 1:
            g += parse_lin(parts[0], n)
        else:
            raise ValueError(f"termino no reconocido en {path}: {t}")
    return n, prods, g


def hull(pts):
    pts = sorted(set(map(tuple, pts)))
    if len(pts) <= 2:
        return pts
    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])
    lo, up = [], []
    for p in pts:
        while len(lo) >= 2 and cross(lo[-2], lo[-1], p) <= 0: lo.pop()
        lo.append(p)
    for p in reversed(pts):
        while len(up) >= 2 and cross(up[-2], up[-1], p) <= 0: up.pop()
        up.append(p)
    return lo[:-1] + up[:-1]


def product_range(a, b):
    """Rango exacto de (a^T x)(b^T x) en [-1,1]^n: u v sobre el zonotopo plano (u, v)."""
    n = len(a)
    V = np.array(list(itertools.product((-1.0, 1.0), repeat=n)))
    H = hull(np.column_stack([V @ a, V @ b]))
    lo, hi = np.inf, -np.inf
    for i in range(len(H)):
        p, q = np.array(H[i]), np.array(H[(i + 1) % len(H)])
        d = q - p
        ts = [0.0, 1.0]
        den = 2 * d[0] * d[1]
        if den != 0:
            t = -(p[0] * d[1] + p[1] * d[0]) / den
            if 0 < t < 1: ts.append(t)
        for t in ts:
            u, v = p + t * d
            lo, hi = min(lo, u * v), max(hi, u * v)
    return lo, hi


def quad_range_exact(Q):
    """min y max de x^T Q x en [-1,1]^n: todo optimo global es punto KKT en el interior relativo de
    alguna cara; en la cara con libres F y fijas B = +-1 se resuelve Q_FF x_F = -Q_FB x_B."""
    n = Q.shape[0]
    lo, hi = np.inf, -np.inf
    for mask in range(1 << n):
        F = [i for i in range(n) if mask >> i & 1]
        B = [i for i in range(n) if not mask >> i & 1]
        XB = np.array(list(itertools.product((-1.0, 1.0), repeat=len(B)))).T if B else np.zeros((0, 1))
        X = np.zeros((n, XB.shape[1]))
        X[B] = XB
        if F:
            QFF, QFB = Q[np.ix_(F, F)], Q[np.ix_(F, B)]
            rhs = -QFB @ XB if B else np.zeros((len(F), 1))
            XF = np.linalg.lstsq(QFF, rhs, rcond=None)[0]
            ok = (np.abs(QFF @ XF - rhs).max(axis=0) <= 1e-9) & (np.abs(XF).max(axis=0) <= 1 + 1e-12)
            X[F] = XF
            X = X[:, ok]
            if X.shape[1] == 0:
                continue
        vals = np.einsum("ij,ik,kj->j", X, Q, X)
        lo, hi = min(lo, vals.min()), max(hi, vals.max())
    return lo, hi


def run_driver(driver, lines):
    out = subprocess.run([driver], input="\n".join(lines) + "\n", capture_output=True, text=True, check=True)
    return [tuple(map(float, l.split())) for l in out.stdout.strip().splitlines()]


def main():
    inst_dir, driver, out_csv = sys.argv[1:4]
    files = sorted(Path(inst_dir).glob("synq_*.bch"))
    manifest = {r["instance"]: r for r in csv.DictReader(open(Path(inst_dir) / "manifest.csv"))}
    rows = []
    for f in files:
        n, prods, g = parse_objective(f)
        Q = sum((np.outer(a, b) + np.outer(b, a)) / 2 for a, b in prods)
        Qm = np.array(eval(manifest[f.name]["Q_efectiva"]))  # comprobacion del parser
        assert np.abs(Q - Qm).max() < 1e-3, f.name
        r = {"instancia": f.stem, "n": n, "K": len(prods),
             "densidad": manifest[f.name]["densidad"], "diagonal": manifest[f.name]["diagonal"],
             "convexa": manifest[f.name]["convexa"], "rango_Q": int(np.linalg.matrix_rank(Q, tol=1e-9))}
        # S (fAF2) y SH por producto
        wS = wSH = 0.0
        for a, b in prods:
            wS += 2 * (np.abs(a).sum() * np.abs(b).sum() - 0.5 * np.abs(a * b).sum())
            lo, hi = product_range(a, b)
            wSH += hi - lo
        r["S"], r["SH"] = wS, wSH
        # driver: por producto k = 0, 6 y Q completa
        lines = []
        for k in (0, 6):
            lines += [f"P {n} {k} " + " ".join(f"{v:.17g}" for v in np.r_[a, b]) for a, b in prods]
        lines += [f"Q {n} {k} " + " ".join(f"{v:.17g}" for v in Q.ravel()) for k in KS_FULL]
        res = run_driver(driver, lines)
        K = len(prods)
        r["RAAp0"] = sum(U - L for L, U, *_ in res[:K])
        r["RAAp6"] = sum(U - L for L, U, *_ in res[K:2 * K])
        for k, (L, U, nb, kmax) in zip(KS_FULL, res[2 * K:]):
            r[f"Qk{k}"] = U - L
            r[f"Qk{k}_bloques"] = int(nb)
            r[f"Qk{k}_tammax"] = int(kmax)
        ev = np.linalg.eigvalsh(Q)
        r["EIG"] = n * (max(ev[-1], 0) - min(ev[0], 0))
        lo, hi = quad_range_exact(Q)
        r["EXACT_lo"], r["EXACT_hi"], r["EXACT"] = lo, hi, hi - lo
        # comprobacion: todas las cotas contienen el rango exacto (anchos >= exacto)
        for key in ["S", "SH", "RAAp0", "RAAp6", "EIG"] + [f"Qk{k}" for k in KS_FULL]:
            assert r[key] >= r["EXACT"] - 1e-9, (f.name, key, r[key], r["EXACT"])
        rows.append(r)
        print(f"{f.stem:34s} SH={wSH:8.3f} Qk0={r['Qk0']:8.3f} Qk6={r['Qk6']:8.3f} EXACT={r['EXACT']:8.3f}", flush=True)
    with open(out_csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)


if __name__ == "__main__":
    main()
