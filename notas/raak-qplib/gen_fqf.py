"""Generador ad hoc de formas cuadraticas generales f(x) = x^T Q x + c^T x en [-1,1]^n (familia 'fqf').

A diferencia de synq (Q = A^T M factorizada a proposito, con cancelacion artificial entre productos), aqui Q se
genera directamente y se escribe en MONOMIOS, que es como aparece en QPLIB y en los modelos reales. La linea base
de Ibex es entonces el encierro termino a termino (Qk0).

Factores (los que, segun las mediciones anteriores, deciden que metodo gana):
  n          10, 20, 40
  densidad   fraccion de entradas fuera de la diagonal no nulas: 0.1, 0.3, 1.0 (valores U[-1, 1])
  diagonal   'cero'  (el k-nomial no puede formar bloques; como QPLIB QBL),
             'pos'   (d_i = U[0.2, 1] * sum_j |Q_ij|, del mismo signo: material para cuadrados),
             'mixta' (como 'pos' pero con signo aleatorio)
  semillas   3
Total 3 x 3 x 3 x 3 = 81 instancias.

Salida: <out>/fqf_n<n>_d<dens>_<diag>_s<seed>.txt (Q por filas), .bch (Minibex, en monomios) y manifest.csv.
Uso: python3 gen_fqf.py <out>
"""
import csv, sys
from pathlib import Path
import numpy as np

NS, DENS, DIAGS, SEEDS = (10, 20, 40), (0.1, 0.3, 1.0), ("cero", "pos", "mixta"), (1, 2, 3)


def gen(n, dens, diag, seed):
    rng = np.random.default_rng(1000003 * n + 1009 * int(dens * 100) + 17 * DIAGS.index(diag) + seed)
    Q = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            if rng.random() < dens:
                Q[i, j] = Q[j, i] = round(rng.uniform(-1, 1), 3)
    if diag != "cero":
        rows = np.abs(Q).sum(axis=1)
        d = rng.uniform(0.2, 1.0, size=n) * np.maximum(rows, 0.5)
        if diag == "mixta": d *= rng.choice([-1.0, 1.0], size=n)
        Q[np.diag_indices(n)] = np.round(d, 3)
    c = np.round(rng.uniform(-1, 1, size=n), 3)
    return Q, c


def write_bch(path, Q, c):
    n = Q.shape[0]
    terms = []
    for i in range(n):
        if Q[i, i]: terms.append(f"{Q[i, i]:+.3f}*x{i + 1}^2")
        for j in range(i + 1, n):
            if Q[i, j]: terms.append(f"{2 * Q[i, j]:+.3f}*x{i + 1}*x{j + 1}")
    terms += [f"{c[i]:+.3f}*x{i + 1}" for i in range(n) if c[i]]
    var = "\n".join(f"x{i + 1} in [-1, 1];" for i in range(n))
    obj = " ".join(terms).lstrip("+")          # Minibex no acepta un '+' unario al comienzo
    Path(path).write_text(f"variables\n{var}\n\nminimize\n" + obj + ";\n\nconstraints\n\nend\n")   # ibexopt exige la seccion


def load(path):
    return np.loadtxt(path, ndmin=2)


def main(out):
    out = Path(out); out.mkdir(parents=True, exist_ok=True)
    rows = []
    for n in NS:
        for dens in DENS:
            for diag in DIAGS:
                for s in SEEDS:
                    Q, c = gen(n, dens, diag, s)
                    name = f"fqf_n{n:02d}_d{int(dens * 100):03d}_{diag}_s{s}"
                    np.savetxt(out / f"{name}.txt", Q, fmt="%.3f")
                    write_bch(out / f"{name}.bch", Q, c)
                    lam = np.linalg.eigvalsh(Q)
                    rows.append({"instancia": name, "archivo": f"{name}.txt", "n_gen": n, "densidad": dens,
                                 "diagonal": diag, "semilla": s, "lambda_min": round(lam[0], 4), "lambda_max": round(lam[-1], 4)})
    with open(out / "manifest.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    print(f"{len(rows)} instancias en {out}")


if __name__ == "__main__":
    main(sys.argv[1])
