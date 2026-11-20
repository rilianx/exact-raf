"""Recalcula solo la columna PRODK de qplib_raiz.csv y fqf_raiz.csv con prodk_fast corregido (orden angular
exacto en sh_range; la version con atan2 podia dar rangos SH demasiado angostos). Informa el cambio.
Uso: python3 recalcular_prodk.py <carpeta con QPLIB_*.qplib>"""
import csv, sys
from pathlib import Path
import numpy as np
import medir_q as M, prodk_fast as PF, gen_fqf
qdir = Path(sys.argv[1])
for csvf, load in (("qplib_raiz.csv", lambda r: M.read_qplib(qdir / f"{r['instancia']}.qplib")),
                   ("fqf_raiz.csv", lambda r: gen_fqf.load(Path("instancias_fqf") / r["archivo"]))):
    rows = list(csv.DictReader(open(csvf)))
    ch = []
    for r in rows:
        x = load(r)
        if isinstance(x, tuple):
            _, _, Q, lb, ub = x; D = (ub - lb) / 2; Q = D[:, None] * Q * D[None, :]
        else: Q = x
        w, nb = PF.prodk(Q)
        ch.append(w / float(r["PRODK"]))
        assert w >= float(r["INNER"]) * (1 - 1e-7), (r["instancia"], w, r["INNER"])
        r["PRODK"], r["n_PRODK"] = w, nb
    with open(csvf, "w", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
    print(f"{csvf}: PRODK nuevo/viejo min {min(ch):.4f} max {max(ch):.4f}; cambian > 0,1 %: {sum(abs(c - 1) > 1e-3 for c in ch)} de {len(ch)}")
