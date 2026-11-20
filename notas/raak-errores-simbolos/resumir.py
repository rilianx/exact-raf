"""Resumen de cajas.csv: ancho de M1 y M2 relativo a M0 (fAF2 + SH hoy), y al rango interior, por profundidad.
Media geometrica: primero sobre las expresiones y cajas de cada instancia, despues sobre instancias."""
import csv, sys
from math import exp, log
gm = lambda xs: exp(sum(map(log, xs)) / len(xs)) if xs else float("nan")
rows = [r for r in csv.DictReader(open(sys.argv[1] if len(sys.argv) > 1 else "cajas.csv"))
        if float(r["M0"]) > 0 and float(r["INNER"]) > 0]
CSP = {"Directkin", "Dietmaier", "Neveu1", "Pramanik", "ponts-geo"} | {r["instancia"] for r in rows if r["instancia"].startswith("Trigexp")}
grupos = [("optimizacion", lambda i: i not in CSP), ("NCSP", lambda i: i in CSP)]
def per_inst(rs, k, ref):
    inst = sorted({r["instancia"] for r in rs})
    return {i: gm([float(r[k]) / float(r[ref]) for r in rs if r["instancia"] == i]) for i in inst}
for gname, sel in grupos:
    print(f"\n## {gname}\n\n| prof | inst | M1/M0 | M2/M0 | INNER/M0 | M0/INNER | M2/INNER |\n|---|---|---|---|---|---|---|")
    for d in ("0", "5", "10", "20"):
        rs = [r for r in rows if r["prof"] == d and sel(r["instancia"])]
        if not rs: continue
        v = {k: per_inst(rs, k, "M0") for k in ("M1", "M2", "INNER")}
        w = {k: per_inst(rs, k, "INNER") for k in ("M0", "M2")}
        print(f"| {d} | {len(v['M1'])} | " + " | ".join(f"{gm(list(v[k].values())):.3f}" for k in ("M1", "M2", "INNER"))
              + " | " + " | ".join(f"{gm(list(w[k].values())):.2f}" for k in ("M0", "M2")) + " |")
print("\n## Por instancia (todas las profundidades), M1/M0 y M2/M0; solo las que cambian (< 0,98)\n")
print("| instancia | M1/M0 | M2/M0 | INNER/M0 |\n|---|---|---|---|")
v = {k: per_inst(rows, k, "M0") for k in ("M1", "M2", "INNER")}
for i in sorted(v["M1"], key=lambda i: v["M2"][i]):
    if v["M2"][i] < 0.98: print(f"| {i} | {v['M1'][i]:.3f} | {v['M2'][i]:.3f} | {v['INNER'][i]:.3f} |")
print(f"\ninstancias: {len(v['M1'])}; con M1 < 0,98 M0: {sum(x < 0.98 for x in v['M1'].values())}; con M2 < 0,98 M0: {sum(x < 0.98 for x in v['M2'].values())}")
