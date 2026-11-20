"""Resumen de nodos.csv: anchos relativos a SH (suma por producto) y al exacto, por nodo e instancia."""
import csv, sys
from collections import defaultdict
from math import exp, log
rows = list(csv.DictReader(open(sys.argv[1] if len(sys.argv) > 1 else "nodos.csv")))
f = lambda r, k: float(r[k])
gm = lambda xs: exp(sum(map(log, xs)) / len(xs)) if xs else float("nan")
M = ["S", "SH", "Qk0", "Qk2", "Qk6", "EIG", "BEST", "EXACT"]
for r in rows:
    r["BEST"] = min(f(r, "SH"), f(r, "Qk6"), f(r, "EIG"))
inst = sorted({r["instancia"] for r in rows}, key=lambda s: [r["instancia"] for r in rows].index(s))
print(f"Nodos suma con >= 2 terminos cuadraticos: {len(rows)} en {len(inst)} instancias; "
      f"comparten simbolos: {sum(r['comparte'] == 'True' for r in rows)}; exacto calculado (m <= 12): {sum(r['exacto'] == 'True' for r in rows)}\n")
print("Por instancia (media geometrica sobre sus nodos), ancho / SH")
print("| instancia | nodos | comparten | K medio | m max | S | Qk0 | Qk2 | Qk6 | EIG | BEST | EXACT |")
print("|---|---|---|---|---|---|---|---|---|---|---|---|")
per = {}
for i in inst:
    rs = [r for r in rows if r["instancia"] == i]
    per[i] = {k: gm([f(r, k) / f(r, "SH") for r in rs if f(r, "SH") > 0]) for k in M}
    ex = "" if all(r["exacto"] == "True" for r in rs) else " (≥)"
    print(f"| {i} | {len(rs)} | {sum(r['comparte']=='True' for r in rs)} | {sum(int(r['K']) for r in rs)/len(rs):.1f} | {max(int(r['m']) for r in rs)} | "
          + " | ".join(f"{per[i][k]:.3f}" for k in M[:1] + M[2:-1]) + f" | {per[i]['EXACT']:.3f}{ex} |")
print("\nMedia geometrica sobre instancias (cada instancia pesa 1), ancho / SH")
for nom, sel in [("todas", inst), ("con simbolos compartidos", [i for i in inst if any(r['comparte']=='True' for r in rows if r['instancia']==i)])]:
    print(f"  {nom} ({len(sel)}): " + ", ".join(f"{k}={gm([per[i][k] for i in sel]):.3f}" for k in M if k != "SH"))
print("\nNodos donde Qk6 < 0,9 SH:", sum(f(r,'Qk6') < 0.9*f(r,'SH') for r in rows), "de", len(rows))
print("Nodos donde Qk6 > SH:", sum(f(r,'Qk6') > f(r,'SH')*(1+1e-9) for r in rows))
print("Nodos donde Qk6 < 0,9 Qk0:", sum(f(r,'Qk6') < 0.9*f(r,'Qk0') for r in rows))
print("Nodos donde EIG < Qk6:", sum(f(r,'EIG') < f(r,'Qk6') for r in rows))
print("Nodos donde SH ya es exacto (SH <= 1,001 EXACT):", sum(f(r,'SH') <= 1.001*f(r,'EXACT') for r in rows))
