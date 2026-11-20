"""Resumen del barrido NCSP (barrido_solve.py): estado, celdas, tiempo y soluciones por instancia y variante;
y en las resueltas por todas, media geometrica de celdas y tiempo relativos a af2."""
import csv, sys
from collections import defaultdict
from math import exp, log
gm = lambda xs: exp(sum(map(log, xs)) / len(xs)) if xs else float("nan")
rows = list(csv.DictReader(open(sys.argv[1])))
V = list(dict.fromkeys(r["variante"] for r in rows))
by = defaultdict(dict)
for r in rows: by[r["instancia"]][r["variante"]] = r
print("| instancia | " + " | ".join(V) + " |\n|---|" + "---|" * len(V))
for i, d in sorted(by.items()):
    cell = lambda r: (f"{r['celdas']} c, {float(r['tiempo']):.1f} s, {r['solution']} sol" + ("" if r["estado"] == "SUCCESS" else f" ({r['estado']}, {r['unknown']} unk, {r['pending']} pend)")) if r else "-"
    print(f"| {i} | " + " | ".join(cell(d.get(v)) for v in V) + " |")
ok = [i for i, d in by.items() if all(v in d and d[v]["estado"] == "SUCCESS" for v in V)]
print(f"\nresueltas por todas: {len(ok)}")
print("| variante | resueltas (de 17) | celdas / af2 (geo) | tiempo / af2 (geo) |\n|---|---|---|---|")
for v in V:
    print(f"| {v} | {sum(by[i].get(v, {}).get('estado') == 'SUCCESS' for i in by)} | "
          f"{gm([float(by[i][v]['celdas']) / float(by[i]['af2']['celdas']) for i in ok]):.3f} | "
          f"{gm([max(1e-3, float(by[i][v]['tiempo'])) / max(1e-3, float(by[i]['af2']['tiempo'])) for i in ok]):.3f} |")
sol = [(i, {v: by[i][v]["solution"] for v in V}) for i in ok if len({by[i][v]["solution"] for v in V}) > 1]
print("\nsoluciones distintas entre variantes (resueltas por todas):", sol if sol else "ninguna")
