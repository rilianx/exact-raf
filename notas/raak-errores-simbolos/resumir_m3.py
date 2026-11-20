"""Resumen de cajas_m3.csv (mismas cajas para todos los modos). Media geometrica: expresiones/cajas de cada
instancia, despues instancias. M2|M3a = min(M2, M3a) (interseccion de encierros validos), idem M2|M3b."""
import csv, sys
from collections import defaultdict
from math import exp, log
gm = lambda xs: exp(sum(map(log, xs)) / len(xs)) if xs else float("nan")
rows = [r for r in csv.DictReader(open(sys.argv[1] if len(sys.argv) > 1 else "cajas_m3.csv"))
        if float(r["M0"]) > 0 and float(r["INNER"]) > 0]
for r in rows:
    for k in ("M0", "M1", "M2", "M3a", "M3b", "INNER"): r[k] = float(r[k])
    r["M2|M3a"] = min(r["M2"], r["M3a"]); r["M2|M3b"] = min(r["M2"], r["M3b"])
M = ["M1", "M2", "M3a", "M3b", "M2|M3a", "M2|M3b", "INNER"]
CSP = {"Directkin", "Dietmaier", "Neveu1", "Pramanik", "ponts-geo"} | {r["instancia"] for r in rows if r["instancia"].startswith("Trigexp")}
def per_inst(rs, k, ref="M0"):
    by = defaultdict(list)
    for r in rs: by[r["instancia"]].append(r[k] / r[ref])
    return {i: gm(v) for i, v in by.items()}
for gname, sel in (("optimizacion", lambda i: i not in CSP), ("NCSP", lambda i: i in CSP)):
    print(f"\n## {gname}: ancho / M0\n\n| prof | inst | " + " | ".join(M) + " | M0/INNER | (M2|M3b)/INNER |\n|---|---|" + "---|" * (len(M) + 2))
    for d in ("0", "5", "10", "20"):
        rs = [r for r in rows if r["prof"] == d and sel(r["instancia"])]
        if not rs: continue
        v = {k: per_inst(rs, k) for k in M}
        print(f"| {d} | {len(v['M1'])} | " + " | ".join(f"{gm(list(v[k].values())):.3f}" for k in M)
              + f" | {gm(list(per_inst(rs, 'M0', 'INNER').values())):.2f} | {gm(list(per_inst(rs, 'M2|M3b', 'INNER').values())):.2f} |")
v = {k: per_inst(rows, k) for k in M}
print("\n## Por instancia (todas las cajas), solo las que cambian (M2|M3b < 0,98)\n\n| instancia | " + " | ".join(M) + " |\n|---|" + "---|" * len(M))
for i in sorted(v["M1"], key=lambda i: v["M2|M3b"][i]):
    if v["M2|M3b"][i] < 0.98: print(f"| {i} | " + " | ".join(f"{v[k][i]:.3f}" for k in M) + " |")
for k in ("M2", "M3a", "M3b", "M2|M3a", "M2|M3b"):
    print(f"instancias con {k} < 0,98 M0: {sum(x < 0.98 for x in v[k].values())} de {len(v['M1'])}")
print("\nM3a peor que M1 (> 1,01) en", sum(1 for r in rows if r["M3a"] > 1.01 * r["M1"]), "de", len(rows), "expresiones/cajas")
tb = sorted(float(r["t_M3b"]) for r in rows); print(f"tiempo del encierro final M3b (Python, SDP+PRODK): mediana {tb[len(tb)//2]:.2f} s, p90 {tb[int(.9*len(tb))]:.2f} s")
