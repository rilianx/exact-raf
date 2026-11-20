"""Resumen de synq_prodk.csv y reales_prodk.csv: reduccion por productos (PRODK) frente al k-nomial (Qk6)."""
import csv
from math import exp, log
gm = lambda xs: exp(sum(map(log, xs)) / len(xs))
M = ["SH", "STAR", "Qk6", "PRODK", "MIN2", "MINALL"]
def load(p):
    rows = list(csv.DictReader(open(p)))
    for r in rows:
        for k in ("SH", "STAR", "EIGP", "Qk0", "Qk6", "EIG", "EXACT", "PRODK"): r[k] = float(r[k])
        r["MIN2"] = min(r["PRODK"], r["Qk6"])
        r["MINALL"] = min(r["PRODK"], r["Qk6"], r["EIG"], r["STAR"], r["EIGP"], r["SH"])
    return rows
S = load("synq_prodk.csv")
dens = lambda r: r["instancia"].split("_")[2]
print("synq, ancho / exacto (media geometrica); MIN2 = min(PRODK, Qk6); MINALL = min de todas")
print("| grupo | N | " + " | ".join(M) + " | PRODK < Qk6 |\n|---|---|" + "---|" * (len(M) + 1))
for nom, rs in [("todas", S)] + [(f"n = {n}", [r for r in S if r["n"] == n]) for n in ("3", "5", "8", "10")] \
        + [(f"d = 0,{d[1]}", [r for r in S if dens(r) == d]) for d in ("d20", "d50", "d90")]:
    print(f"| {nom} | {len(rs)} | " + " | ".join(f"{gm([r[k] / r['EXACT'] for r in rs]):.3f}" for k in M)
          + f" | {sum(r['PRODK'] < r['Qk6'] for r in rs)} |")
R = load("reales_prodk.csv")
inst = list(dict.fromkeys(r["instancia"] for r in R))
per = {i: {k: gm([r[k] / r["SH"] for r in R if r["instancia"] == i]) for k in M + ["EXACT"]} for i in inst}
print(f"\nreales ({len(R)} nodos, {len(inst)} instancias), ancho / SH; media geometrica por instancia")
print("| | " + " | ".join(M + ["EXACT"]) + " |\n|---|" + "---|" * (len(M) + 1))
print("| media | " + " | ".join(f"{gm([per[i][k] for i in inst]):.3f}" for k in M + ["EXACT"]) + " |")
print("\ninstancias donde PRODK o Qk6 bajan de 0,99 SH:")
print("| instancia | " + " | ".join(M + ["EXACT"]) + " |\n|---|" + "---|" * (len(M) + 1))
for i in inst:
    if min(per[i]["PRODK"], per[i]["Qk6"]) < 0.99:
        print(f"| {i} | " + " | ".join(f"{per[i][k]:.3f}" for k in M + ["EXACT"]) + " |")
print(f"\nreales: nodos con PRODK < Qk6: {sum(r['PRODK'] < r['Qk6'] * (1 - 1e-9) for r in R)}; "
      f"PRODK > Qk6: {sum(r['PRODK'] > r['Qk6'] * (1 + 1e-9) for r in R)}; PRODK > SH: {sum(r['PRODK'] > r['SH'] * (1 + 1e-9) for r in R)} de {len(R)}")
