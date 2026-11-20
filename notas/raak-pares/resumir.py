"""Resumen de synq_pares.csv y reales_pares.csv: ancho / SH original y / exacto (medias geometricas)."""
import csv
from math import exp, log
gm = lambda xs: exp(sum(map(log, xs)) / len(xs))
M = ["SH", "STAR", "EIGP", "EIGSQ", "Qk0", "Qk6", "EIG", "PARES", "TODO", "EXACT"]
def load(p):
    rows = list(csv.DictReader(open(p)))
    for r in rows:
        for k in M[:7] + ["EXACT"]: r[k] = float(r[k])
        for k in ("STAR", "EIGP", "EIGSQ"):  # validez: ninguna descomposicion por debajo del exacto
            assert r[k] >= r["EXACT"] * (1 - 1e-9) - 1e-9 or r.get("exacto") == "False", (r["instancia"], k)
        r["PARES"] = min(r["SH"], r["STAR"], r["EIGP"])                  # solo reescrituras + SH
        r["TODO"] = min(r["PARES"], r["Qk6"], r["EIG"])                  # + k-nomial y espectral
    return rows
def tabla(titulo, grupos):
    print(f"\n{titulo}\n| grupo | N | " + " | ".join(M) + " |\n|---|---|" + "---|" * len(M))
    for nom, rs, ref in grupos:
        print(f"| {nom} | {len(rs)} | " + " | ".join(f"{gm([r[k] / r[ref] for r in rs]):.3f}" for k in M) + " |")
S = load("synq_pares.csv")
tabla("synq (78 instancias), ancho / SH original", [("todas", S, "SH")] + [(f"n = {n}", [r for r in S if r["n"] == n], "SH") for n in ("3", "5", "8", "10")])
tabla("synq, ancho / exacto", [("todas", S, "EXACT")] + [(f"n = {n}", [r for r in S if r["n"] == n], "EXACT") for n in ("3", "5", "8", "10")])
R = load("reales_pares.csv")
inst = list(dict.fromkeys(r["instancia"] for r in R))
print(f"\nreales: {len(R)} nodos en {len(inst)} instancias. Por instancia (media geometrica de sus nodos), ancho / SH original")
print("| instancia | nodos | " + " | ".join(M) + " |\n|---|---|" + "---|" * len(M))
per = []
for i in inst:
    rs = [r for r in R if r["instancia"] == i]
    v = {k: gm([r[k] / r["SH"] for r in rs]) for k in M}; per.append(v)
    print(f"| {i} | {len(rs)} | " + " | ".join(f"{v[k]:.3f}" for k in M) + " |")
print("| **media (instancias)** | | " + " | ".join(f"**{gm([v[k] for v in per]):.3f}**" for k in M) + " |")
for k in ("STAR", "EIGP"):
    print(f"reales: nodos con {k} < 0,99 SH: {sum(r[k] < 0.99 * r['SH'] for r in R)}; con {k} > 1,01 SH: {sum(r[k] > 1.01 * r['SH'] for r in R)}")
print(f"reales: nodos con PARES < 0,99 SH: {sum(r['PARES'] < 0.99 * r['SH'] for r in R)}; PARES <= 1,001 exacto: {sum(r['PARES'] <= 1.001 * r['EXACT'] for r in R)} de {len(R)}")
