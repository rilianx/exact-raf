"""Resumen de qplib_raiz.csv y fqf_raiz.csv: anchos relativos al termino a termino (Qk0 = Ibex con monomios)
y al rango interior (INNER; exacto si n <= 12). Medias geometricas."""
import csv, sys
from math import exp, log
gm = lambda xs: exp(sum(map(log, xs)) / len(xs)) if xs else float("nan")
M = ["Qk0", "Qk2", "Qk3", "Qk6", "PRODK", "BLOQ", "EIG", "SDP", "INNER"]


def load(p):
    rows = list(csv.DictReader(open(p)))
    for r in rows:
        for k in M:
            if k in r: r[k] = float(r[k])
        r["BLOQ"] = min(r["Qk6"], r["PRODK"])   # mejor de los dos metodos por bloques
    return rows


def tabla(titulo, grupos, ref):
    print(f"\n{titulo}\n| grupo | N | " + " | ".join(M) + " |\n|---|---|" + "---|" * len(M))
    for nom, rs in grupos:
        if rs: print(f"| {nom} | {len(rs)} | " + " | ".join(f"{gm([r[k] / r[ref] for r in rs]):.3f}" for k in M) + " |")


Q = load("qplib_raiz.csv")
gq = [("todas", Q)] + [(t, [r for r in Q if r["tipo"] == t]) for t in sorted({r["tipo"] for r in Q})] \
   + [("diagonal cero", [r for r in Q if r["nnz_diag"] == "0"]), ("diagonal no nula", [r for r in Q if r["nnz_diag"] != "0"])]
tabla("QPLIB (58), ancho / Qk0 (termino a termino = Ibex)", gq, "Qk0")
tabla("QPLIB, ancho / INNER (rango interior: cota inferior del exacto)", gq, "INNER")
print("\nQPLIB: instancias con BLOQ < 0,95 Qk0:", sum(r["BLOQ"] < 0.95 * r["Qk0"] for r in Q),
      "| PRODK < Qk6:", sum(r["PRODK"] < r["Qk6"] for r in Q), "| SDP < BLOQ:", sum(r["SDP"] < r["BLOQ"] for r in Q), "de", len(Q))

F = load("fqf_raiz.csv")
gf = [("todas", F)] + [(f"n = {n}", [r for r in F if r["n_gen"] == n]) for n in ("10", "20", "40")] \
   + [(f"densidad {d}", [r for r in F if r["densidad"] == d]) for d in ("0.1", "0.3", "1.0")] \
   + [(f"diagonal {d}", [r for r in F if r["diagonal"] == d]) for d in ("cero", "pos", "mixta")]
tabla("fqf (81), ancho / Qk0", gf, "Qk0")
tabla("fqf, ancho / INNER (exacto con n = 10)", gf, "INNER")
print("\nfqf por diagonal y densidad, BLOQ / Qk0 (y PRODK < Qk6):")
for dg in ("cero", "pos", "mixta"):
    print("  " + dg + ": " + ", ".join(
        f"d={d}: {gm([r['BLOQ'] / r['Qk0'] for r in F if r['diagonal'] == dg and r['densidad'] == d]):.3f} "
        f"({sum(r['PRODK'] < r['Qk6'] for r in F if r['diagonal'] == dg and r['densidad'] == d)}/9)" for d in ("0.1", "0.3", "1.0")))
