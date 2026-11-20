"""Tabla final del prototipo en C++ (proto_cpp.csv): precision (ancho / M00 y / INNER) y tiempo por expresion y
caja (evaluacion + encierro final), todo en C++. INNER viene de ../raak-errores-simbolos/cajas_m3.csv (mismas cajas).
Comprueba que ningun encierro quede por debajo del rango interior."""
import csv, statistics as st
from collections import defaultdict
from math import exp, log
gm = lambda xs: exp(sum(map(log, xs)) / len(xs))
k4 = lambda r: (r["instancia"], r["expr"], str(r["caminata"]), str(r["prof"]))
INNER = {k4(r): float(r["INNER"]) for r in csv.DictReader(open("../raak-errores-simbolos/cajas_m3.csv"))}
rows = []
for r in csv.DictReader(open("proto_cpp.csv")):
    k = k4(r)
    if k not in INNER or INNER[k] <= 0 or float(r["M00"]) <= 0: continue
    x = {c: float(v) for c, v in r.items() if c not in ("instancia", "expr", "caminata", "prof")}
    x["inst"], x["INNER"], x["prof"] = r["instancia"], INNER[k], int(r["prof"])
    for c in ("M00", "M0", "M3c_tt", "M3c_kn6", "M3c_pc", "M3c_full"):
        assert x[c] >= x["INNER"] * (1 - 1e-6), (k, c, x[c], x["INNER"])
    rows.append(x)
CSP = {"Directkin", "Dietmaier", "Neveu1", "Pramanik", "ponts-geo"} | {r["inst"] for r in rows if r["inst"].startswith("Trigexp")}
def per(k, ref):
    by = defaultdict(list)
    for r in rows: by[r["inst"]].append(r[k] / r[ref])
    return {i: gm(v) for i, v in by.items()}
METH = [("M00", "fAF2 de Ibex (formula S)", ["t_M00"]),
        ("M0", "fAF2 + SH", ["t_M0"]),
        ("M3c_tt", "grado 2 + encierro barato", ["t_M3c_eval", "t_tt"]),
        ("M3c_kn6", "grado 2 + k-nomial k=6", ["t_M3c_eval", "t_tt", "t_kn6"]),
        ("M3c_pc", "grado 2 + bloques producto (pares + completo)", ["t_M3c_eval", "t_tt", "t_pc"]),
        ("M3c_full", "grado 2 + bloques producto (completo)", ["t_M3c_eval", "t_tt", "t_full"])]
print(f"{len(rows)} expresiones/cajas, {len({r['inst'] for r in rows})} instancias; ningun encierro por debajo del rango interior\n")
print("| metodo | ancho/M00 optim. | ancho/M00 NCSP | inst. con mejora > 2 % (de 58) | ancho/INNER optim. | ancho/INNER NCSP | tiempo por expresion, mediana / p90 (µs) | tiempo / M00 (mediana) |")
print("|---|---|---|---|---|---|---|---|")
for m, desc, parts in METH:
    v = per(m, "M00"); w = per(m, "INNER")
    t = [sum(r[p] for p in parts) for r in rows]; q = [sum(r[p] for p in parts) / r["t_M00"] for r in rows if r["t_M00"] > 0]
    ts = sorted(t)
    print(f"| {desc} | {gm([x for i, x in v.items() if i not in CSP]):.3f} | {gm([x for i, x in v.items() if i in CSP]):.3f} | "
          f"{sum(x < 0.98 for x in v.values())} | {gm([x for i, x in w.items() if i not in CSP]):.2f} | {gm([x for i, x in w.items() if i in CSP]):.2f} | "
          f"{st.median(ts):.2f} / {ts[int(.9 * len(ts))]:.1f} | {st.median(q):.2f} |")
print("\nPor profundidad (ancho / M00, todas las instancias juntas; media por instancia):")
print("| prof | M0 | M3c_tt | M3c_kn6 | M3c_pc | M3c_full |\n|---|---|---|---|---|---|")
for d in (0, 5, 10, 20):
    sub = [r for r in rows if int(r["prof"]) == d]
    def pd(k):
        by = defaultdict(list)
        for r in sub: by[r["inst"]].append(r[k] / r["M00"])
        return gm([gm(v) for v in by.values()])
    print(f"| {d} | " + " | ".join(f"{pd(k):.3f}" for k in ("M0", "M3c_tt", "M3c_kn6", "M3c_pc", "M3c_full")) + " |")
