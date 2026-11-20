"""Tabla conjunta de precision y tiempo, con M00 (fAF2 de Ibex por omision) como referencia.
Precision: ancho / M00 y ancho / INNER, media geometrica (expresiones/cajas de cada instancia, luego instancias),
todas las profundidades (0, 5, 10, 20). Tiempo: mediana por expresion/caja (tiempos.csv: caminata 1, prof 0 y 10;
Python salvo el k-nomial en C++) del costo total del metodo (evaluacion + encierro final), y su razon con M00."""
import csv, statistics as st
from collections import defaultdict
from math import exp, log, isnan
gm = lambda xs: exp(sum(map(log, xs)) / len(xs))
k4 = lambda r: (r["instancia"], r["expr"], r["caminata"], r["prof"])
A = {k4(r): r for r in csv.DictReader(open("cajas_m3.csv"))}
C = {k4(r): r for r in csv.DictReader(open("cajas_m3c.csv"))}
Z = {k4(r): r for r in csv.DictReader(open("cajas_m00.csv"))}
rows = []
for k, a in A.items():
    if k not in C or k not in Z: continue
    x = {"inst": k[0], "M00": float(Z[k]["M00"]), "INNER": float(a["INNER"])}
    for m in ("M0", "M1", "M2"): x[m] = float(a[m])
    for m in ("tt", "prodk", "kn6", "sdp"): x["M3c_" + m] = float(C[k]["M3c_" + m])
    x["M2∩M3c_prodk"] = min(x["M2"], x["M3c_prodk"])
    if x["M00"] > 0 and x["INNER"] > 0: rows.append(x)
CSP = {"Directkin", "Dietmaier", "Neveu1", "Pramanik", "ponts-geo"} | {r["inst"] for r in rows if r["inst"].startswith("Trigexp")}
def per(rs, k, ref):
    by = defaultdict(list)
    for r in rs: by[r["inst"]].append(r[k] / r[ref])
    return {i: gm(v) for i, v in by.items()}
# tiempos
T = list(csv.DictReader(open("tiempos.csv")))
T00 = {(r["instancia"], r["expr"], r["prof"]): float(r["t_M00"]) for r in csv.DictReader(open("cajas_m00.csv")) if r["caminata"] == "1"}
def ftime(r, parts):
    tot = 0.0
    for p in parts:
        v = r.get(p, "")
        if v in ("", None): return None
        v = float(v)
        if isnan(v): return None
        tot += v
    return tot
METH = [("M00", "fAF2 de Ibex (formula S)", ["t_M00"]),
        ("M0", "fAF2 + SH (parches 0016/0017)", ["t_M0"]),
        ("M1", "errores como simbolos + SH", ["t_M1"]),
        ("M2", "M1 + sumas de productos agregadas", ["t_M2"]),
        ("M3c_tt", "grado 2 + encierro barato", ["t_M3c_eval", "t_tt"]),
        ("M3c_kn6", "grado 2 + k-nomial k=6 (C++)", ["t_M3c_eval", "t_tt", "t_homog", "t_kn6_medio_lote"]),
        ("M3c_prodk", "grado 2 + bloques producto (Python)", ["t_M3c_eval", "t_tt", "t_homog", "t_prodk"]),
        ("M2∩M3c_prodk", "interseccion de M2 y la anterior", ["t_M2", "t_M3c_eval", "t_tt", "t_homog", "t_prodk"]),
        ("M3c_sdp", "grado 2 + SDP de Shor", ["t_M3c_eval", "t_tt", "t_homog", "t_sdp"])]
print(f"Precision: {len(rows)} expresiones/cajas, {len({r['inst'] for r in rows})} instancias. Tiempo: {len(T)} expresiones/cajas.\n")
print("| metodo | qué es | ancho/M00 optim. | ancho/M00 NCSP | inst. con mejora > 2 % (de 58) | ancho/INNER optim. | ancho/INNER NCSP | tiempo mediano (µs) | tiempo / M00 |")
print("|---|---|---|---|---|---|---|---|---|")
for m, desc, parts in METH:
    v = per(rows, m, "M00")
    po = gm([x for i, x in v.items() if i not in CSP]); pn = gm([x for i, x in v.items() if i in CSP])
    nb = sum(x < 0.98 for x in v.values())
    w = per(rows, m, "INNER")
    io = gm([x for i, x in w.items() if i not in CSP]); inn = gm([x for i, x in w.items() if i in CSP])
    ts, rs_ = [], []
    for r in T:
        r = dict(r); r["t_M00"] = T00.get((r["instancia"], r["expr"], r["prof"]), "")
        t = ftime(r, parts)
        if t is not None and r["t_M00"] != "":
            ts.append(t); rs_.append(t / float(r["t_M00"]))
    print(f"| {m} | {desc} | {po:.3f} | {pn:.3f} | {nb} | {io:.2f} | {inn:.2f} | {1e6 * st.median(ts):.0f} | {st.median(rs_):.2f} |")
print("\nNotas: tiempo de prodk y sdp solo donde se calcularon (m <= 40 y m <= 60); en el resto el metodo cae al encierro barato.")
