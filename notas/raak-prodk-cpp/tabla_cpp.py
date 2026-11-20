"""Precision final de las expresiones con M3 centrado + bloques producto en C++ (variantes de soportes), contra
Ibex por omision (M00), y tiempo del encierro final en C++. Une: cajas_m00.csv (M00), cajas_m3.csv (INNER, M0..M2),
cajas_m3c.csv (M3c_tt, M3c_kn6, M3c_sdp), formas_meta.csv + salidas del driver C++ (out_v_*.txt).
El encierro final de cada variante es min(M3c_tt, ancho C++)."""
import csv, statistics as st, sys
from collections import defaultdict
from math import exp, log
D = Path = None
base = "../raak-errores-simbolos/"
gm = lambda xs: exp(sum(map(log, xs)) / len(xs))
k4 = lambda r: (r["instancia"], r["expr"], str(r["caminata"]), str(r["prof"]))
Z = {k4(r): float(r["M00"]) for r in csv.DictReader(open(base + "cajas_m00.csv"))}
A = {k4(r): r for r in csv.DictReader(open(base + "cajas_m3.csv"))}
C = {k4(r): r for r in csv.DictReader(open(base + "cajas_m3c.csv"))}
meta = list(csv.DictReader(open("formas_meta.csv")))
VAR = [("0_0", "solo pares"), ("0_12", "pares + completo (m <= 12)"), ("40_0", "pares + triples"), ("40_12", "pares + triples + completo")]
outs = {v: [l.split() for l in open(f"{sys.argv[1]}/out_v_{v}.txt")] for v, _ in VAR}
cpp = {v: {k4(m): (float(o[1]) - float(o[0]), float(o[3])) for m, o in zip(meta, outs[v])} for v, _ in VAR}
rows = []
for k, a in A.items():
    if k not in C or k not in Z or Z[k] <= 0 or float(a["INNER"]) <= 0: continue
    x = {"inst": k[0], "M00": Z[k], "INNER": float(a["INNER"]), "M0": float(a["M0"]), "M3c_tt": float(C[k]["M3c_tt"]),
         "M3c_kn6": float(C[k]["M3c_kn6"]), "M3c_sdp": float(C[k]["M3c_sdp"])}
    for v, _ in VAR:
        w = cpp[v].get(k)
        x["C_" + v] = min(x["M3c_tt"], w[0]) if w else x["M3c_tt"]
        x["t_" + v] = w[1] if w else 0.0
    rows.append(x)
CSP = {"Directkin", "Dietmaier", "Neveu1", "Pramanik", "ponts-geo"} | {r["inst"] for r in rows if r["inst"].startswith("Trigexp")}
def per(k, ref):
    by = defaultdict(list)
    for r in rows: by[r["inst"]].append(r[k] / r[ref])
    return {i: gm(v) for i, v in by.items()}
print(f"{len(rows)} expresiones/cajas, {len({r['inst'] for r in rows})} instancias\n")
print("| encierro final de M3 centrado | ancho/M00 optim. | ancho/M00 NCSP | inst. con mejora > 2 % | ancho/INNER optim. | ancho/INNER NCSP | tiempo C++ del encierro (mediana / p90, µs) |")
print("|---|---|---|---|---|---|---|")
def line(name, key, tkey=None):
    v = per(key, "M00"); w = per(key, "INNER")
    po = gm([x for i, x in v.items() if i not in CSP]); pn = gm([x for i, x in v.items() if i in CSP])
    io = gm([x for i, x in w.items() if i not in CSP]); inn = gm([x for i, x in w.items() if i in CSP])
    if tkey:
        t = sorted(r[tkey] for r in rows if r[tkey] > 0); ts = f"{st.median(t):.1f} / {t[int(.9 * len(t))]:.0f}"
    else: ts = "—"
    print(f"| {name} | {po:.3f} | {pn:.3f} | {sum(x < 0.98 for x in v.values())} | {io:.2f} | {inn:.2f} | {ts} |")
line("(fAF2 + SH, M0, referencia)", "M0")
line("encierro barato (tt)", "M3c_tt")
line("k-nomial k = 6 (C++)", "M3c_kn6")
for v, name in VAR: line("bloques producto C++: " + name, "C_" + v, "t_" + v)
line("SDP (Python, referencia)", "M3c_sdp")
