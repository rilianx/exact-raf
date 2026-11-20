"""Une cajas_m3.csv (M0, M1, M2, INNER) con cajas_m3c.csv (M3 original 'o' y centrado 'c', con encierros finales
tt, prodk, kn6, sdp) y resume ancho / M0 (media geometrica por instancia y luego entre instancias). Comprueba
que ningun encierro quede por debajo del rango interior."""
import csv
from collections import defaultdict
from math import exp, log
gm = lambda xs: exp(sum(map(log, xs)) / len(xs)) if xs else float("nan")
A = {(r["instancia"], r["expr"], r["caminata"], r["prof"]): r for r in csv.DictReader(open("cajas_m3.csv"))}
rows = []
for r in csv.DictReader(open("cajas_m3c.csv")):
    a = A.get((r["instancia"], r["expr"], r["caminata"], r["prof"]))
    if a is None or float(a["M0"]) <= 0 or float(a["INNER"]) <= 0: continue
    x = {k: float(a[k]) for k in ("M0", "M1", "M2", "INNER")}
    for k, v in r.items():
        if k.startswith("M3"): x[k] = float(v)
    for k in [k for k in x if k.startswith("M3")]:
        assert x[k] >= x["INNER"] * (1 - 1e-6), (r["instancia"], r["expr"], r["prof"], k, x[k], x["INNER"])
    # combinaciones baratas: interseccion con M2 (encierros validos)
    x["M2|M3c_tt"] = min(x["M2"], x["M3c_tt"]); x["M2|M3c_kn6"] = min(x["M2"], x["M3c_kn6"]); x["M2|M3c_sdp"] = min(x["M2"], x["M3c_sdp"])
    x["inst"], x["prof"] = r["instancia"], r["prof"]
    rows.append(x)
M = ["M1", "M2", "M3o_tt", "M3c_tt", "M3c_prodk", "M3c_kn6", "M3c_sdp", "M2|M3c_tt", "M2|M3c_kn6", "M2|M3c_sdp"]
CSP = {"Directkin", "Dietmaier", "Neveu1", "Pramanik", "ponts-geo"} | {r["inst"] for r in rows if r["inst"].startswith("Trigexp")}
def per(rs, k, ref="M0"):
    by = defaultdict(list)
    for r in rs: by[r["inst"]].append(r[k] / r[ref])
    return {i: gm(v) for i, v in by.items()}
print(f"filas: {len(rows)}; instancias: {len({r['inst'] for r in rows})}; ninguna por debajo del rango interior")
for g, sel in (("optimizacion", lambda i: i not in CSP), ("NCSP", lambda i: i in CSP)):
    print(f"\n## {g}: ancho / M0\n\n| prof | inst | " + " | ".join(M) + " | M0/INNER | (M2|M3c_tt)/INNER | (M2|M3c_sdp)/INNER |\n|---|---|" + "---|" * (len(M) + 3))
    for d in ("0", "5", "10", "20"):
        rs = [r for r in rows if r["prof"] == d and sel(r["inst"])]
        if not rs: continue
        v = {k: per(rs, k) for k in M}
        print(f"| {d} | {len(v['M1'])} | " + " | ".join(f"{gm(list(v[k].values())):.3f}" for k in M)
              + " | " + " | ".join(f"{gm(list(per(rs, k, 'INNER').values())):.2f}" for k in ("M0", "M2|M3c_tt", "M2|M3c_sdp")) + " |")
v = {k: per(rows, k) for k in M}
print("\n## Por instancia (todas las cajas), las que cambian con M2|M3c_tt o M2|M3c_sdp (< 0,98)\n\n| instancia | " + " | ".join(M) + " |\n|---|" + "---|" * len(M))
for i in sorted(v["M1"], key=lambda i: v["M2|M3c_sdp"][i]):
    if min(v["M2|M3c_tt"][i], v["M2|M3c_sdp"][i]) < 0.98: print(f"| {i} | " + " | ".join(f"{v[k][i]:.3f}" for k in M) + " |")
for k in M: print(f"instancias con {k} < 0,98 M0: {sum(x < 0.98 for x in v[k].values())} de {len(v['M1'])}")
print("\nM3c_tt peor que M1 (> 1,01):", sum(r["M3c_tt"] > 1.01 * r["M1"] for r in rows), "| M3o_tt peor que M1:", sum(r["M3o_tt"] > 1.01 * r["M1"] for r in rows), "de", len(rows))
