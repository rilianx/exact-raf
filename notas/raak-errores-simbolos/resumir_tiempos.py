"""Resumen de tiempos.csv: medianas y percentil 90 por componente (microsegundos, Python salvo kn6 en C++),
y relativos a la evaluacion M0 de la misma expresion."""
import csv, math, statistics as st
rows = list(csv.DictReader(open("tiempos.csv")))
K = ["t_M0", "t_M1", "t_M2", "t_M3c_eval", "t_tt", "t_homog", "t_prodk", "t_kn6_medio_lote", "t_sdp"]
print(f"{len(rows)} expresiones/cajas (caminata 1, prof 0 y 10); {sum('m' in r and r['m'] != '' for r in rows)} con parte cuadratica\n")
print("| componente | N | mediana (us) | p90 (us) | mediana / t_M0 |\n|---|---|---|---|---|")
for k in K:
    xs = [(float(r[k]), float(r["t_M0"])) for r in rows if r.get(k) not in (None, "") and not math.isnan(float(r[k]))]
    if not xs: continue
    v = sorted(x for x, _ in xs); q = sorted(x / t0 for x, t0 in xs)
    print(f"| {k[2:]} | {len(v)} | {1e6 * st.median(v):.0f} | {1e6 * v[int(0.9 * len(v))]:.0f} | {st.median(q):.2f} |")
ms = sorted(int(r["m"]) for r in rows if r.get("m") not in (None, ""))
print(f"\ntamano m de la forma homogeneizada: mediana {st.median(ms):.0f}, p90 {ms[int(0.9 * len(ms))]}, max {ms[-1]}")
print("prodk se mide solo con m <= 40 y sdp con m <= 60 (como en la medicion)")
