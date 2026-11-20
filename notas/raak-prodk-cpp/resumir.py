"""Resumen de comparacion.csv: fidelidad C++ frente a Python (anchos) y tiempos por tamano m."""
import csv, statistics as st
from math import exp, log
rows = list(csv.DictReader(open("comparacion.csv")))
both = [r for r in rows if r.get("py_ancho") not in (None, "")]
rat = [float(r["cpp_ancho"]) / float(r["py_ancho"]) for r in both if float(r["py_ancho"]) > 0]
print(f"formas: {len(rows)}; con Python (m <= 40): {len(both)}")
print(f"ancho C++ / Python: media geometrica {exp(sum(map(log, rat)) / len(rat)):.4f}, min {min(rat):.3f}, max {max(rat):.3f}; "
      f"difieren > 0,1 %: {sum(abs(x - 1) > 1e-3 for x in rat)}; C++ mas angosto > 0,1 %: {sum(x < 1 - 1e-3 for x in rat)}")
print("\n| m | formas | C++ mediana (µs) | C++ p90 (µs) | Python mediana (µs) | Python / C++ |\n|---|---|---|---|---|---|")
for lo, hi in ((2, 5), (6, 10), (11, 20), (21, 40), (41, 200)):
    rs = [r for r in rows if lo <= int(r["m"]) <= hi]
    if not rs: continue
    c = sorted(float(r["cpp_us"]) for r in rs)
    p = [float(r["py_us"]) for r in rs if r.get("py_us") not in (None, "")]
    pr = [float(r["py_us"]) / float(r["cpp_us"]) for r in rs if r.get("py_us") not in (None, "") and float(r["cpp_us"]) > 0]
    print(f"| {lo}–{hi} | {len(rs)} | {st.median(c):.1f} | {c[int(0.9 * len(c))]:.1f} | "
          + (f"{st.median(p):.0f} | {st.median(pr):.0f} |" if p else "— | — |"))
c = sorted(float(r["cpp_us"]) for r in rows)
print(f"\ntodas: C++ mediana {st.median(c):.1f} µs, p90 {c[int(0.9 * len(c))]:.1f} µs, max {c[-1]:.0f} µs")
