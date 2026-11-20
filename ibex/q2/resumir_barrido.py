"""Resumen de un barrido de barrido_opt.py: por conjunto, instancias resueltas por variante, y en las resueltas
por todas las variantes, media geometrica de celdas y CPU relativas a af2 (Ibex por omision). Comprueba la
coherencia de f*: ninguna cota inferior por encima de la mejor cota superior (con tolerancia relativa 1e-6)."""
import csv, sys
from collections import defaultdict
from math import exp, log
gm = lambda xs: exp(sum(map(log, xs)) / len(xs)) if xs else float("nan")
rows = list(csv.DictReader(open(sys.argv[1])))
V = list(dict.fromkeys(r["variante"] for r in rows))
by = defaultdict(dict)
for r in rows: by[(r["conjunto"], r["instancia"])][r["variante"]] = r
def grupo(c, i):
    if i.startswith("synq"): return "synq"
    if i.startswith("fqf"): return "fqf"
    return "ibex medium+hard y coconut"
bad = 0
for k, d in by.items():
    ubs = [float(r["f_ub"]) for r in d.values() if r["f_ub"] not in ("", "inf", "-inf")]
    if not ubs: continue
    ub = min(ubs)
    for v, r in d.items():
        if r["f_lb"] not in ("", "-inf") and float(r["f_lb"]) > ub + 1e-6 * max(1, abs(ub)):
            bad += 1; print("INCOHERENTE", k, v, r["f_lb"], ub)
print(f"coherencia de f*: {bad} casos con cota inferior por encima de la mejor cota superior\n")
for g in ("ibex medium+hard y coconut", "synq", "fqf"):
    ks = [k for k in by if grupo(*k) == g and all(v in by[k] for v in V)]
    if not ks: continue
    solved = {v: sum(by[k][v]["estado"] == "SUCCESS" for k in ks) for v in V}
    allk = [k for k in ks if all(by[k][v]["estado"] == "SUCCESS" for v in V)]
    allk = [k for k in allk if float(by[k]["af2"]["cpu_s"]) >= 0.05]
    print(f"## {g}: {len(ks)} instancias; resueltas por todas (cpu af2 >= 0,05 s): {len(allk)}\n")
    print("| variante | resueltas | celdas / af2 (geo) | CPU / af2 (geo) | celdas por s |\n|---|---|---|---|---|")
    for v in V:
        c = gm([float(by[k][v]["celdas"]) / float(by[k]["af2"]["celdas"]) for k in allk])
        t = gm([max(1e-3, float(by[k][v]["cpu_s"])) / max(1e-3, float(by[k]["af2"]["cpu_s"])) for k in allk])
        cps = gm([float(by[k][v]["celdas"]) / max(1e-3, float(by[k][v]["cpu_s"])) for k in allk])
        print(f"| {v} | {solved[v]} | {c:.3f} | {t:.3f} | {cps:.0f} |")
    print()
