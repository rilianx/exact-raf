"""Resumen de atribucion.csv: fraccion del radio del encierro (M1) por tipo de operacion.
Media aritmetica de las fracciones: primero en cada instancia, despues entre instancias. Ademas, la misma tabla
restringida a las expresiones/cajas con sobreestimacion relevante (M0 / INNER >= 1,5, de cajas.csv)."""
import csv
from collections import defaultdict
A = list(csv.DictReader(open("atribucion.csv")))
C = {(r["instancia"], r["expr"], r["caminata"], r["prof"]): float(r["M0"]) / float(r["INNER"])
     for r in csv.DictReader(open("cajas.csv")) if float(r["INNER"]) > 0}
K = ["lineal", "mul_lin", "mul_comp", "sqr_lin", "sqr_comp", "recip", "pow", "exp", "log", "sqrt", "sincos", "abs", "otro"]
CSP = {"Directkin", "Dietmaier", "Neveu1", "Pramanik", "ponts-geo"} | {r["instancia"] for r in A if r["instancia"].startswith("Trigexp")}
def mean_frac(rs):
    by = defaultdict(list)
    for r in rs: by[r["instancia"]].append(r)
    per = [{k: sum(float(r[k]) for r in v) / len(v) for k in K} for v in by.values()]
    return len(per), {k: sum(p[k] for p in per) / len(per) for k in K}
def tabla(titulo, filtro):
    print(f"\n{titulo}\n| conjunto | prof | inst | " + " | ".join(K) + " |\n|---|---|---|" + "---|" * len(K))
    for gname, sel in (("optimizacion", lambda i: i not in CSP), ("NCSP", lambda i: i in CSP)):
        for d in ("0", "5", "10", "20"):
            rs = [r for r in A if r["prof"] == d and sel(r["instancia"]) and filtro(r)]
            if not rs: continue
            n, m = mean_frac(rs)
            print(f"| {gname} | {d} | {n} | " + " | ".join(f"{100 * m[k]:.0f}" for k in K) + " |")
tabla("Fraccion del radio (%) por origen, todas las expresiones", lambda r: True)
key = lambda r: (r["instancia"], r["expr"], r["caminata"], r["prof"])
tabla("Fraccion del radio (%) por origen, solo expresiones con M0/INNER >= 1,5", lambda r: C.get(key(r), 0) >= 1.5)
print("\nOrigen dominante por instancia (todas las cajas), con la sobreestimacion media M0/INNER:")
by = defaultdict(list)
for r in A: by[r["instancia"]].append(r)
cnt = defaultdict(int)
for i, v in sorted(by.items()):
    m = {k: sum(float(r[k]) for r in v) / len(v) for k in K if k != "lineal"}
    top = sorted(m, key=lambda k: -m[k])[:2]
    ov = [C[key(r)] for r in v if key(r) in C]
    ovm = sum(ov) / len(ov) if ov else float("nan")
    cnt[top[0]] += 1
    print(f"  {i:16s} M0/INNER={ovm:6.2f}  lineal={100 * sum(float(r['lineal']) for r in v) / len(v):3.0f}%  "
          + ", ".join(f"{k}={100 * m[k]:.0f}%" for k in top))
print("\nInstancias por origen dominante (sin contar la parte lineal):", dict(sorted(cnt.items(), key=lambda t: -t[1])))
