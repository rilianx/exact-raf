"""Tabla del solver para el paper (results/table_solver.tex) y estadisticas en texto.
Cada conjunto: S (af2) de un barrido, SH y E (parche 0017) de otro, y opcionalmente SH con el cuadrado de S."""
import csv, glob, math, collections, os, sys
RES = "/home/user/ibex-local/results"
SETS = [  # nombre, tlim, dir S, dir SH/E, dir S+ (S con cuadrado de Chebyshev; None = sin cuadrados, S+ = S)
    # (nombre, tlim, dir S, dir SH, dir S+, dir E con el parche 0019: un recorrido con dos desplazamientos)
    ("synq", 300, "opt-sinteticas", "opt-synq-v17-t300", None, "opt-synq-v19-t300"),
    ("synt", 60, "opt-fun-t60", "opt-fun-v17-t60", None, "opt-fun-v19-t60"),
    ("synp", 60, "opt-synp-t60", "opt-synp-v17-t60", "opt-synp-splus-t60", "opt-synp-v19-t60"),
    ("Ibex medium+hard", 200, "opt-ibexmh-t200", "opt-ibexmh-fast-t200", "opt-ibexmh-splus-t200", "opt-ibexmh-v19-t200"),
]
def load(d):
    R = collections.defaultdict(dict)
    for f in glob.glob(os.path.join(RES, d, "rows", "*.csv")):
        for r in csv.reader(open(f)):
            if len(r) < 8 or r[0] == "instancia": continue
            R[r[0].split("/")[-1]][r[1].split("_")[1]] = dict(st=r[2], c=int(r[3]), lb=float(r[4]), ub=float(r[5]), t=float(r[6]))
    return R
gm = lambda x: math.exp(sum(map(math.log, x)) / len(x)) if x else float("nan")
lines = []; out = []
for name, tl, ds, dn, dq, de in SETS:
    S, N = load(ds), load(dn); Q = load(dq) if dq else S
    Ex = load(de)
    for i in N:
        if "exact" in Ex.get(i, {}): N[i]["exact"] = Ex[i]["exact"]
        else: N[i].pop("exact", None)
    ok = lambda r: r["st"] == "SUCCESS"
    inst = [i for i in S if "af2" in S[i] and "af2" in Q.get(i, {}) and "sh" in N.get(i, {}) and "exact" in N.get(i, {})]
    solved = [sum(ok(X[i][v]) for i in inst) for X, v in [(S, "af2"), (Q, "af2"), (N, "sh"), (N, "exact")]]
    allok = [i for i in inst if ok(S[i]["af2"]) and ok(Q[i]["af2"]) and ok(N[i]["sh"]) and ok(N[i]["exact"])]
    base = lambda i: Q[i]["af2"]
    rc = lambda X, v: gm([max(X[i][v]["c"], 1) / max(base(i)["c"], 1) for i in allok])
    rt = lambda X, v: gm([max(X[i][v]["t"], 1e-2) / max(base(i)["t"], 1e-2) for i in allok])
    es = [N[i]["exact"]["c"] / max(N[i]["sh"]["c"], 1) for i in allok]
    less, more = sum(x < 0.99 for x in es), sum(x > 1.01 for x in es)
    lines.append(f"{name} & {len(inst)} & {tl} & {solved[0]} & {solved[1]} & {solved[2]} & {solved[3]} & {rc(S,'af2'):.2f} & {rc(N,'sh'):.2f} & {rc(N,'exact'):.2f} & {rt(N,'sh'):.2f} & {rt(N,'exact'):.2f} & {gm(es):.3f} & {less}/{more} \\\\")
    to = [i for i in inst if not any(ok(X[i][v]) for X, v in [(Q, "af2"), (N, "sh"), (N, "exact")])]
    gap = lambda X, v: gm([max(X[i][v]["ub"] - X[i][v]["lb"], 1e-12) / max(base(i)["ub"] - base(i)["lb"], 1e-12) for i in to if math.isfinite(base(i)["ub"]) and math.isfinite(X[i][v]["ub"])])
    out.append(f"{name}: inst={len(inst)} resueltas S/S+/SH/E={solved} comunes={len(allok)} nodos/S+: S={rc(S,'af2'):.3f} SH={rc(N,'sh'):.3f} E={rc(N,'exact'):.3f} "
               f"cpu/S+: SH={rt(N,'sh'):.3f} E={rt(N,'exact'):.3f} E/SH={gm(es):.3f} menos/mas={less}/{more} no resueltas={len(to)} gap SH={gap(N,'sh'):.3f} E={gap(N,'exact'):.3f}")
tex = [r"\begin{tabular}{lrrrrrrrrrrrrr}", r"\toprule",
       r" & & & \multicolumn{4}{c}{Solved} & \multicolumn{3}{c}{Nodes / S$^+$} & \multicolumn{2}{c}{CPU / S$^+$} & \multicolumn{2}{c}{E vs SH} \\",
       r"\cmidrule(lr){4-7}\cmidrule(lr){8-10}\cmidrule(lr){11-12}\cmidrule(lr){13-14}",
       r"Set & Inst. & $t_{\max}$ & S & S$^+$ & SH & E & S & SH & E & SH & E & nodes & $<$/$>$ \\", r"\midrule"] + lines + [r"\bottomrule", r"\end{tabular}"]
open("/home/user/raa-affine-arithmetic/results/table_solver.tex", "w").write("\n".join(tex) + "\n")
print("\n".join(out))
