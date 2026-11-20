"""Tabla del solver para el paper (results/table_solver.tex) y estadisticas en texto.
Cada conjunto: S (af2) de un barrido, SH y E (parche 0017) de otro, y opcionalmente SH con el cuadrado de S."""
import csv, glob, math, collections, os, sys
RES = "/home/user/ibex-local/results"
SETS = [  # nombre, tlim, dir S, dir SH/E, dir SH con cuadrado S (o None)
    ("synq", 300, "opt-sinteticas", "opt-synq-v17-t300", None),
    ("synt", 60, "opt-fun-t60", "opt-fun-v17-t60", "opt-fun-sqaf2-t60"),
    ("synp", 60, "opt-synp-t60", "opt-synp-v17-t60", "opt-synp-sqaf2-t60"),
    ("Ibex medium+hard", 200, "opt-ibexmh-t200", "opt-ibexmh-fast-t200", "opt-ibexmh-sqaf2-t200"),
]
def load(d):
    R = collections.defaultdict(dict)
    for f in glob.glob(os.path.join(RES, d, "rows", "*.csv")):
        for r in csv.reader(open(f)):
            if len(r) < 8 or r[0] == "instancia": continue
            R[r[0]][r[1].split("_")[1]] = dict(st=r[2], c=int(r[3]), lb=float(r[4]), ub=float(r[5]), t=float(r[6]))
    return R
gm = lambda x: math.exp(sum(map(math.log, x)) / len(x)) if x else float("nan")
lines = []; out = []
for name, tl, ds, dn, dq in SETS:
    S, N = load(ds), load(dn); Q = load(dq) if dq else None
    inst = [i for i in S if "af2" in S[i] and "sh" in N.get(i, {}) and "exact" in N.get(i, {}) and (Q is None or "sh" in Q.get(i, {}))]
    ok = lambda r: r["st"] == "SUCCESS"
    solved = [sum(ok(X[i][v]) for i in inst) for X, v in [(S, "af2"), (N, "sh"), (N, "exact")]]
    allok = [i for i in inst if ok(S[i]["af2"]) and ok(N[i]["sh"]) and ok(N[i]["exact"]) and (Q is None or ok(Q[i]["sh"]))]
    rc = lambda X, v: gm([max(X[i][v]["c"], 1) / max(S[i]["af2"]["c"], 1) for i in allok])
    rt = lambda X, v: gm([max(X[i][v]["t"], 1e-2) / max(S[i]["af2"]["t"], 1e-2) for i in allok])
    es = [N[i]["exact"]["c"] / max(N[i]["sh"]["c"], 1) for i in allok]
    less, more = sum(x < 0.99 for x in es), sum(x > 1.01 for x in es)
    q = f"{rc(Q,'sh'):.2f}" if Q else "--"
    lines.append(f"{name} & {len(inst)} & {tl} & {solved[0]} & {solved[1]} & {solved[2]} & {rc(N,'sh'):.2f} & {rc(N,'exact'):.2f} & {q} & {rt(N,'sh'):.2f} & {rt(N,'exact'):.2f} & {gm(es):.3f} & {less}/{more} \\\\")
    to = [i for i in inst if not any(ok(X[i][v]) for X, v in [(S, "af2"), (N, "sh"), (N, "exact")])]
    gap = lambda X, v: gm([max(X[i][v]["ub"] - X[i][v]["lb"], 1e-12) / max(S[i]["af2"]["ub"] - S[i]["af2"]["lb"], 1e-12) for i in to if math.isfinite(S[i]["af2"]["ub"]) and math.isfinite(X[i][v]["ub"])])
    bad = sum(1 for i in allok if max(S[i]['af2']['ub'], N[i]['sh']['ub'], N[i]['exact']['ub']) < min(S[i]['af2']['lb'], N[i]['sh']['lb'], N[i]['exact']['lb']) - 1e-6 * max(1, abs(S[i]['af2']['ub'])))
    out.append(f"{name}: inst={len(inst)} resueltas S/SH/E={solved} comunes={len(allok)} celdas SH={rc(N,'sh'):.3f} E={rc(N,'exact'):.3f} SHsqS={q} cpu SH={rt(N,'sh'):.3f} E={rt(N,'exact'):.3f} E/SH={gm(es):.3f} menos/mas={less}/{more} "
               f"no resueltas={len(to)} gap SH={gap(N,'sh'):.3f} E={gap(N,'exact'):.3f} inconsistencias={bad}")
tex = [r"\begin{tabular}{lrrrrrrrrrrrr}", r"\toprule",
       r" & & & \multicolumn{3}{c}{Solved} & \multicolumn{3}{c}{Nodes / S} & \multicolumn{2}{c}{CPU / S} & \multicolumn{2}{c}{E vs SH} \\",
       r"\cmidrule(lr){4-6}\cmidrule(lr){7-9}\cmidrule(lr){10-11}\cmidrule(lr){12-13}",
       r"Set & Inst. & $t_{\max}$ & S & SH & E & SH & E & SH$^\dagger$ & SH & E & nodes & $<$/$>$ \\", r"\midrule"] + lines + [r"\bottomrule", r"\end{tabular}"]
open("/home/user/raa-affine-arithmetic/results/table_solver.tex", "w").write("\n".join(tex) + "\n")
print("\n".join(out))
