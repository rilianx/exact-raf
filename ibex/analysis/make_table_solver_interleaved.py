"""Tabla 5 desde la corrida intercalada (todas las variantes en la misma sesion): results/table_solver.tex.
Variantes: af2 (S), splus (S+), sh (SH), exact (E). En synq y synt no hay cuadrados: S+ = S."""
import csv, glob, math, collections, os, sys
RES = sys.argv[1] if len(sys.argv) > 1 else "/home/user/ibex-local/results/interleaved"
RES_ALL = sys.argv[2] if len(sys.argv) > 2 else "/home/user/ibex-local/results/interleaved_all"
SETS = [("synq", "synq", 300), ("synt", "synt", 60), ("synp", "synp", 60), ("COCONUT", "all", 60)]
# En el conjunto Ibex (COCONUT) solo se cuentan las instancias con al menos un producto en que el radio no lineal
# de S supera al de E (contador s_mas_ancho del parche 0021, corrida de 5 s con E): en las demas S, SH y E coinciden.
import re
SDIAG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "results", "sdiag")
SWIDER = set()
for f in glob.glob(os.path.join(SDIAG, "*.txt")):
    m = re.search(r"s_mas_ancho=(\d+)", open(f).read())
    if m and int(m.group(1)) > 0: SWIDER.add(os.path.basename(f)[:-4])
R = collections.defaultdict(lambda: collections.defaultdict(dict))
for f in glob.glob(os.path.join(RES, "rows", "*.csv")) + glob.glob(os.path.join(RES_ALL, "rows", "*.csv")):
    for r in csv.reader(open(f)):
        st, inst = r[0].split("/")
        if st.startswith("ref") or st == "ibexmh": continue
        R[st][inst][r[1].split("_")[1]] = dict(st=r[2], c=int(r[3] or 0), lb=float(r[4] or "nan"), ub=float(r[5] or "nan"), t=float(r[6] or "nan"))
gm = lambda x: math.exp(sum(map(math.log, x)) / len(x)) if x else float("nan")
ok = lambda r: r["st"] == "SUCCESS"
lines, out = [], []
for name, key, tl in SETS:
    D = R[key]
    for i in D:
        if "splus" not in D[i] and "af2" in D[i]: D[i]["splus"] = D[i]["af2"]
    V = ["af2", "splus", "sh", "exact"]
    inst = [i for i in D if all(v in D[i] for v in V)]
    ntot = len(inst)
    if key == "all": inst = [i for i in inst if i in SWIDER]
    ninst = f"{len(inst)}/{ntot}" if key == "all" else f"{len(inst)}"
    solved = [sum(ok(D[i][v]) for i in inst) for v in V]
    allok = [i for i in inst if all(ok(D[i][v]) for v in V)]
    rc = lambda v: gm([max(D[i][v]["c"], 1) / max(D[i]["splus"]["c"], 1) for i in allok])
    rt = lambda v: gm([max(D[i][v]["t"], 1e-2) / max(D[i]["splus"]["t"], 1e-2) for i in allok])
    es = [max(D[i]["exact"]["c"], 1) / max(D[i]["sh"]["c"], 1) for i in allok]
    less, more = sum(x < 0.99 for x in es), sum(x > 1.01 for x in es)
    lines.append(f"{name} & {ninst} & {tl} & {solved[0]} & {solved[1]} & {solved[2]} & {solved[3]} & {rc('af2'):.2f} & {rc('sh'):.2f} & {rc('exact'):.2f} & {rt('sh'):.2f} & {rt('exact'):.2f} & {gm(es):.3f} & {less}/{more} \\\\")
    to = [i for i in inst if not any(ok(D[i][v]) for v in ["splus", "sh", "exact"])]
    gap = lambda v: gm([max(D[i][v]["ub"] - D[i][v]["lb"], 1e-12) / max(D[i]["splus"]["ub"] - D[i]["splus"]["lb"], 1e-12) for i in to if all(math.isfinite(D[i][w]["ub"]) for w in ("splus", v))])
    bad = sum(1 for i in allok if max(D[i][v]["ub"] for v in V) < min(D[i][v]["lb"] for v in V) - 1e-6 * max(1, abs(D[i]["splus"]["ub"])))
    tpn = lambda v: gm([(D[i][v]["t"] / max(D[i][v]["c"], 1)) / (D[i]["splus"]["t"] / max(D[i]["splus"]["c"], 1)) for i in allok if D[i]["splus"]["t"] > 0.5])
    out.append(f"{name}: inst={len(inst)} resueltas S/S+/SH/E={solved} comunes={len(allok)} nodos/S+ S={rc('af2'):.3f} SH={rc('sh'):.3f} E={rc('exact'):.3f} "
               f"cpu/S+ S={rt('af2'):.3f} SH={rt('sh'):.3f} E={rt('exact'):.3f} t/nodo/S+ SH={tpn('sh'):.3f} E={tpn('exact'):.3f} E/SH={gm(es):.3f} menos/mas={less}/{more} "
               f"no resueltas={len(to)} gap SH={gap('sh'):.3f} E={gap('exact'):.3f} inconsistencias={bad}")
ref = [float(r[6]) for f in glob.glob(os.path.join(RES, "rows", "ref_*")) for r in csv.reader(open(f))]
tex = [r"\begin{tabular}{lrrrrrrrrrrrrr}", r"\toprule",
       r" & & & \multicolumn{4}{c}{Solved} & \multicolumn{3}{c}{Nodes / S$^+$} & \multicolumn{2}{c}{CPU / S$^+$} & \multicolumn{2}{c}{E vs SH} \\",
       r"\cmidrule(lr){4-7}\cmidrule(lr){8-10}\cmidrule(lr){11-12}\cmidrule(lr){13-14}",
       r"Set & Inst. & $t_{\max}$ & S & S$^+$ & SH & E & S & SH & E & SH & E & nodes & $<$/$>$ \\", r"\midrule"] + lines + [r"\bottomrule", r"\end{tabular}"]
open("/home/user/raa-affine-arithmetic/results/table_solver.tex", "w").write("\n".join(tex) + "\n")
print("\n".join(out))
