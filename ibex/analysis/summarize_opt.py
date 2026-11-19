"""Resumen de barridos ibexopt: media geometrica de celdas y CPU vs af2 (instancias resueltas por todas),
resueltas, y comparacion exact vs sh. Uso: python summarize_opt.py dir_rows [dir_rows ...]"""
import csv, glob, math, sys, collections
V = ["art_af2_k0", "art_raa_k0", "art_sh_k0", "art_exact_k0"]
for d in sys.argv[1:]:
    R = collections.defaultdict(dict)
    for f in glob.glob(d + "/*.csv"):
        for r in csv.reader(open(f)):
            if len(r) < 8 or r[0] == "instancia": continue
            R[r[0]][r[1]] = dict(st=r[2], cells=int(r[3]), lb=float(r[4]), ub=float(r[5]), cpu=float(r[6]))
    inst = [i for i in R if all(v in R[i] for v in V)]
    print(f"== {d}: {len(inst)} instancias completas")
    allok = [i for i in inst if all(R[i][v]["st"] == "SUCCESS" for v in V)]
    gm = lambda xs: math.exp(sum(map(math.log, xs)) / len(xs)) if xs else float("nan")
    print(f"   resueltas por todas: {len(allok)}")
    print("   variante   resueltas  celdas/af2  cpu/af2  (gm, resueltas por todas)")
    for v in V:
        ns = sum(R[i][v]["st"] == "SUCCESS" for i in inst)
        c = gm([max(R[i][v]["cells"], 1) / max(R[i]["art_af2_k0"]["cells"], 1) for i in allok])
        t = gm([max(R[i][v]["cpu"], 1e-3) / max(R[i]["art_af2_k0"]["cpu"], 1e-3) for i in allok])
        print(f"   {v:13s} {ns:5d}   {c:9.3f}  {t:8.3f}")
    # exact vs sh
    c = [R[i]["art_exact_k0"]["cells"] / R[i]["art_sh_k0"]["cells"] for i in allok]
    print(f"   exact/sh celdas gm={gm(c):.3f}  menos={sum(x<0.99 for x in c)} igual={sum(abs(x-1)<=0.01 for x in c)} mas={sum(x>1.01 for x in c)}")
    # gap en no resueltas (TIME_OUT por todas): ancho [lb,ub] relativo a af2
    to = [i for i in inst if all(R[i][v]["st"] != "SUCCESS" for v in V)]
    if to:
        print(f"   no resueltas por ninguna: {len(to)}; gap (ub-lb) gm vs af2:",
              " ".join(f"{v.split('_')[1]}={gm([max(R[i][v]['ub']-R[i][v]['lb'],1e-12)/max(R[i]['art_af2_k0']['ub']-R[i]['art_af2_k0']['lb'],1e-12) for i in to]):.3f}" for v in V))
    # consistencia f*
    bad = [i for i in allok if max(R[i][v]["ub"] for v in V) < min(R[i][v]["lb"] for v in V) - 1e-6*max(1,abs(R[i]["art_af2_k0"]["ub"]))]
    print(f"   inconsistencias f*: {len(bad)}")
