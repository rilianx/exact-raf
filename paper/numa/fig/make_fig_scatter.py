"""Figura de dispersion del solver: nodos por instancia, (a) SH contra S+, (b) E contra SH.
Instancias resueltas por las cuatro variantes; COCONUT: las 110 en que S pierde en algun producto (parche 0021).
    python make_fig_scatter.py [kindle]"""
import csv, glob, os, re, sys, collections
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
KINDLE = len(sys.argv) > 1 and sys.argv[1] == "kindle"
REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", ".."))
RES = ["/home/user/ibex-local/results/interleaved", "/home/user/ibex-local/results/interleaved_all"]
SW = set()
for f in glob.glob(os.path.join(REPO, "ibex", "results", "sdiag", "*.txt")):
    m = re.search(r"s_mas_ancho=(\d+)", open(f).read())
    if m and int(m.group(1)) > 0: SW.add(os.path.basename(f)[:-4])
R = collections.defaultdict(lambda: collections.defaultdict(dict))
for d in RES:
    for f in glob.glob(os.path.join(d, "rows", "*.csv")):
        for r in csv.reader(open(f)):
            st, inst = r[0].split("/")
            if st.startswith("ref") or st == "ibexmh": continue
            R[st][inst][r[1].split("_")[1]] = (r[2], max(int(r[3] or 0), 1))
SETS = [("synq", "synq", "o", "white"), ("synt", "synt", "^", "0.55"), ("synp", "synp", "s", "black"), ("COCONUT", "all", "x", "black")]
V = ["af2", "splus", "sh", "exact"]
plt.rcParams.update({"font.size": 11 if KINDLE else 9, "font.family": "serif", "mathtext.fontset": "stix"})
fig, ax = plt.subplots(2, 1, figsize=(3.4, 6.4)) if KINDLE else plt.subplots(1, 2, figsize=(6.8, 3.1))
pairs = [("splus", "sh", "nodes of S$^+$", "nodes of SH", "(a) SH versus S$^+$"), ("sh", "exact", "nodes of SH", "nodes of E", "(b) E versus SH")]
lo, hi = 1, 1
for a, (vx, vy, lx, ly, tt) in zip(ax, pairs):
    for name, key, mk, fc in SETS:
        D = R[key]
        for i in D:
            if "splus" not in D[i] and "af2" in D[i]: D[i]["splus"] = D[i]["af2"]
        I = [i for i in D if all(v in D[i] for v in V) and all(D[i][v][0] == "SUCCESS" for v in V) and (key != "all" or i in SW)]
        xs = [D[i][vx][1] for i in I]; ys = [D[i][vy][1] for i in I]
        hi = max(hi, max(xs + ys))
        kw = dict(marker=mk, s=16, linewidths=0.8, label=f"{name} ({len(I)})", zorder=3)
        if mk == "x": a.scatter(xs, ys, color=fc, **kw)
        else: a.scatter(xs, ys, facecolors=fc, edgecolors="black", **kw)
for a, (vx, vy, lx, ly, tt) in zip(ax, pairs):
    L = [0.7, hi * 1.5]
    a.plot(L, L, "k-", lw=0.7, zorder=1)
    for f, ls in [(0.5, ":"), (2.0, ":")]: a.plot(L, [f * L[0], f * L[1]], color="0.5", ls=ls, lw=0.7, zorder=1)
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(*L); a.set_ylim(*L); a.set_aspect("equal")
    a.set_xlabel(lx); a.set_ylabel(ly); a.set_title(tt, fontsize=10 if KINDLE else 9)
    a.grid(True, which="major", color="0.9", lw=0.5, zorder=0)
ax[0].legend(loc="upper left", fontsize=7.5 if not KINDLE else 8.5, frameon=True, framealpha=0.95, handletextpad=0.3)
fig.tight_layout()
name = "solver_scatter_kindle" if KINDLE else "solver_scatter"
out = os.path.dirname(os.path.abspath(__file__))
fig.savefig(os.path.join(out, name + ".pdf")); fig.savefig(os.path.join(out, name + ".png"), dpi=150)
