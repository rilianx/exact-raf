"""Table 2 of the paper (overestimation on the test functions of Skalna and Hladik) from
results/sh_benchmark.json (output of experiments/sh_benchmark.py): writes results/table_sh.tex.
E is in bold when it is narrower than SH."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "..", "results")
NAMES = {"Dixon (n=5)": r"Dixon ($n=5$)", "Goldstein-Price": "Goldstein--Price", "Rosenbrock (n=4)": r"Rosenbrock ($n=4$)"}
out = [r"\begin{tabular}{llrrr}", r"\toprule", r"Function & Box & S & SH & E \\", r"\midrule"]
for r in json.load(open(os.path.join(RES, "sh_benchmark.json"))):
    num, name = r["function"].split(" ", 1)
    lo, hi = r["box"]
    s, sh, e = (r[k + "_product"]["overestimation"] for k in ("faf2", "skalna", "exact"))
    es = f"{e:.1f}"
    if round(e, 1) < round(sh, 1): es = r"\textbf{" + es + "}"
    out.append(f"${num[0]}_{num[1:]}$ {NAMES.get(name, name)} & $[{lo:g},{hi:g}]^{r['n']}$ & {s:.1f} & {sh:.1f} & {es} \\\\")
out += [r"\bottomrule", r"\end{tabular}"]
open(os.path.join(RES, "table_sh.tex"), "w").write("\n".join(out) + "\n")
print("\n".join(out))
