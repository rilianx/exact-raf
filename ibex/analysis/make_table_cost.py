"""Table 3 of the paper (time per product, microseconds) from ibex/results/bench_mult_0019.csv
(output of the bench_mult program built by ibex/build.sh): writes results/table_cost.tex."""
import csv, os
HERE = os.path.dirname(os.path.abspath(__file__))
rows = {(r["n"], r["errores"]): r for r in csv.DictReader(open(os.path.join(HERE, "..", "results", "bench_mult_0019.csv")))}
ns = sorted({int(n) for n, _ in rows})
f = lambda x: f"{float(x):.2f}"
out = [r"\begin{tabular}{rrrrrrr}", r"\toprule",
       r" & \multicolumn{3}{c}{Without errors} & \multicolumn{3}{c}{With errors} \\",
       r"\cmidrule(lr){2-4}\cmidrule(lr){5-7}", r"$n$ & S & SH & E & S & SH & E \\", r"\midrule"]
for n in ns:
    a, b = rows[(str(n), "no")], rows[(str(n), "si")]
    out.append(f"{n} & " + " & ".join(f(r[k]) for r in (a, b) for k in ("af2_us", "sh_us", "exact_us")) + r" \\")
out += [r"\bottomrule", r"\end{tabular}"]
open(os.path.join(HERE, "..", "..", "results", "table_cost.tex"), "w").write("\n".join(out) + "\n")
