"""Tables and statistics of the multiplication rules (results/product_rules.json and
results/random_products.json, produced by product_rules.py and random_products.py).

    python experiments/make_tables_rules.py     # prints the statistics, writes results/table_rules_*.tex

Gap(M) = (W_M - W_OPT,M) / W_OPT,M * 100, where W_OPT,M is the exact width of the range of
f(x) - L_M(x) over the box, L_M being the linear part of method M (for AF2 and TM, the
linear part of the four rules, as in the stored records).
"""
import json
import os

import numpy as np
from scipy import stats

from common import RESULTS

# methods shown in the tables of the paper (M = old ZCAA and TM are kept in the records only)
METHODS = [("af2", "AF2"), ("faf2", "S"), ("skalna", "SH"), ("exact", "E")]
SIZE = {"small": "small", "medium": "medium", "large": "large"}


def gap(r, m):
    wo = r.get(f"{m}_opt_width", r["opt_width"])
    return (r[f"{m}_width"] - wo) / wo * 100


def table(records, caption, label):
    keys = list(dict.fromkeys((r["function"], r["size"]) for r in records))
    lines = [r"\begin{tabular}{ll" + "r" * len(METHODS) + "r}", r"\toprule",
             r" & & \multicolumn{" + str(len(METHODS)) + r"}{c}{Mean gap (\%)} & E $<$ SH \\",
             r"\cmidrule(lr){3-" + str(2 + len(METHODS)) + "}",
             "Expression & Size & " + " & ".join(n for _, n in METHODS) + r" & (runs) \\", r"\midrule"]
    prev = None
    for f, size in keys:
        sel = [r for r in records if r["function"] == f and r["size"] == size]
        cells = [f"{np.mean([gap(r, m) for r in sel]):.1f}" for m, _ in METHODS]
        strict = sum(r["exact_width"] < r["skalna_width"] * (1 - 1e-9) for r in sel)
        if prev is not None and f != prev:
            lines.append(r"\midrule")
        lines.append(f"{f if f != prev else ''} & {SIZE[size]} & " + " & ".join(cells) + f" & {strict}/{len(sel)} \\\\")
        prev = f
    lines += [r"\bottomrule", r"\end{tabular}"]
    return "\n".join([r"\begin{table}[t]", r"\centering", f"\\caption{{{caption}}}", f"\\label{{{label}}}",
                      r"\small", r"\setlength{\tabcolsep}{4pt}"] + lines + [r"\end{table}"])


def combined_table(groups):
    """One row per expression (gaps averaged over the three box sizes), several record sets."""
    lines = [r"\begin{table}[t]", r"\centering",
             r"\caption{Controlled products and constraints from a benchmark of nonlinear systems: mean gap (\%) of each "
             r"method with respect to the optimal residual width for its own linear part, over all boxes "
             r"(100 per size for the controlled cases and Trigexp1, 50 for the other constraints), number of boxes "
             r"where E is strictly narrower than SH, and largest width reduction of E with respect to SH.}",
             r"\label{tab:expressions}", r"\small", r"\setlength{\tabcolsep}{4pt}",
             r"\begin{tabular}{l" + "r" * len(METHODS) + "rr}", r"\toprule",
             r" & \multicolumn{" + str(len(METHODS)) + r"}{c}{Mean gap (\%)} & \multicolumn{2}{c}{E vs SH} \\",
             r"\cmidrule(lr){2-" + str(1 + len(METHODS)) + r"}\cmidrule(lr){" + str(2 + len(METHODS)) + "-"
             + str(3 + len(METHODS)) + "}",
             "Expression & " + " & ".join(n for _, n in METHODS) + r" & E $<$ SH & max red. \\", r"\midrule"]
    for k, (records, functions) in enumerate(groups):
        if k:
            lines.append(r"\midrule")
        for f in functions:
            sel = [r for r in records if r["function"] == f]
            cells = [f"{np.mean([gap(r, m) for r in sel]):.1f}" for m, _ in METHODS]
            strict = sum(r["exact_width"] < r["skalna_width"] * (1 - 1e-9) for r in sel)
            red = max((r["skalna_width"] - r["exact_width"]) / r["skalna_width"] * 100 for r in sel)
            lines.append(f"{f} & " + " & ".join(cells) + f" & {strict}/{len(sel)} & {max(red, 0.0):.1f}\\% \\\\")
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(lines)


def stats_block(name, records):
    e = np.array([r["exact_width"] for r in records])
    s = np.array([r["skalna_width"] for r in records])
    z = np.array([r["zcaa_width"] for r in records])
    red = (s - e) / s * 100
    strict = red > 1e-7
    p = stats.wilcoxon(s[strict], e[strict], alternative="greater").pvalue if strict.sum() > 0 else float("nan")
    print(f"[{name}] runs {len(records)}; E narrower than SH in {strict.sum()} ({strict.mean():.0%}), equal otherwise; "
          f"width reduction E vs SH: mean {red.mean():.2f}%, mean over strict {red[strict].mean():.2f}%, "
          f"max {red.max():.1f}%; Wilcoxon p = {p:.1e}")
    red_m = (z - s) / z * 100
    print(f"[{name}] SH narrower than M by {red_m.mean():.1f}% on average; M narrower than SH in "
          f"{np.sum(z < s * (1 - 1e-9))} runs; E never wider than M or SH: "
          f"{bool(np.all(e <= z * (1 + 1e-9)) and np.all(e <= s * (1 + 1e-9)))}")
    f2 = np.array([r["faf2_width"] for r in records])
    print(f"[{name}] S = fAF2 never narrower than M: {bool(np.all(z <= f2 * (1 + 1e-9)))}")


def random_table(rows):
    lines = [r"\begin{table}[t]", r"\centering",
             r"\caption{Single products of random forms with $n$ variables: mean width of the enclosures of rules S "
             r"and SH relative to E, and maximum of SH/E. $\rho$ scales the error radii of the operands; "
             r"500 products per row.}",
             r"\label{tab:random}", r"\small",
             r"\begin{tabular}{rrrrr}", r"\toprule",
             r" & & \multicolumn{2}{c}{Mean width / E} & \\", r"\cmidrule(lr){3-4}",
             r"$n$ & $\rho$ & S & SH & max SH/E \\", r"\midrule"]
    prev = None
    for r in rows:
        if r["n"] > 20:
            continue
        if prev is not None and r["n"] != prev:
            lines.append(r"\midrule")
        q = r["ratio_to_exact"]
        lines.append(f"{r['n'] if r['n'] != prev else ''} & {r['rho']:.1f} & {q['faf2']:.3f} & "
                     f"{q['skalna']:.3f} & {r['max_skalna_over_exact']:.3f} \\\\")
        prev = r["n"]
    lines += [r"\bottomrule", r"\end{tabular}", r"\end{table}"]
    return "\n".join(lines)


if __name__ == "__main__":
    with open(os.path.join(RESULTS, "product_rules.json")) as fh:
        data = json.load(fh)
    ctrl = data["controlled"] + data.get("mixed", [])
    for name, recs in (("controlled+mixed", ctrl), ("instances", data["instances"])):
        stats_block(name, recs)
    caption = (r"Controlled products (Cases 1--4) and the same expressions with mixed signs (1m--4m): mean gap of "
               r"each method with respect to the optimal residual width for its own linear part, "
               r"over 100 boxes per row, and number of boxes where E is strictly narrower than SH.")
    with open(os.path.join(RESULTS, "table_rules_controlled.tex"), "w") as fh:
        fh.write(table(ctrl, caption, "tab:controlled"))
    caption = (r"Constraints from a benchmark of nonlinear systems: mean gap of each method (\%) and number of boxes "
               r"where E is strictly narrower than SH (100 boxes per row for Trigexp1, 50 for the others).")
    with open(os.path.join(RESULTS, "table_rules_instances.tex"), "w") as fh:
        fh.write(table(data["instances"], caption, "tab:instances"))
    with open(os.path.join(RESULTS, "table_rules_expressions.tex"), "w") as fh:
        fh.write(combined_table([(ctrl, ["Case 1", "Case 1m", "Case 2", "Case 3", "Case 4"]),
                                 (data["instances"], ["Trigexp1", "Trigexp2 (1)", "Trigexp2 (2)", "Pramanik"])]))
    with open(os.path.join(RESULTS, "random_products.json")) as fh:
        with open(os.path.join(RESULTS, "table_rules_random.tex"), "w") as out:
            out.write(random_table(json.load(fh)))
    print("written results/table_rules_{controlled,instances,expressions,random}.tex")
