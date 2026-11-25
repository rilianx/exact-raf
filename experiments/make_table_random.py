"""Table 1 of the paper (single products of random forms) from results/random_products.json
(output of experiments/random_products.py): writes results/table_rules_random.tex.

    python experiments/make_table_random.py
"""
import json
import os

from common import RESULTS


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
    with open(os.path.join(RESULTS, "random_products.json")) as fh:
        with open(os.path.join(RESULTS, "table_rules_random.tex"), "w") as out:
            out.write(random_table(json.load(fh)))
    print("written results/table_rules_random.tex")
