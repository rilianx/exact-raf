"""Range overestimation on the test functions of Skalna and Hladik (2017, Sect. 5.2, Table 3).

    IBEXOPT=/path/to/ibexopt python experiments/sh_benchmark.py     # writes results/sh_benchmark.json

For each function and box, the range of f is enclosed with each multiplication rule and compared
with the exact range
computed by IbexOpt, with the measure of Skalna and Hladik (Eq. 21):
    O(f_rule, f_true) = (1 - wid(f_true) / wid(f_rule)) * 100.
Integer powers are evaluated in two ways: as repeated products with the rule itself
('product', as in the multiplication formulae compared by Skalna and Hladik) and with the
univariate Chebyshev envelope, the same for all rules ('chebyshev').
"""
import json
import os
import re
import subprocess
import tempfile

from common import RESULTS, ibexopt_path
from raf import set_box, raf_variables, set_product, set_square

RULES = ["faf2", "skalna", "exact"]

# name, number of variables, box (same interval for every variable), expression in x1..xn
FUNCTIONS = [
    ("f1 Beale", 2, (-5, 5),
     "(1.5 - x1 + x1*x2)^2 + (2.25 - x1 + x1*x2^2)^2 + (2.625 - x1 + x1*x2^3)^2"),
    ("f2 Booth", 2, (-10, 10), "(x1 + 2*x2 - 7)^2 + (2*x1 + x2 - 5)^2"),
    ("f3 Colville", 4, (-10, 10),
     "100*(x1^2 - x2)^2 + (x1 - 1)^2 + (x3 - 1)^2 + 90*(x3^2 - x4)^2"
     " + 10.1*((x2 - 1)^2 + (x4 - 1)^2) + 19.8*(x2 - 1)*(x4 - 1)"),
    ("f4 Dixon (n=5)", 5, (-5, 5),
     "(x1 - 1)^2 + 2*(2*x2^2 - x1)^2 + 3*(2*x3^2 - x2)^2 + 4*(2*x4^2 - x3)^2 + 5*(2*x5^2 - x4)^2"),
    ("f5 Goldstein-Price", 2, (-2, 2),
     "(1 + (x1 + x2 + 1.0)^2*(19 - 14*x1 + 13*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2))"
     "*(30 + (2*x1 - 3*x2)^2*(18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2))"),
    ("f6 Leon", 2, (-1.2, 1.2), "(1 - x1)^2 + 100*(x2 - x1^2)^2"),
    ("f7 Matyas", 2, (-10, 10), "0.26*(x1^2 + x2^2) + 0.48*x1*x2"),
    ("f8 Rosenbrock (n=4)", 4, (-5, 10),
     "100*(x2 - x1^2)^2 + (x1 - 1)^2 + 100*(x3 - x2^2)^2 + (x2 - 1)^2 + 100*(x4 - x3^2)^2 + (x3 - 1)^2"),
    ("f9 Six-hump camel", 2, (-5, 5), "4*x1^2 - 2.1*x1^4 + x1^6/3 + x1*x2 - 4*x2^2 + 4*x2^4"),
]


def evaluate(expr, n, values):
    env = {f"x{i + 1}": v for i, v in enumerate(values)}
    return eval(expr.replace("^", "**"), {}, env)


def ibex_range(expr, n, box, prec=1e-9, timeout=60):
    """Exact range [min f, max f] over the box, from two IbexOpt runs."""
    iexpr = re.sub(r"x(\d+)", lambda m: f"x[{int(m.group(1)) - 1}]", expr)
    dom = "(" + " ; ".join(f"[{a!r}, {b!r}]" for a, b in box) + ")"
    out = []
    for sign in ("", "-"):
        with tempfile.NamedTemporaryFile("w", suffix=".txt", delete=False) as fh:
            fh.write(f"variables\nx[{n}] in {dom};\n\nminimize {sign}({iexpr});")
            name = fh.name
        try:
            res = subprocess.run([ibexopt_path(), "-a", str(prec), "-r", str(prec), "-t", str(timeout), name],
                                 capture_output=True, text=True).stdout
        finally:
            for p in (name, name[:-4] + ".cov"):
                if os.path.exists(p):
                    os.remove(p)
        m = re.search(r"f\* in\s*\[(.*),(.*)\]", res)
        if not m:
            raise RuntimeError(res)
        lo, hi = float(m.group(1)), float(m.group(2))
        out.append(0.5 * (lo + hi))       # f* is enclosed to the requested precision
    return out[0], -out[1]


if __name__ == "__main__":
    records = []
    print(f"{'function':22s} {'true range':>26s} " + " ".join(f"{r:>8s}" for r in RULES)
          + " | " + " ".join(f"{r:>8s}" for r in RULES) + "   (overestimation %: powers as products | Chebyshev)")
    for name, n, (a, b), expr in FUNCTIONS:
        box = [(float(a), float(b))] * n
        set_box(box)
        lo, hi = ibex_range(expr, n, box)
        rec = {"function": name, "n": n, "box": [a, b], "true_range": [lo, hi]}
        for square in ("product", "chebyshev"):
            set_square(square)
            for rule in RULES:
                set_product(rule)
                form = evaluate(expr, n, raf_variables(n))
                flo, fhi = form.get_range()
                rec[f"{rule}_{square}"] = {"range": [flo, fhi], "overestimation": (1 - (hi - lo) / (fhi - flo)) * 100}
        set_product("exact")
        set_square("chebyshev")
        records.append(rec)
        print(f"{name:22s} [{lo:11.4g}, {hi:11.4g}] "
              + " ".join(f"{rec[r + '_product']['overestimation']:8.1f}" for r in RULES) + " | "
              + " ".join(f"{rec[r + '_chebyshev']['overestimation']:8.1f}" for r in RULES))
    with open(os.path.join(RESULTS, "sh_benchmark.json"), "w") as fh:
        json.dump(records, fh, indent=1)
