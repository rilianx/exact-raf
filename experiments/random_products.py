"""Single products of random reduced affine forms: residual width of each rule relative to the
exact one, and running time per product.

    python experiments/random_products.py      # writes results/random_products.json

The box is [-1, 1]^n, so that the variables coincide with the noise symbols. The linear
coefficients are standard normal; the centers are normal with standard deviation 2; the
residual radii are |N(0, 1)| times rho, where rho controls how much error the operands carry.
"""
import json
import os
import random
import time

import numpy as np

from common import RESULTS
from raf import set_box, RAFForm, mul_faf2, mul_skalna, mul_exact

RULES = {"faf2": mul_faf2, "skalna": mul_skalna, "exact": mul_exact}
DIMS = (2, 5, 10, 20, 50)
RHOS = (0.0, 0.1, 1.0)
TRIALS = 500


def rand_form(n, rho):
    c = random.gauss(0, 2)
    r = abs(random.gauss(0, 1)) * rho
    return RAFForm([random.gauss(0, 1) for _ in range(n)], c - r, c + r)


if __name__ == "__main__":
    random.seed(2024)
    out = []
    print(f"{'n':>3s} {'rho':>4s} " + " ".join(f"{m + '/exact':>12s}" for m in ("faf2", "skalna"))
          + "  skalna>exact  max skalna/exact   time per product (us): " + " ".join(RULES))
    for n in DIMS:
        set_box([(-1.0, 1.0)] * n)
        for rho in RHOS:
            pairs = [(rand_form(n, rho), rand_form(n, rho)) for _ in range(TRIALS)]
            widths, times = {}, {}
            for name, mul in RULES.items():
                t0 = time.perf_counter()
                res = [mul(v, w) for v, w in pairs]
                times[name] = (time.perf_counter() - t0) / TRIALS * 1e6
                widths[name] = np.array([u.residual_width() for u in res])
            we = widths["exact"]
            ratios = {m: float(np.mean(widths[m] / we)) for m in ("faf2", "skalna")}
            strict = int(np.sum(widths["skalna"] > we * (1 + 1e-9)))
            mx = float(np.max(widths["skalna"] / we))
            out.append({"n": n, "rho": rho, "trials": TRIALS, "ratio_to_exact": ratios,
                        "skalna_strictly_wider": strict, "max_skalna_over_exact": mx, "time_us": times})
            print(f"{n:3d} {rho:4.1f} " + " ".join(f"{ratios[m]:12.3f}" for m in ("faf2", "skalna"))
                  + f"  {strict:5d}/{TRIALS}  {mx:12.3f}        "
                  + " ".join(f"{times[m]:7.1f}" for m in RULES))
    with open(os.path.join(RESULTS, "random_products.json"), "w") as fh:
        json.dump(out, fh, indent=1)
