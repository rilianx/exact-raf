"""A simple interval branch-and-bound on the test functions of Skalna and Hladik (2017, Sect. 5.3),
with each multiplication rule (integer powers as repeated products, as in that paper).

    python experiments/sh_bnb.py      # writes results/sh_bnb.json

Best-first search; the lower bound of a box is the lower end of the enclosure of f computed
with the rule; the upper bound is the value of f at the midpoint. A box is discarded when its
lower bound exceeds the best upper bound minus EPS, and is accepted (not split) when
f(mid) - lower bound < EPS, which mirrors the stopping rule of Skalna and Hladik (they use
EPS = 1e-12 in C++; we use a larger EPS because of the cost of Python). The box is bisected
along its widest side.
"""
import heapq
import json
import os
import time

from common import RESULTS
from zcaa import set_box, zcaa_variables, set_product, set_square
from sh_benchmark import FUNCTIONS, evaluate

RULES = ["faf2", "zcaa", "skalna", "exact"]
EPS = 1e-4
MAX_ITERS = 200000


def bnb(expr, n, box):
    def lower(b):
        set_box(b)
        return evaluate(expr, n, zcaa_variables(n)).get_range()[0]

    def mid_value(b):
        return evaluate(expr, n, [0.5 * (lo + hi) for lo, hi in b])

    ub = mid_value(box)
    heap = [(lower(box), 0, box)]
    count, iters, accepted = 1, 0, 0
    while heap and iters < MAX_ITERS:
        lb, _, b = heapq.heappop(heap)
        if lb > ub - EPS:
            continue
        iters += 1
        fm = mid_value(b)
        ub = min(ub, fm)
        if fm - lb < EPS:
            accepted += 1
            continue
        k = max(range(n), key=lambda i: b[i][1] - b[i][0])
        m = 0.5 * (b[k][0] + b[k][1])
        for half in ((b[k][0], m), (m, b[k][1])):
            c = list(b)
            c[k] = half
            l2 = lower(c)
            if l2 <= ub - EPS:
                heapq.heappush(heap, (l2, count, c))
                count += 1
    return {"iterations": iters, "boxes_left": len(heap) + accepted, "best": ub,
            "finished": not heap or iters < MAX_ITERS}


if __name__ == "__main__":
    set_square("product")
    out = []
    print(f"{'function':22s} " + " ".join(f"{r:>16s}" for r in RULES) + "   (iterations / time s)")
    for name, n, (a, b), expr in FUNCTIONS:
        rec = {"function": name}
        for rule in RULES:
            set_product(rule)
            t0 = time.perf_counter()
            res = bnb(expr, n, [(float(a), float(b))] * n)
            res["time_s"] = time.perf_counter() - t0
            rec[rule] = res
        out.append(rec)
        print(f"{name:22s} " + " ".join(f"{rec[r]['iterations']:9d}/{rec[r]['time_s']:6.1f}" for r in RULES), flush=True)
        with open(os.path.join(RESULTS, "sh_bnb.json"), "w") as fh:
            json.dump(out, fh, indent=1)
    set_product("zcaa")
    set_square("chebyshev")
