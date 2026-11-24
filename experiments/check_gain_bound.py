"""How often the bound of Corollary 1, wid SH - wid E <= 2 min(|x0| yr, |y0| xr), is attained.
Random products as in Table random: n uniform in {1,...,6}, noise coefficients N(0,1), centers N(0, 2^2),
error radii |N(0,1)|. Run from the repository root:  python experiments/check_gain_bound.py [samples] [seed]"""
import sys, os, math
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import numpy as np
from zcaa.zcaa import _uv_range, _sorted_generators, _uv_range_sorted

def widths(x0, a, xr, y0, b, yr):
    rlo, rhi = _uv_range(list(zip(a, b)) + [(xr, 0.0), (0.0, yr)])
    sh = rhi - rlo + 2 * (abs(y0) * xr + abs(x0) * yr)
    sg = _sorted_generators(zip(a, b)); lo, hi = math.inf, -math.inf
    for s in (-1, 1):
        for t in (-1, 1):
            l, h = _uv_range_sorted(sg, xr * s, yr * t); k = x0 * yr * t + y0 * xr * s
            lo, hi = min(lo, l + k), max(hi, h + k)
    return sh, hi - lo

N = int(sys.argv[1]) if len(sys.argv) > 1 else 100000
rng = np.random.default_rng(int(sys.argv[2]) if len(sys.argv) > 2 else 1)
att, ratios = 0, []
for _ in range(N):
    n = int(rng.integers(1, 7))
    a, b = rng.normal(size=n), rng.normal(size=n)
    x0, y0 = rng.normal(0, 2), rng.normal(0, 2)
    xr, yr = abs(rng.normal()), abs(rng.normal())
    sh, e = widths(x0, list(a), xr, y0, list(b), yr)
    bound = 2 * min(abs(x0) * yr, abs(y0) * xr)
    r = (sh - e) / bound
    ratios.append(r)
    att += abs(r - 1) < 1e-9
ratios = np.array(ratios)
print(f"{N} products: bound attained in {100 * att / N:.1f}%, median ratio {np.median(ratios):.3f}, "
      f"mean ratio {ratios.mean():.3f}, max ratio {ratios.max():.6f}")
for n_lo, n_hi in ((1, 1), (2, 3), (4, 6)):
    pass
