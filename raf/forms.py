"""Reduced affine forms written in the variables.

A form is f(x) = sum_i a_i x_i + r, with r in [lo, hi], over the box set with
raf.set_box: with x_i = mid_i + rad_i eps_i, it is the reduced affine form whose noise symbols
are the variables and whose error term is the residual interval. Arithmetic is in floating
point, without outward rounding. The product is selected with set_product: 'faf2' (S),
'skalna' (SH) or 'exact' (E).
"""
from .box import get_box


class RAFForm:
    """f(x) = sum(slopes[i] * x_i) + [lo, hi]."""
    __slots__ = ("slopes", "lo", "hi")

    def __init__(self, slopes, lo, hi=None):
        self.slopes = slopes if isinstance(slopes, list) else list(slopes)
        self.lo = float(lo)
        self.hi = float(lo if hi is None else hi)

    # ---- ranges -------------------------------------------------------------
    def linear_range(self):
        lo = hi = 0.0
        for s, (a, b) in zip(self.slopes, get_box()):
            if s > 0:
                lo += s * a
                hi += s * b
            elif s < 0:
                lo += s * b
                hi += s * a
        return lo, hi

    def get_range(self):
        """Enclosure of the range of the form over the box."""
        lo, hi = self.linear_range()
        return lo + self.lo, hi + self.hi

    def residual_width(self):
        """Vertical distance between the two parallel envelopes."""
        return self.hi - self.lo

    def __call__(self, x):
        """Lower and upper envelopes evaluated at the point x."""
        lin = sum(s * xi for s, xi in zip(self.slopes, x))
        return lin + self.lo, lin + self.hi

    # ---- affine operations (exact) ------------------------------------------
    def __add__(self, other):
        if isinstance(other, RAFForm):
            return RAFForm([s + o for s, o in zip(self.slopes, other.slopes)],
                           self.lo + other.lo, self.hi + other.hi)
        if isinstance(other, (int, float)):
            return RAFForm(self.slopes, self.lo + other, self.hi + other)
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return RAFForm([-s for s in self.slopes], -self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def scale(self, c):
        c = float(c)
        lo, hi = (c * self.lo, c * self.hi) if c >= 0 else (c * self.hi, c * self.lo)
        return RAFForm([s * c for s in self.slopes], lo, hi)

    # ---- nonlinear operations -----------------------------------------------
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return self.scale(other)
        if isinstance(other, RAFForm):
            return PRODUCTS[_PRODUCT](self, other)
        return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return self.scale(1.0 / other)
        raise NotImplementedError("division by a form is not implemented")

    def __pow__(self, k):
        if _SQUARE == "product" and isinstance(k, int) and k >= 1:
            out = self                     # x^k as repeated multiplication with the selected rule
            for _ in range(k - 1):
                out = out * self
            return out
        if k == 2:
            return chebyshev_square(self)
        from .unary import power
        return power(self, k)

    def __repr__(self):
        terms = " + ".join(f"{a:g}*x{i + 1}" for i, a in enumerate(self.slopes) if a != 0) or "0"
        return f"RAFForm({terms} + [{self.lo:g}, {self.hi:g}])"


def raf_variables(n=None):
    """Forms of the variables x_1..x_n of the current box."""
    n = len(get_box()) if n is None else n
    out = []
    for i in range(n):
        s = [0.0] * n
        s[i] = 1.0
        out.append(RAFForm(s, 0.0, 0.0))
    return out


def chebyshev_square(form):
    """Square of a form: Chebyshev slope a+b over its range [a, b] and residual
    [-(a+b)^2/4, -ab] (Chebyshev linearization)."""
    a, b = form.get_range()
    if b - a < 1e-15:
        v = a * a
        return RAFForm([0.0] * len(form.slopes), v, v)
    m = a + b
    rlo, rhi = -m * m / 4.0, -a * b
    elo, ehi = (m * form.lo, m * form.hi) if m >= 0 else (m * form.hi, m * form.lo)
    return RAFForm([s * m for s in form.slopes], rlo + elo, rhi + ehi)


def _noise(v):
    box = get_box()
    m = [0.5 * (p + q) for p, q in box]
    rho = [0.5 * (q - p) for p, q in box]
    v0 = sum(a * mi for a, mi in zip(v.slopes, m)) + 0.5 * (v.lo + v.hi)
    return v0, [a * r for a, r in zip(v.slopes, rho)], 0.5 * (v.hi - v.lo), m


def _assemble(v, w, nlo, nhi):
    """Form v0 w0 + sum_i (v0 b_i + w0 a_i)(x_i - m_i) + [nlo, nhi]."""
    v0, _, _, m = _noise(v)
    w0, _, _, _ = _noise(w)
    slopes = [v0 * b + w0 * a for a, b in zip(v.slopes, w.slopes)]
    c = v0 * w0 - sum(s * mi for s, mi in zip(slopes, m))
    return RAFForm(slopes, c + nlo, c + nhi)


def _sorted_generators(gens):
    """Nonzero generators turned to the half-plane u > 0 (or u = 0, v > 0), sorted by angle."""
    import math
    g = []
    for gu, gv in gens:
        if gu < 0 or (gu == 0 and gv < 0):
            gu, gv = -gu, -gv
        if gu != 0 or gv != 0:
            g.append((gu, gv))
    g.sort(key=lambda t: math.atan2(t[1], t[0]))
    return g, sum(t[0] for t in g), sum(t[1] for t in g)


def _uv_range_sorted(sg, c0=0.0, c1=0.0):
    """Range of (u + c0)(v + c1) over the zonotope of the sorted generators sg (see _uv_range). O(K)."""
    g, su, sv = sg
    if not g:
        return c0 * c1, c0 * c1
    pu, pv = c0 - su, c1 - sv
    lo = hi = pu * pv
    for sign in (2.0, -2.0):
        for gu, gv in g:
            du, dv = sign * gu, sign * gv
            qu, qv = pu + du, pv + dv
            val = qu * qv
            if val < lo:
                lo = val
            elif val > hi:
                hi = val
            a = du * dv                                 # (pu + t du)(pv + t dv) = pu pv + b t + a t^2
            if a != 0:
                t = -(pu * dv + pv * du) / (2.0 * a)
                if 0.0 < t < 1.0:
                    val = (pu + t * du) * (pv + t * dv)
                    lo, hi = min(lo, val), max(hi, val)
            pu, pv = qu, qv
    return lo, hi


def _uv_range(gens, c0=0.0, c1=0.0):
    """Exact range of (u + c0)(v + c1) for (u, v) in the zonotope sum_k e_k g_k, e in [-1, 1]^K.

    The only critical point of (u + c0)(v + c1) is a saddle, so the extrema lie on the boundary
    of the zonotope, a convex polygon whose edges are the generators sorted by angle
    (Skalna and Hladik 2017, Prop. 4). O(K log K)."""
    return _uv_range_sorted(_sorted_generators(gens), c0, c1)


def nonlinear_enclosure(rule, x0, a, xr, y0, b, yr):
    """Interval [lo, hi] enclosing the nonlinear part N of the product of the reduced affine forms
    x = x0 + sum_i a_i e_i + xr ex and y = y0 + sum_i b_i e_i + yr ey (Eq. (3) of the paper):
        N = x0 yr ey + y0 xr ex + (a.e + xr ex)(b.e + yr ey),
    with rule 'faf2' (S), 'skalna' (SH) or 'exact' (E). The product is then
        x0 y0 + sum_i (x0 b_i + y0 a_i) e_i + [lo, hi],
    that is, center x0 y0 + (lo + hi)/2 and error term (hi - lo)/2. This function is self-contained
    (plain lists of floats) and is the reference to port the multiplications to another library."""
    a, b = list(a), list(b)
    if rule == "faf2":          # standard multiplication (fAF2 in Ibex: AF2 with re-centering)
        sx, sy = sum(abs(t) for t in a), sum(abs(t) for t in b)
        c = 0.5 * sum(ai * bi for ai, bi in zip(a, b))
        r = abs(y0) * xr + abs(x0) * yr + (sx + xr) * (sy + yr) - 0.5 * sum(abs(ai * bi) for ai, bi in zip(a, b))
        return c - r, c + r
    if rule == "skalna":        # Skalna and Hladik (2017, Eq. 10): exact range of R, plus |y0| xr + |x0| yr
        rlo, rhi = _uv_range(list(zip(a, b)) + [(xr, 0.0), (0.0, yr)])
        e = abs(y0) * xr + abs(x0) * yr
        return rlo - e, rhi + e
    if rule == "exact":         # exact range of N (Algorithm 1): N is affine in ex and in ey, so its extrema
        sg = _sorted_generators(zip(a, b))      # are attained at ex, ey in {-1, 1}; generators sorted once
        lo, hi = float("inf"), float("-inf")
        for s in (-1.0, 1.0):                   # opposite corners give the same range (Lemma 3): two walks
            ulo, uhi = _uv_range_sorted(sg, xr, s * yr)
            k = abs(x0 * s * yr + y0 * xr)
            lo, hi = min(lo, ulo - k), max(hi, uhi + k)
        return lo, hi
    raise ValueError(f"unknown rule {rule!r}")


def _mul(rule, v, w):
    v0, a, vr, _ = _noise(v)
    w0, b, wr, _ = _noise(w)
    lo, hi = nonlinear_enclosure(rule, v0, a, vr, w0, b, wr)
    return _assemble(v, w, lo, hi)


def mul_faf2(v, w):
    """Standard multiplication S (Skalna and Hladik 2017, Eq. 5), the rule of the fAF2 model of Ibex."""
    return _mul("faf2", v, w)


def mul_skalna(v, w):
    """Multiplication SH of Skalna and Hladik (2017, Eq. 10)."""
    return _mul("skalna", v, w)


def mul_exact(v, w):
    """Exact multiplication E (Algorithm 1 of the paper): O(n log n)."""
    return _mul("exact", v, w)


PRODUCTS = {"faf2": mul_faf2, "skalna": mul_skalna, "exact": mul_exact}
_PRODUCT = "exact"


def set_product(name):
    """Select the multiplication rule used by RAFForm.__mul__ (one of PRODUCTS)."""
    global _PRODUCT
    if name not in PRODUCTS:
        raise ValueError(f"unknown product {name!r}; choose one of {sorted(PRODUCTS)}")
    _PRODUCT = name


def get_product():
    return _PRODUCT


_SQUARE = "chebyshev"


def set_square(mode):
    """Integer powers: 'chebyshev' (univariate envelope, default) or 'product' (repeated products)."""
    global _SQUARE
    if mode not in ("chebyshev", "product"):
        raise ValueError(mode)
    _SQUARE = mode
