"""Zero-centered affine arithmetic (ZCAA).

A form is f(x) = sum_i a_i x_i + r, with r in [lo, hi], over the box set with
zcaa.set_box. The residual includes the intercept and is not assumed to be
centered. Arithmetic is in floating point, without outward rounding.
"""
from .box import get_box


class ZCAAForm:
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
        if isinstance(other, ZCAAForm):
            return ZCAAForm([s + o for s, o in zip(self.slopes, other.slopes)],
                           self.lo + other.lo, self.hi + other.hi)
        if isinstance(other, (int, float)):
            return ZCAAForm(self.slopes, self.lo + other, self.hi + other)
        return NotImplemented

    def __radd__(self, other):
        return self + other

    def __neg__(self):
        return ZCAAForm([-s for s in self.slopes], -self.hi, -self.lo)

    def __sub__(self, other):
        return self + (-other)

    def __rsub__(self, other):
        return (-self) + other

    def scale(self, c):
        c = float(c)
        lo, hi = (c * self.lo, c * self.hi) if c >= 0 else (c * self.hi, c * self.lo)
        return ZCAAForm([s * c for s in self.slopes], lo, hi)

    # ---- nonlinear operations -----------------------------------------------
    def __mul__(self, other):
        if isinstance(other, (int, float)):
            return self.scale(other)
        if isinstance(other, ZCAAForm):
            return PRODUCTS[_PRODUCT](self, other)
        return NotImplemented

    def __rmul__(self, other):
        return self * other

    def __truediv__(self, other):
        if isinstance(other, (int, float)):
            return self.scale(1.0 / other)
        raise NotImplementedError("division by a form is not implemented")

    def __pow__(self, k):
        if k == 2:
            return pow_zcaa(self)
        from .unary import power
        return power(self, k)

    def __repr__(self):
        terms = " + ".join(f"{a:g}*x{i + 1}" for i, a in enumerate(self.slopes) if a != 0) or "0"
        return f"ZCAAForm({terms} + [{self.lo:g}, {self.hi:g}])"


def zcaa_variables(n=None):
    """ZCAA forms of the variables x_1..x_n of the current box."""
    n = len(get_box()) if n is None else n
    out = []
    for i in range(n):
        s = [0.0] * n
        s[i] = 1.0
        out.append(ZCAAForm(s, 0.0, 0.0))
    return out


def pow_zcaa(form):
    """Square of a form: Chebyshev slope a+b over its range [a, b] and residual
    [-(a+b)^2/4, -ab] (Section 3.2 of the paper)."""
    a, b = form.get_range()
    if b - a < 1e-15:
        v = a * a
        return ZCAAForm([0.0] * len(form.slopes), v, v)
    m = a + b
    rlo, rhi = -m * m / 4.0, -a * b
    elo, ehi = (m * form.lo, m * form.hi) if m >= 0 else (m * form.hi, m * form.lo)
    return ZCAAForm([s * m for s in form.slopes], rlo + elo, rhi + ehi)


def mul_zcaa(v, w):
    """Product of two ZCAA forms (Section 3.4, Eq. 17 of the paper).

    The midpoints of the residuals are absorbed into the linear part, every
    quadratic and bilinear monomial is enclosed by its own envelope, and the
    residual of the cross terms (Eq. 16) is evaluated at the four vertices.
    """
    box = get_box()
    a, b = v.slopes, w.slopes
    lv, hv, lw, hw = v.lo, v.hi, w.lo, w.hi
    mv, mw = 0.5 * (lv + hv), 0.5 * (lw + hw)
    nz = [i for i in range(len(a)) if a[i] != 0 or b[i] != 0]

    # 1. linear part (midpoints of the residuals absorbed)
    slopes = [mw * ai + mv * bi for ai, bi in zip(a, b)]
    elo = ehi = 0.0

    # 2. quadratic terms a_i b_i x_i^2 (Chebyshev envelope)
    for i in nz:
        c = a[i] * b[i]
        if c == 0:
            continue
        p, q = box[i]
        if q - p < 1e-15:
            elo += c * p * p
            ehi += c * p * p
            continue
        m = p + q
        rlo, rhi = -m * m / 4.0, -p * q
        slopes[i] += c * m
        if c > 0:
            elo += c * rlo
            ehi += c * rhi
        else:
            elo += c * rhi
            ehi += c * rlo

    # 3. bilinear terms (a_i b_j + a_j b_i) x_i x_j (tangent plane at the center)
    for k, i in enumerate(nz):
        pi, qi = box[i]
        for j in nz[k + 1:]:
            c = a[i] * b[j] + a[j] * b[i]
            if c == 0:
                continue
            pj, qj = box[j]
            mi, mj = 0.5 * (pj + qj), 0.5 * (pi + qi)
            r = (qi - pi) * (qj - pj) / 4.0
            shift = mi * mj
            slopes[i] += c * mi
            slopes[j] += c * mj
            if c > 0:
                elo += c * (-r - shift)
                ehi += c * (r - shift)
            else:
                elo += c * (r - shift)
                ehi += c * (-r - shift)

    # 4. residual of the cross terms, Eq. (16), evaluated at the four vertices
    glo, ghi = float("inf"), float("-inf")
    for rv in (lv, hv):
        dv = rv - mv
        for rw in (lw, hw):
            dw = rw - mw
            lo = hi = rv * rw
            for i in nz:
                c = dw * a[i] + dv * b[i]
                if c > 0:
                    lo += c * box[i][0]
                    hi += c * box[i][1]
                elif c < 0:
                    lo += c * box[i][1]
                    hi += c * box[i][0]
            glo = min(glo, lo)
            ghi = max(ghi, hi)

    return ZCAAForm(slopes, elo + glo, ehi + ghi)


# ---- other multiplication rules on the same representation ---------------------
#
# In noise coordinates, x_i = m_i + rho_i e_i with e_i in [-1, 1], and a form
# v = sum_i a_i x_i + [lo, hi] is the reduced affine form
#     v = v0 + sum_i alpha_i e_i + v_r xi,     v0 = a.m + mid([lo, hi]), alpha_i = a_i rho_i,
# with v_r = rad([lo, hi]) and xi in [-1, 1] (Skalna and Hladik, Numer. Algor. 2017, Eq. 1).
# The product of v and w = w0 + sum_i beta_i e_i + w_r eta is
#     v w = v0 w0 + sum_i (v0 beta_i + w0 alpha_i) e_i + N,
#     N = v0 w_r eta + w0 v_r xi + (alpha.e + v_r xi)(beta.e + w_r eta)        (*)
# The rules below share this linear part, which is also the linear part of mul_zcaa,
# and differ only in the interval that encloses N.

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
    return ZCAAForm(slopes, c + nlo, c + nhi)


def _uv_range(gens, c0=0.0, c1=0.0):
    """Exact range of (u + c0)(v + c1) for (u, v) in the zonotope sum_k e_k g_k, e in [-1, 1]^K.

    The only critical point of uv is a saddle, so the extrema lie on the boundary of the
    zonotope, a convex polygon whose edges are the generators sorted by angle
    (Skalna and Hladik 2017, Prop. 4). O(K log K)."""
    g = []
    for gu, gv in gens:
        if gu < 0 or (gu == 0 and gv < 0):
            gu, gv = -gu, -gv
        if gu != 0 or gv != 0:
            g.append((gu, gv))
    if not g:
        return c0 * c1, c0 * c1
    import math
    g.sort(key=lambda t: math.atan2(t[1], t[0]))
    pu, pv = c0 - sum(t[0] for t in g), c1 - sum(t[1] for t in g)
    lo = hi = pu * pv
    for sign in (2.0, -2.0):
        for gu, gv in g:
            du, dv = sign * gu, sign * gv
            qu, qv = pu + du, pv + dv
            val = qu * qv
            lo, hi = min(lo, val), max(hi, val)
            a, b = du * dv, pu * dv + pv * du          # (pu + t du)(pv + t dv) = pu pv + b t + a t^2
            if a != 0:
                t = -b / (2.0 * a)
                if 0.0 < t < 1.0:
                    val = (pu + t * du) * (pv + t * dv)
                    lo, hi = min(lo, val), max(hi, val)
            pu, pv = qu, qv
    return lo, hi


def mul_faf2(v, w):
    """'Standard' multiplication of reduced affine forms (Skalna and Hladik 2017, Eq. 5),
    the rule of the fAF2 model of the Ibex affine plugin."""
    v0, a, vr, _ = _noise(v)
    w0, b, wr, _ = _noise(w)
    sx, sy = sum(abs(t) for t in a), sum(abs(t) for t in b)
    c = 0.5 * sum(ai * bi for ai, bi in zip(a, b))
    r = abs(w0) * vr + abs(v0) * wr + (sx + vr) * (sy + wr) - 0.5 * sum(abs(ai * bi) for ai, bi in zip(a, b))
    return _assemble(v, w, c - r, c + r)


def mul_skalna(v, w):
    """Chebyshev minimum-error multiplication of Skalna and Hladik (2017, Eq. 10): exact range of
    R = (alpha.e + v_r xi)(beta.e + w_r eta), plus |w0| v_r + |v0| w_r for the other terms of (*)."""
    v0, a, vr, _ = _noise(v)
    w0, b, wr, _ = _noise(w)
    rlo, rhi = _uv_range(list(zip(a, b)) + [(vr, 0.0), (0.0, wr)])
    e = abs(w0) * vr + abs(v0) * wr
    return _assemble(v, w, rlo - e, rhi + e)


def mul_exact(v, w):
    """Exact range of N in (*). N is affine in xi and in eta, so its extrema are attained with
    xi, eta in {-1, 1} (the four vertices of mul_zcaa); for each vertex, N is a product of two
    shifted linear forms in e, whose exact range follows from _uv_range. O(n log n)."""
    v0, a, vr, _ = _noise(v)
    w0, b, wr, _ = _noise(w)
    gens = list(zip(a, b))
    nlo, nhi = float("inf"), float("-inf")
    for xi in (-1.0, 1.0):
        for eta in (-1.0, 1.0):
            lo, hi = _uv_range(gens, vr * xi, wr * eta)
            k = v0 * wr * eta + w0 * vr * xi
            nlo, nhi = min(nlo, lo + k), max(nhi, hi + k)
    return _assemble(v, w, nlo, nhi)


PRODUCTS = {"zcaa": mul_zcaa, "faf2": mul_faf2, "skalna": mul_skalna, "exact": mul_exact}
_PRODUCT = "zcaa"


def set_product(name):
    """Select the multiplication rule used by ZCAAForm.__mul__ (one of PRODUCTS)."""
    global _PRODUCT
    if name not in PRODUCTS:
        raise ValueError(f"unknown product {name!r}; choose one of {sorted(PRODUCTS)}")
    _PRODUCT = name


def get_product():
    return _PRODUCT
