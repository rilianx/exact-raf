"""Sampling checks of the multiplication rules S (faf2), SH (skalna) and E (exact).

    python -m pytest tests          # or: python tests/test_soundness.py

These are numerical checks, not proofs (the Python arithmetic has no outward rounding); the
rigorous implementation is the one in Ibex (ibex/).
"""
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))

from raf import set_box, raf_variables, RAFForm, mul_faf2, mul_skalna, mul_exact, set_product, set_square, PRODUCTS  # noqa: E402
from sh_benchmark import FUNCTIONS, evaluate                                                                         # noqa: E402

POINTS = 300


def _rand_box(n):
    box = []
    for _ in range(n):
        lo = random.uniform(-5, 5)
        box.append((lo, lo + random.uniform(0, 4)))
    return box


def _rand_form(n):
    c, w = random.uniform(-3, 3), random.uniform(0, 2) * random.random()
    return RAFForm([random.gauss(0, 1) for _ in range(n)], c - w * random.random(), c + w * random.random())


def test_products_enclose():
    """For every rule, (v(x) + r) (w(x) + s) lies between the envelopes of the product, for
    random points x of the box and random values r, s of the residuals of the operands."""
    random.seed(1)
    rules = {"faf2": mul_faf2, "skalna": mul_skalna, "exact": mul_exact}
    for _ in range(500):
        n = random.randint(1, 6)
        box = _rand_box(n)
        set_box(box)
        v, w = _rand_form(n), _rand_form(n)
        prods = {k: f(v, w) for k, f in rules.items()}
        for _ in range(50):
            x = [random.uniform(lo, hi) for lo, hi in box]
            r, s = random.uniform(v.lo, v.hi), random.uniform(w.lo, w.hi)
            val = (sum(a * xi for a, xi in zip(v.slopes, x)) + r) * (sum(b * xi for b, xi in zip(w.slopes, x)) + s)
            for k, p in prods.items():
                lo, hi = p(x)
                tol = 1e-9 * (1 + abs(val))
                assert lo - tol <= val <= hi + tol, (k, n, x)


def test_products_inclusions():
    """Single product: E inside SH inside S, with the same linear part (Section 2)."""
    random.seed(4)
    for _ in range(2000):
        n = random.randint(1, 6)
        set_box(_rand_box(n))
        v, w = _rand_form(n), _rand_form(n)
        e, s, f = mul_exact(v, w), mul_skalna(v, w), mul_faf2(v, w)
        for p, q in ((e, s), (s, f)):
            assert all(abs(a - b) <= 1e-9 * (1 + abs(a)) for a, b in zip(p.slopes, q.slopes))
            tol = 1e-9 * (1 + abs(q.lo) + abs(q.hi))
            assert q.lo - tol <= p.lo and p.hi <= q.hi + tol


def test_test_functions_enclose():
    """The test functions of Skalna and Hladik, evaluated with every rule and both treatments of
    the powers, enclose the values of the functions at random points of the box."""
    random.seed(2)
    try:
        for name, n, (a, b), expr in FUNCTIONS:
            box = [(float(a), float(b))] * n
            set_box(box)
            for square in ("product", "chebyshev"):
                set_square(square)
                for rule in PRODUCTS:
                    set_product(rule)
                    form = evaluate(expr, n, raf_variables(n))
                    for _ in range(POINTS):
                        x = [random.uniform(a, b) for _ in range(n)]
                        val = evaluate(expr, n, x)
                        lo, hi = form(x)
                        tol = 1e-7 * max(1.0, abs(val))
                        assert lo - tol <= val <= hi + tol, (name, square, rule, x)
    finally:
        set_product("exact")
        set_square("chebyshev")


def test_exact_product_mixed_signs():
    # (x1 + x2)(x1 - x2) = x1^2 - x2^2 over [-1, 1]^2: exact residual [-1, 1], fAF2 width 6
    set_box([(-1, 1), (-1, 1)])
    v, w = RAFForm([1.0, 1.0], 0, 0), RAFForm([1.0, -1.0], 0, 0)
    e, f = mul_exact(v, w), mul_faf2(v, w)
    assert abs(e.lo + 1) < 1e-12 and abs(e.hi - 1) < 1e-12
    assert abs(f.residual_width() - 6) < 1e-12


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
