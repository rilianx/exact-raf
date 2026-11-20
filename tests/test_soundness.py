"""Sampling checks of the enclosures and regression tests.

    python -m pytest tests          # or: python tests/test_soundness.py

For random boxes, every test expression is evaluated at random points of the box,
and the value must lie between the lower and upper envelopes of ZCAA and AF2.
This is a numerical check, not a proof (the arithmetic has no outward rounding).
"""
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from zcaa import (set_box, zcaa_variables, af2_variables, tm_variables, AF2, ZCAAForm, mul_zcaa,
                  mul_faf2, mul_skalna, mul_exact, set_product, PRODUCTS)   # noqa: E402
from benchmarks.functions import ALL, SIZES                                     # noqa: E402

BOXES, POINTS = 20, 200


def _box(f, size):
    wmin, wmax = SIZES[size]
    lo, hi = f.domain
    out = []
    for _ in range(3):
        w = min(random.uniform(wmin, wmax), hi - lo)
        s = random.uniform(lo, hi - w)
        out.append((s, s + w))
    return out


def test_enclosures():
    random.seed(0)
    for f in ALL:
        for size in f.sizes:
            for _ in range(BOXES):
                box = _box(f, size)
                set_box(box)
                r = f(*zcaa_variables(3))
                a = f(*af2_variables(box))
                t = f(*tm_variables(3))
                for _ in range(POINTS):
                    x = [random.uniform(lo, hi) for lo, hi in box]
                    v = f(*x)
                    tol = 1e-7 * max(1.0, abs(v))
                    lo, hi = r(x)
                    assert lo - tol <= v <= hi + tol, (f.name, "ZCAA", box, x)
                    lo, hi = a.envelopes(x, box)
                    assert lo - tol <= v <= hi + tol, (f.name, "AF2", box, x)
                    lo, hi = t(x)
                    assert lo - tol <= v <= hi + tol, (f.name, "TM", box, x)


def test_zcaa_never_looser_than_af2_and_tm():
    random.seed(1)
    for f in ALL:
        for size in f.sizes:
            for _ in range(BOXES):
                box = _box(f, size)
                set_box(box)
                wr = f(*zcaa_variables(3)).residual_width()
                wa = f(*af2_variables(box)).residual_width()
                wt = f(*tm_variables(3)).residual_width()
                assert wr <= wa * (1 + 1e-7) + 1e-9, (f.name, size, wr, wa)
                assert wr <= wt * (1 + 1e-7) + 1e-9, (f.name, size, wr, wt)


def test_af2_independent_error_symbols():
    # (x1^2 - 0.5)(0.5 - x2^2) over [-1, 1]^2 ranges over [-0.25, 0.25]; the product of
    # the two e(+-) terms must be enclosed in [-1, 1], not in [0, 1].
    x1 = AF2(0.0, {1: 1.0})
    x2 = AF2(0.0, {2: 1.0})
    lo, hi = ((x1 ** 2 - 0.5) * (0.5 - x2 ** 2)).get_range()
    assert lo <= -0.25 + 1e-12 and hi >= 0.25 - 1e-12


def test_af2_negative_scaling_swaps_one_sided_terms():
    x = AF2(1.0, {1: 1.0}, s=0.0, p=0.5, n=0.0)        # range [0, 2.5]
    lo, hi = (-2 * x).get_range()
    assert abs(lo + 5.0) < 1e-12 and abs(hi) < 1e-12


def test_af2_add_constant_is_not_in_place():
    x = AF2(1.0, {1: 1.0})
    _ = x + 3.0
    assert x.base == 1.0


def test_lemma_example():
    # Example of Section 3.4: v = x1 + 2 x2 + [-1, 1], w = 2 x1 + 3 x2 + [0, 2], x in [0, 2]^2
    set_box([(0, 2), (0, 2)])
    u = mul_zcaa(ZCAAForm([1.0, 2.0], -1, 1), ZCAAForm([2.0, 3.0], 0, 2))
    assert [round(s, 12) for s in u.slopes] == [12.0, 21.0]
    assert abs(u.lo + 38) < 1e-12 and abs(u.hi - 18) < 1e-12


def test_product_rules_enclose():
    """Every multiplication rule gives valid envelopes on the test expressions."""
    random.seed(3)
    try:
        for rule in PRODUCTS:
            set_product(rule)
            for f in ALL:
                for size in f.sizes:
                    for _ in range(5):
                        box = _box(f, size)
                        set_box(box)
                        r = f(*zcaa_variables(3))
                        for _ in range(POINTS):
                            x = [random.uniform(lo, hi) for lo, hi in box]
                            v = f(*x)
                            tol = 1e-7 * max(1.0, abs(v))
                            lo, hi = r(x)
                            assert lo - tol <= v <= hi + tol, (f.name, rule, box, x)
    finally:
        set_product("zcaa")


def test_product_rules_inclusions():
    """Single product: exact in skalna and in zcaa, zcaa and skalna in faf2, same linear part."""
    random.seed(4)
    for _ in range(2000):
        n = random.randint(1, 6)
        box = []
        for _ in range(n):
            lo = random.uniform(-5, 5)
            box.append((lo, lo + random.uniform(0, 4)))
        set_box(box)

        def rand_form():
            c, w = random.uniform(-3, 3), random.uniform(0, 2) * random.random()
            return ZCAAForm([random.gauss(0, 1) for _ in range(n)], c - w * random.random(), c + w * random.random())

        v, w = rand_form(), rand_form()
        e, s, z, f = mul_exact(v, w), mul_skalna(v, w), mul_zcaa(v, w), mul_faf2(v, w)
        for p, q in ((e, s), (e, z), (z, f), (s, f)):
            assert all(abs(a - b) <= 1e-9 * (1 + abs(a)) for a, b in zip(p.slopes, q.slopes))
            tol = 1e-9 * (1 + abs(q.lo) + abs(q.hi))
            assert q.lo - tol <= p.lo and p.hi <= q.hi + tol


def test_exact_product_mixed_signs():
    # (x1 + x2)(x1 - x2) = x1^2 - x2^2 over [-1, 1]^2: exact residual [-1, 1], AF2 and fAF2 width 6
    set_box([(-1, 1), (-1, 1)])
    v, w = ZCAAForm([1.0, 1.0], 0, 0), ZCAAForm([1.0, -1.0], 0, 0)
    e, f = mul_exact(v, w), mul_faf2(v, w)
    assert abs(e.lo + 1) < 1e-12 and abs(e.hi - 1) < 1e-12
    assert abs(f.residual_width() - 6) < 1e-12


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_") and callable(fn):
            fn()
            print("ok", name)
