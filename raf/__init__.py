"""Reduced affine forms written in the variables, with the multiplications S (faf2), SH (skalna)
and E (exact) of the paper "Exact multiplication of reduced affine forms"."""
from .box import set_box, get_box
from .forms import (RAFForm, raf_variables, chebyshev_square, nonlinear_enclosure, mul_faf2, mul_skalna, mul_exact,
                    set_product, get_product, set_square, PRODUCTS)
from .unary import sin, cos, exp, log, sqrt, sinh, power

__all__ = ["set_box", "get_box", "RAFForm", "raf_variables", "chebyshev_square", "nonlinear_enclosure", "mul_faf2", "mul_skalna",
           "mul_exact", "set_product", "get_product", "set_square", "PRODUCTS",
           "sin", "cos", "exp", "log", "sqrt", "sinh", "power"]
