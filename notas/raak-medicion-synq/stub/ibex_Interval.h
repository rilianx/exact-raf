// Interval minimo, solo para compilar ibex_KNomial.h fuera de Ibex. El driver usa el modo 1
// (dobles con redondeo dirigido), asi que este tipo no se ejecuta.
#pragma once
#include <algorithm>
namespace ibex {
struct Interval {
  double l, u;
  Interval(double a = 0) : l(a), u(a) {}
  Interval(double a, double b) : l(a), u(b) {}
  double ub() const { return u; }
  double mid() const { return 0.5 * (l + u); }
  double rad() const { return 0.5 * (u - l); }
  Interval& operator+=(const Interval& o) { l += o.l; u += o.u; return *this; }
  Interval& operator-=(const Interval& o) { l -= o.u; u -= o.l; return *this; }
};
inline Interval operator+(Interval a, const Interval& b) { return a += b; }
inline Interval operator*(const Interval& a, const Interval& b) {
  double p[4] = {a.l * b.l, a.l * b.u, a.u * b.l, a.u * b.u};
  return Interval(*std::min_element(p, p + 4), *std::max_element(p, p + 4));
}
inline Interval operator*(double a, const Interval& b) { return Interval(a) * b; }
}
