// Prueba del cuadrado exacto (parche 0020): expresiones con cuadrados de subexpresiones con error,
// evaluadas en modo SH sobre cajas aleatorias; contencion en puntos muestreados (intervalos) y anchos.
#include "ibex.h"
#include "ibex_Affine.h"
#include "ibex_KNomial.h"
#include "ibex_ExactMult.h"
#include <random>
#include <cstdio>
using namespace ibex; using namespace std;
static mt19937_64 rng(7);
static double U(double a, double b) { return uniform_real_distribution<double>(a, b)(rng); }
static string num(double v) { char b[64]; snprintf(b, sizeof b, "%.17g", v); return v < 0 ? "(" + string(b) + ")" : string(b); }
static string L(int n) { string s = num(U(-1, 1)); for (int i = 1; i <= n; i++) s += "+" + num(U(-1, 1)) + "*x(" + to_string(i) + ")"; return "(" + s + ")"; }
int main(int argc, char** argv) {
  int N = argc > 1 ? atoi(argv[1]) : 3000; long viol = 0, cases = 0; double sumw = 0;
  for (int c = 0; c < N; c++) {
    int n = 1 + c % 3, fam = c % 5; string e;
    switch (fam) {
      case 0: e = "(" + L(n) + "^2)^2"; break;
      case 1: e = "(" + L(n) + "*" + L(n) + "+" + num(U(-2, 2)) + ")^2"; break;
      case 2: e = "(sin(" + L(n) + ")+" + num(U(-1, 1)) + ")^2"; break;
      case 3: e = "(exp(" + L(n) + ")+" + num(U(-2, 2)) + ")^2"; break;
      default: e = "(" + L(n) + "^2+" + L(n) + ")^2*" + L(n); break;
    }
    Function f(n == 1 ? "x[1]" : (n == 2 ? "x[2]" : "x[3]"), e.c_str());
    IntervalVector B(n); double rs = pow(10.0, U(-3, 0.5));
    for (int i = 0; i < n; i++) { double m = U(-2, 2), r = rs * U(0.2, 1); B[i] = Interval(m - r, m + r); }
    affine_mode() = AFFINE_SH;
    Affine2Eval ev(f);
    pair<Domain*, TemplateDomain<AffineMain<AF_fAF2> >*> res = ev.eval(B, AffineVarMainVector<AF_fAF2>(B));
    const AffineMain<AF_fAF2>& af = res.second->i();
    affine_mode() = AFFINE_AF2;
    if (!af.is_active()) continue;
    cases++; sumw += af.err();
    for (int s = 0; s < 40; s++) {
      Vector p(n); for (int i = 0; i < n; i++) p[i] = (s < 20) ? (U(0, 1) < .5 ? B[i].lb() : B[i].ub()) : B[i].lb() + U(0, 1) * B[i].diam();
      Interval F = f.eval(IntervalVector(p)), E(af.mid());
      for (int i = 0; i < n && i < af.noise_count(); i++) E += Interval(af.noise(i)) * ((Interval(p[i]) - B[i].mid()) / Interval(B[i].rad()));
      E += Interval(-af.err(), af.err());
      if (!E.intersects(F)) { viol++; break; }
    }
    printf("W %d %.17g\n", c, af.err());
  }
  ExactMultStats& st = exact_mult_stats();
  fprintf(stderr, "casos=%ld violaciones=%ld cuadrados_exactos=%ld mas_angostos=%ld suma_err=%.6g\n", cases, viol, st.sq_exact, st.sq_tighter, sumw);
}
