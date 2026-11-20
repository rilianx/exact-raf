#include "ibex.h"
#include "ibex_ExactMult_interval.h"
#include "ibex_ExactMult_fast.h"
#include <random>
#include <chrono>
#include <cstdio>
using namespace ibex;
int main(int argc, char** argv) {
  std::mt19937_64 R(argc > 2 ? atoi(argv[2]) : 7); std::uniform_real_distribution<double> U(-1, 1);
  const int mode = atoi(argv[1]);
  if (mode == 0) { // volcar casos + resultados para verificacion exacta en Python
    for (int c = 0; c < (argc > 3 ? atoi(argv[3]) : 3000); c++) {
      int n = 1 + R() % 12; std::vector<double> a(n), b(n);
      double sc = std::pow(10.0, (int)(R() % 13) - 6);
      for (int i = 0; i < n; i++) { a[i] = U(R) * sc; b[i] = U(R); }
      if (c % 3 == 0) for (int i = 1; i < n; i++) { a[i] = a[0] * (1 + 1e-15 * (R() % 5)); b[i] = b[0] * (1 - 1e-15 * (R() % 3)); } // casi paralelos
      if (c % 5 == 0) for (int i = 0; i < n; i++) { if (R()%2) a[i] = 0; }
      double x0 = U(R) * 3, y0 = U(R) * 3, ex = std::fabs(U(R)) * (c % 4 ? 1 : 0), ey = std::fabs(U(R));
      for (int e = 0; e < 2; e++) {
        Interval F = exact_mult_fast::raf_mult_enclosure_fast(a.data(), b.data(), n, x0, ex, y0, ey, e);
        Interval I = raf_mult_enclosure(a, b, x0, ex, y0, ey, e);
        printf("%d %d %a %a %a %a", e, n, x0, ex, y0, ey);
        for (int i = 0; i < n; i++) printf(" %a %a", a[i], b[i]);
        printf(" | %a %a %a %a\n", F.lb(), F.ub(), I.lb(), I.ub());
      }
    }
    return 0;
  }
  // tiempos
  printf("n,errores,sh_interval_us,sh_fast_us,exact_interval_us,exact_fast_us\n");
  for (int n : {2, 5, 10, 20, 50, 100}) for (int we = 0; we < 2; we++) {
    const int P = 2000; std::vector<std::vector<double>> A(P, std::vector<double>(n)), B = A; std::vector<double> X0(P), Y0(P), EX(P), EY(P);
    for (int p = 0; p < P; p++) { for (int i = 0; i < n; i++) { A[p][i] = U(R); B[p][i] = U(R); } X0[p] = U(R)*3; Y0[p] = U(R)*3; EX[p] = we ? std::fabs(U(R)) : 0; EY[p] = we ? std::fabs(U(R)) : 0; }
    double t[4]; double sink = 0; int reps = std::max(1, 20000 / (n * 10));
    for (int m = 0; m < 4; m++) {
      auto t0 = std::chrono::steady_clock::now();
      for (int r = 0; r < reps; r++) for (int p = 0; p < P; p++) {
        bool ex = m >= 2; Interval I;
        if (m % 2 == 0) I = raf_mult_enclosure(A[p], B[p], X0[p], EX[p], Y0[p], EY[p], ex);
        else I = exact_mult_fast::raf_mult_enclosure_fast(A[p].data(), B[p].data(), n, X0[p], EX[p], Y0[p], EY[p], ex);
        sink += I.ub();
      }
      t[m] = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count() / (reps * P);
    }
    printf("%d,%s,%.3f,%.3f,%.3f,%.3f\n", n, we ? "si" : "no", t[0], t[1], t[2], t[3]);
    if (sink == 42) puts("");
  }
}
