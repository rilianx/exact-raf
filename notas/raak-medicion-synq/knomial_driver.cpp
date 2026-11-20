// Driver de la medicion synq: aplica la reescritura k-nomial de raak/affine/ibex_KNomial.h
// (JCGeraldo/ibex-lib, rama raak-experimentos) a una Q simetrica cualquiera o a un producto.
// Entrada (stdin), una linea por trabajo:
//   Q n k  q11 q12 ... qnn        -> cota [L,U] de x^T Q x en [-1,1]^n (Q completa)
//   P n k  a1..an b1..bn          -> knomial_quad_bound del producto (a^T x)(b^T x)
// Salida: L U nbloques tam_max_bloque tiempo_us (tiempo solo para "Q")
#include <chrono>
#include <cstdio>
#include <iostream>
#include <string>
#include "ibex_KNomial.h"
using namespace ibex;

// Misma cota que knomial_quad_bound en modo 1, pero con Q general:
// [L,U] = sum_b c_b [0, s_b^2] + sum_i R_ii [0,1] + sum_{i<j} 2 R_ij [-1,1],  R = Q - sum_b c_b u_b u_b^T.
static void bound_full(const std::vector<double>& Q0, int n, int k, double& L, double& U, int& nb, int& kmax) {
  std::vector<double> Q = Q0, w(n, 2.0);
  std::vector<KNomialBlock> blocks;
  if (k >= 2) knomial_rewrite(Q, w, n, k, blocks);
  nb = blocks.size(); kmax = 0;
  const int old = std::fegetround(); std::fesetround(FE_UPWARD);
  double u = 0, negl = 0;
  for (auto& b : blocks) {
    double s = 0; int sz = 0;
    for (int i = 0; i < n; i++) if (b.u[i] != 0) { s += std::fabs(b.u[i]); sz++; }
    kmax = std::max(kmax, sz);
    if (b.c > 0) u += s * s; else negl += s * s;
  }
  for (int i = 0; i < n; i++)
    for (int j = i; j < n; j++) {
      double q = 0.5 * (Q0[i * n + j] + Q0[j * n + i]);
      double tu = 0, nt = 0;
      for (auto& b : blocks) if (b.u[i] != 0 && b.u[j] != 0) { double si = b.c * b.u[i]; tu += si * b.u[j]; nt += (-si) * b.u[j]; }
      double rup = q + nt, rlo = -(tu - q);
      if (i == j) { if (rup > 0) u += rup; if (rlo < 0) negl += -rlo; }
      else { double m2 = 2.0 * std::max(std::fabs(rlo), std::fabs(rup)); u += m2; negl += m2; }
    }
  std::fesetround(old);
  L = -negl; U = u;
}

int main() {
  knomial_bound_mode() = 1;
  std::string tag;
  while (std::cin >> tag) {
    int n, k; std::cin >> n >> k;
    double L, U; int nb = 0, kmax = 0;
    double us = 0;
    if (tag == "Q") {
      std::vector<double> Q(n * n);
      for (auto& v : Q) std::cin >> v;
      auto t0 = std::chrono::steady_clock::now();
      bound_full(Q, n, k, L, U, nb, kmax);
      us = std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count();
    } else {
      std::vector<double> a(n), b(n);
      for (auto& v : a) std::cin >> v;
      for (auto& v : b) std::cin >> v;
      double sh, rad;
      knomial_max_k() = k;
      nb = knomial_quad_bound(a, b, n, k, sh, rad);
      L = sh - rad; U = sh + rad;
    }
    std::printf("%.17g %.17g %d %d %.3f\n", L, U, nb, kmax, us);
  }
}
