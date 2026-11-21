/* Version rapida de raf_mult_enclosure: recorrido en double + cota a priori del redondeo.
 *
 * Cota (eps = 2^-52 cubre cualquier modo de redondeo; m generadores no nulos; c el desplazamiento):
 *   U = |c0| + 3 Su, V = |c1| + 3 Sv (Su = sum |g_u|, Sv = sum |g_v|) acotan |u|, |v| en todo vertice,
 *   suma parcial y punto de arista, verdaderos y calculados. Cada vertice calculado se obtiene con a
 *   lo mas 3m+2 sumas, asi que |p^ - p| <= (3m+2) eps U (idem en v). En cada arista
 *   [p^, p^ + d] (d = +-2g exacto) los valores phi(p^), phi(q^) y P - B^2/(4A) (con tau* en (0,1),
 *   lo que da |B| < 2|A|, |A| <= UV, |B| <= 2UV) se evaluan con error relativo de a lo mas ~12 eps UV,
 *   y un error en la prueba tau* in (0,1) cambia el valor en <= 12 eps UV.
 *   Total por valor: 2 (3m+2) eps UV + O(eps^2) + 24 eps UV + el desplazamiento y kappa.
 *   Se usa M = (8m + 64) eps (UV + |x0| ey + |y0| ex), con holgura, y se redondea hacia afuera
 *   con nextafter despues de restar/sumar M.
 * El predicado angular usa un filtro en double y recurre a __float128 solo si el signo es dudoso. */
#ifndef __IBEX_EXACT_MULT_FAST_H__
#define __IBEX_EXACT_MULT_FAST_H__
#include "ibex_Interval.h"
#include <algorithm>
#include <cmath>
#include <vector>
#include <cstdio>
#include <cstdlib>
namespace ibex { namespace exact_mult_fast {

struct Gen { double u, v; };

inline bool angle_less(const Gen& a, const Gen& b) {
	// filtro: los productos tienen error relativo <= eps salvo subdesborde; un producto nulo es exacto
	// si un factor es 0. Si |d| supera la cota del error, su signo es el exacto.
	const double p = a.u * b.v, q = a.v * b.u, d = p - q;
	const double ap = std::fabs(p), aq = std::fabs(q);
	const bool okp = (a.u == 0.0 || b.v == 0.0) || ap > 1e-280;
	const bool okq = (a.v == 0.0 || b.u == 0.0) || aq > 1e-280;
	if (okp && okq && std::fabs(d) > 4.0 * 2.220446049250313e-16 * (ap + aq)) return d > 0;
	const __float128 c = (__float128) a.u * (__float128) b.v - (__float128) a.v * (__float128) b.u;
	return c > 0;
}

struct Buf { std::vector<Gen> g; };
inline Buf& buf() { static thread_local Buf b; return b; }

/* Recorre el borde del zonotopo (generadores g ordenados) una sola vez y devuelve, para cada desplazamiento
 * (cu[k], cv[k]), k < K, los extremos lo[k], hi[k] de (u+cu)(v+cv) sobre el borde. Como el poligono es
 * cerrado, el extremo q de la ultima arista es el vertice de partida: no hace falta evaluarlo antes. */
template<int K>
inline void walk(const std::vector<Gen>& g, const double* cu, const double* cv, double* lo, double* hi) {
	double su = 0, sv = 0;
	for (const Gen& x : g) { su += x.u; sv += x.v; }
	double pu = -su, pv = -sv;
	for (int k = 0; k < K; k++) { lo[k] = INFINITY; hi[k] = -INFINITY; }
	if (g.empty()) {
		for (int k = 0; k < K; k++) lo[k] = hi[k] = cu[k] * cv[k];
		return;
	}
	for (int sgn = 0; sgn < 2; sgn++) {
		const double f = sgn == 0 ? 2.0 : -2.0;
		for (const Gen& x : g) {
			const double du = f * x.u, dv = f * x.v, A = du * dv;
			for (int k = 0; k < K; k++) {
				const double a = pu + cu[k], b = pv + cv[k];
				const double qa = a + du, qb = b + dv;
				const double vq = qa * qb;
				lo[k] = std::min(lo[k], vq); hi[k] = std::max(hi[k], vq);
				if (A != 0.0) {
					const double B = a * dv + b * du;
					const double t = -B / (2.0 * A);
					if (t > 0.0 && t < 1.0) {
						const double vs = a * b - B * B / (4.0 * A);
						lo[k] = std::min(lo[k], vs); hi[k] = std::max(hi[k], vs);
					}
				}
			}
			pu += du; pv += dv;
		}
	}
}

inline Interval outward(double lo, double hi, double M) {
	return Interval(std::nextafter(lo - M, -INFINITY), std::nextafter(hi + M, INFINITY));
}

inline Interval raf_mult_enclosure_fast(const double* a, const double* b, size_t n,
                                        double x0, double ex, double y0, double ey, bool exact) {
	const double eps = 2.220446049250313e-16;
	std::vector<Gen>& g = buf().g;
	g.clear();
	double Su = 0, Sv = 0;
	for (size_t i = 0; i < n; i++) {
		double u = a[i], v = b[i];
		if (u == 0.0 && v == 0.0) continue;
		if (u < 0.0 || (u == 0.0 && v < 0.0)) { u = -u; v = -v; }
		g.push_back(Gen{u, v}); Su += u; Sv += std::fabs(v);
	}
	const bool coupled = (std::fabs(x0) * ey > 0.0) && (std::fabs(y0) * ex > 0.0);
	if (!exact || !coupled) {
		if (ex > 0) { g.push_back(Gen{ex, 0.0}); Su += ex; }
		if (ey > 0) { g.push_back(Gen{0.0, ey}); Sv += ey; }
		std::sort(g.begin(), g.end(), angle_less);
		const double z = 0.0;
		double lo, hi;
		walk<1>(g, &z, &z, &lo, &hi);
		const double m = (double) g.size();
		const double U = 3.0 * Su * (1 + 8 * eps), V = 3.0 * Sv * (1 + 8 * eps);
		const double lin = (std::fabs(y0) * ex + std::fabs(x0) * ey) * (1 + 4 * eps);
		const double M = (8.0 * m + 64.0) * eps * (U * V + lin) * (1 + 8 * eps);
		return outward(lo - lin, hi + lin, M);
	}
	std::sort(g.begin(), g.end(), angle_less);
	// Esquinas (1,1) y (1,-1): la esquina (-s,-t) da el mismo rango de uv que (s,t) (Z = -Z) con kappa
	// de signo opuesto, asi que cada par aporta [min - |kappa|, max + |kappa|]. Un solo recorrido.
	const double cu[2] = {ex, ex}, cv[2] = {ey, -ey};
	const double k0 = std::fabs(x0 * ey + y0 * ex), k1 = std::fabs(-x0 * ey + y0 * ex);
	double L[2], H[2];
	walk<2>(g, cu, cv, L, H);
	const double lo = std::min(L[0] - k0, L[1] - k1), hi = std::max(H[0] + k0, H[1] + k1);
	const double m = (double) g.size();
	const double U = (ex + 3.0 * Su) * (1 + 8 * eps), V = (ey + 3.0 * Sv) * (1 + 8 * eps);
	const double K = (std::fabs(x0) * ey + std::fabs(y0) * ex) * (1 + 4 * eps);
	const double M = (8.0 * m + 64.0) * eps * (U * V + K) * (1 + 8 * eps);
	return outward(lo, hi, M);
}
}} // namespaces
#endif
