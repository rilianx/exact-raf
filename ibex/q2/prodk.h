// Reduccion avida por bloques producto (PRODK) para acotar z^T H z en [-1,1]^m.
// Port en C++ de notas/raak-qplib/prodk_fast.py (misma busqueda y mismos candidatos), en dobles.
//
//   z^T H z = sum_p (u_p^T z)(v_p^T z) + z^T R z,   rango en  sum_p SH(u_p, v_p) + TT(R)
//
// SH(u, v): rango exacto de (u^T z)(v^T z) (Skalna-Hladik): zonotopo plano con orden angular por producto cruz
// (no atan2) y u, v normalizados por separado. TT(R): diagonal R_ii [0,1], fuera de la diagonal 2 R_ij [-1,1].
//
// Soportes: pares con entrada no nula (o ambas diagonales no nulas), triples con >= 2 aristas (si m <= 40) y el
// soporte completo (si m <= 12). Candidatos por soporte: par espectral (a+b)(a-b), autovector dominante como
// cuadrado y filas de factor comun, cada uno escalado por 1, 3/4, 1/2, 1/4.
//
// Busqueda avida perezosa (lazy = true, por omision): tras extraer un bloque, los soportes que tocan sus
// variables solo se marcan y se recalculan cuando llegan al tope del monticulo. lazy = false recalcula al tiro
// (como prodk_fast.py). Las dos dan cotas validas.
//
// Sin reservas de memoria en el bucle: los soportes tienen a lo mas MAXS = 12 simbolos y se usan arreglos fijos.
//
// Nota de rigor: version para medir precision y costo, sin redondeo dirigido. Para Ibex vale el esquema de
// knomial_quad_bound (parches 0012-0013): la cota es valida para bloques cualesquiera, asi que basta recalcular
// el residuo R y los rangos SH con aritmetica de intervalos al final.
#pragma once
#include <algorithm>
#include <cmath>
#include <queue>
#include <vector>

namespace prodk {

typedef std::vector<double> Vec;
static const int MAXS = 12;

struct Pair { int s; double u[MAXS], v[MAXS]; };   // coordenadas locales del soporte

// ---------------------------------------------------------------- SH
inline void sh_range(const double* u0, const double* v0, int n, double& lo, double& hi) {
	double su = 0, sv = 0;
	for (int i = 0; i < n; i++) { su = std::max(su, std::fabs(u0[i])); sv = std::max(sv, std::fabs(v0[i])); }
	lo = hi = 0;
	if (su == 0 || sv == 0) return;
	double ga[MAXS + 2], gb[MAXS + 2]; int ord[MAXS + 2];
	int k = 0; double sa = 0, sb = 0;
	for (int i = 0; i < n; i++) {
		double a = u0[i] / su, b = v0[i] / sv;
		if (a == 0 && b == 0) continue;
		if (!(b > 0 || (b == 0 && a > 0))) { a = -a; b = -b; }
		ga[k] = a; gb[k] = b; ord[k] = k; k++; sa += a; sb += b;
	}
	std::sort(ord, ord + k, [&](int p, int q) { return ga[p] * gb[q] - gb[p] * ga[q] > 0; });
	double X[2 * (MAXS + 2)], Y[2 * (MAXS + 2)];
	double px = -sa, py = -sb;
	for (int i = 0; i < k; i++) { X[i] = px; Y[i] = py; px += 2 * ga[ord[i]]; py += 2 * gb[ord[i]]; }
	for (int i = 0; i < k; i++) { X[k + i] = px; Y[k + i] = py; px -= 2 * ga[ord[i]]; py -= 2 * gb[ord[i]]; }
	lo = INFINITY; hi = -INFINITY;
	const int nv = 2 * k;
	for (int i = 0; i < nv; i++) {
		const int j = (i + 1 == nv) ? 0 : i + 1;
		const double x0 = X[i], y0 = Y[i], dx = X[j] - x0, dy = Y[j] - y0;
		double f = x0 * y0; lo = std::min(lo, f); hi = std::max(hi, f);
		const double den = 2 * dx * dy;
		if (den != 0) {
			const double t = -(x0 * dy + y0 * dx) / den;
			if (t > 0 && t < 1) { f = (x0 + t * dx) * (y0 + t * dy); lo = std::min(lo, f); hi = std::max(hi, f); }
		}
	}
	lo *= su * sv; hi *= su * sv;
}
inline void sh_range(const Pair& p, double& lo, double& hi) { sh_range(p.u, p.v, p.s, lo, hi); }
inline double sh_width(const Pair& p) { double lo, hi; sh_range(p, lo, hi); return hi - lo; }

// ---------------------------------------------------------------- utilidades
inline double tt(const double* R, int s) {
	double d = 0, a = 0;
	for (int i = 0; i < s; i++) for (int j = 0; j < s; j++) {
		const double v = std::fabs(R[i * s + j]);
		if (i == j) d += v; else a += v;
	}
	return d + 2 * a;
}

// Jacobi ciclico para simetricas s x s (s <= MAXS). V por columnas: V[i*s + k].
inline void jacobi(const double* A0, int s, double* lam, double* V) {
	double A[MAXS * MAXS];
	std::copy(A0, A0 + s * s, A);
	std::fill(V, V + s * s, 0.0);
	for (int i = 0; i < s; i++) V[i * s + i] = 1.0;
	double nrm = 0;
	for (int i = 0; i < s * s; i++) nrm += A[i] * A[i];
	for (int sweep = 0; sweep < 60; sweep++) {
		double off = 0;
		for (int p = 0; p < s; p++) for (int q = p + 1; q < s; q++) off += A[p * s + q] * A[p * s + q];
		if (off <= 1e-30 * nrm) break;
		for (int p = 0; p < s; p++) for (int q = p + 1; q < s; q++) {
			const double apq = A[p * s + q];
			if (std::fabs(apq) <= 1e-300) continue;
			const double theta = (A[q * s + q] - A[p * s + p]) / (2 * apq);
			const double t = (theta >= 0 ? 1.0 : -1.0) / (std::fabs(theta) + std::sqrt(theta * theta + 1));
			const double c = 1 / std::sqrt(t * t + 1), sn = t * c;
			for (int k = 0; k < s; k++) {
				const double akp = A[k * s + p], akq = A[k * s + q];
				A[k * s + p] = c * akp - sn * akq; A[k * s + q] = sn * akp + c * akq;
			}
			for (int k = 0; k < s; k++) {
				const double apk = A[p * s + k], aqk = A[q * s + k];
				A[p * s + k] = c * apk - sn * aqk; A[q * s + k] = sn * apk + c * aqk;
			}
			for (int k = 0; k < s; k++) {
				const double vkp = V[k * s + p], vkq = V[k * s + q];
				V[k * s + p] = c * vkp - sn * vkq; V[k * s + q] = sn * vkp + c * vkq;
			}
		}
	}
	for (int i = 0; i < s; i++) lam[i] = A[i * s + i];
}

// candidatos sobre RS (s x s): devuelve cuantos escribio en out (cada candidato es un solo par)
inline int candidates(const double* RS, int s, Pair* out) {
	double lam[MAXS], V[MAXS * MAXS];
	jacobi(RS, s, lam, V);
	double mx = 0;
	for (int k = 0; k < s; k++) mx = std::max(mx, std::fabs(lam[k]));
	const double sc = std::max(1e-300, mx) * 1e-12;
	int ip = 0, in = 0, ia = 0;
	for (int k = 1; k < s; k++) {
		if (lam[k] > lam[ip]) ip = k;
		if (lam[k] < lam[in]) in = k;
		if (std::fabs(lam[k]) > std::fabs(lam[ia])) ia = k;
	}
	int n = 0;
	if (lam[ip] > sc && lam[in] < -sc) {
		const double fa = std::sqrt(lam[ip]), fb = std::sqrt(-lam[in]);
		Pair& p = out[n++]; p.s = s;
		for (int i = 0; i < s; i++) {
			const double a = fa * V[i * s + ip], b = fb * V[i * s + in];
			p.u[i] = a + b; p.v[i] = a - b;
		}
	}
	if (std::fabs(lam[ia]) > sc) {
		const double fa = std::sqrt(std::fabs(lam[ia])), sg = lam[ia] < 0 ? -1.0 : 1.0;
		Pair& p = out[n++]; p.s = s;
		for (int i = 0; i < s; i++) { p.u[i] = fa * V[i * s + ia]; p.v[i] = sg * p.u[i]; }
	}
	for (int t = 0; t < s; t++) {
		Pair& p = out[n]; p.s = s;
		double sa = 0;
		for (int j = 0; j < s; j++) { p.v[j] = 2 * RS[t * s + j]; p.u[j] = 0.0; }
		p.v[t] = RS[t * s + t]; p.u[t] = 1.0;
		for (int j = 0; j < s; j++) sa += std::fabs(p.v[j]);
		if (sa > 0) n++;
	}
	return n;
}

// mejor candidato en el soporte S (tamano s) de R (m x m): ganancia y par (en best.s = 0 si no hay)
inline double best_on(const Vec& R, int m, const int* S, int s, Pair& best) {
	double RS[MAXS * MAXS], W[MAXS * MAXS];
	for (int i = 0; i < s; i++) for (int j = 0; j < s; j++) RS[i * s + j] = R[S[i] * m + S[j]];
	best.s = 0;
	const double base = tt(RS, s);
	if (base == 0) return 0.0;
	static const double TS[4] = {1.0, 0.75, 0.5, 0.25};
	Pair cand[MAXS + 2];
	const int nc = candidates(RS, s, cand);
	double gain = 0.0;
	for (int c = 0; c < nc; c++) {
		const Pair& p0 = cand[c];
		const double w0 = sh_width(p0);            // SH(t u, v) = t SH(u, v) para t > 0: un solo calculo
		double M[MAXS * MAXS];
		for (int i = 0; i < s; i++) for (int j = 0; j < s; j++) M[i * s + j] = 0.5 * (p0.u[i] * p0.v[j] + p0.v[i] * p0.u[j]);
		for (double t : TS) {
			for (int i = 0; i < s * s; i++) W[i] = RS[i] - t * M[i];
			const double g = base - tt(W, s) - t * w0;
			if (g > gain) {
				gain = g; best = p0;
				for (int i = 0; i < s; i++) best.u[i] *= t;
			}
		}
	}
	return gain;
}

// ---------------------------------------------------------------- reduccion
struct Result { double L, U; int blocks; };

inline Result prodk_interval(const Vec& H, int m, bool lazy = true, int max_triples_m = 40, int full_max = MAXS) {
	Vec R = H;
	auto nz = [&](int i, int j) { return H[i * m + j] != 0.0; };
	// soportes en un arreglo plano: inicio y tamano de cada uno
	std::vector<int> sv, sstart, ssize;
	auto add = [&](std::initializer_list<int> l) { sstart.push_back((int) sv.size()); ssize.push_back((int) l.size()); sv.insert(sv.end(), l); };
	for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++)
		if (nz(i, j) || (nz(i, i) && nz(j, j))) add({i, j});
	if (m <= max_triples_m)
		for (int i = 0; i < m; i++) for (int j = i + 1; j < m; j++) for (int k = j + 1; k < m; k++)
			if ((int) nz(i, j) + (int) nz(i, k) + (int) nz(j, k) >= 2) add({i, j, k});
	if (m <= full_max && m > 2) {
		sstart.push_back((int) sv.size()); ssize.push_back(m);
		for (int i = 0; i < m; i++) sv.push_back(i);
	}
	const int ns = (int) sstart.size();
	std::vector<std::vector<int>> by_var(m);
	for (int s = 0; s < ns; s++) for (int t = 0; t < ssize[s]; t++) by_var[sv[sstart[s] + t]].push_back(s);
	std::vector<int> ver(ns, 0), done(ns, -1), mark(ns, -1);
	std::vector<Pair> cache(ns);
	typedef std::pair<double, std::pair<int, int>> HI;              // (ganancia, (-soporte, version))
	std::priority_queue<HI> heap;
	auto push = [&](int s) {
		const double g = best_on(R, m, &sv[sstart[s]], ssize[s], cache[s]);
		done[s] = ver[s];
		if (cache[s].s) heap.push(HI(g, std::make_pair(-s, ver[s])));
	};
	for (int s = 0; s < ns; s++) push(s);
	const double scale = std::max(1.0, tt(H.data(), m)) * 1e-12;
	double Lsh = 0, Ush = 0;
	int nblocks = 0;
	for (int it = 0; it < 1000000; it++) {
		if (heap.empty()) break;
		const HI top = heap.top();
		const int s = -top.second.first;
		if (top.second.second != ver[s]) {
			heap.pop();
			if (lazy && done[s] != ver[s]) push(s);
			continue;
		}
		if (top.first <= scale) break;
		heap.pop();
		const int* S = &sv[sstart[s]];
		const int ss = ssize[s];
		const Pair p = cache[s];
		for (int i = 0; i < ss; i++) for (int j = 0; j < ss; j++)
			R[S[i] * m + S[j]] -= 0.5 * (p.u[i] * p.v[j] + p.v[i] * p.u[j]);
		double lo, hi; sh_range(p, lo, hi); Lsh += lo; Ush += hi; nblocks++;
		for (int i = 0; i < ss; i++) for (int t : by_var[S[i]]) if (mark[t] != it) {
			mark[t] = it; ver[t]++;
			if (!lazy) push(t);
		}
		if (lazy) push(s);
	}
	double dlo = 0, dhi = 0, off = 0;
	for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {
		const double v = R[i * m + j];
		if (i == j) { dlo += std::min(v, 0.0); dhi += std::max(v, 0.0); } else off += std::fabs(v);
	}
	return {Lsh + dlo - off, Ush + dhi + off, nblocks};
}

} // namespace prodk
