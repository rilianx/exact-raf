// Prototipo en C++ (fuera de Ibex) para comparar, con precision y tiempo, en las mismas cintas y cajas:
//   M00  fAF2 de Ibex por omision: producto con la formula (5) de Skalna-Hladik (regla 'af2')
//   M0   fAF2 con producto SH (rango exacto de la parte cuadratica, errores como generadores propios)
//   M3c  formas afines de grado 2 (c + a^T z + z^T Q z, z = variables y simbolos de error), resto de grado 3-4
//        centrado, con cuatro encierros finales: tt (barato), kn6 (k-nomial k=6 de ibex_KNomial.h),
//        pc (bloques producto: pares + soporte completo, m <= 12) y full (pares + triples + completo).
// Funciones (iguales en todos los modos): min-range para exp, log, sqrt, 1/x; secante con residuo analitico
// (puntos criticos) para sin, cos, abs y x^k. Medicion en dobles, sin redondeo dirigido.
//
// Entrada: cintas de emitir_cintas.py. Salida (CSV por stdout): una fila por (expresion, caja).
// Uso: aa2 <cintas.txt> [REPS=3]
#include <algorithm>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include "ibex_KNomial.h"
#include "prodk.h"

using std::vector;
typedef vector<double> Vec;

// ------------------------------------------------------------------ linealizacion de funciones
struct Lin { double alpha, d1, d2; bool ok; };

static Lin lin_crit(double (*f)(double, double), double prm, double lo, double hi, double alpha, const vector<double>& crit) {
	double d1 = INFINITY, d2 = -INFINITY;
	auto upd = [&](double t) { const double d = f(t, prm) - alpha * t; d1 = std::min(d1, d); d2 = std::max(d2, d); };
	upd(lo); upd(hi);
	for (double t : crit) if (t > lo && t < hi) upd(t);
	return {alpha, d1, d2, true};
}
static double f_sin(double t, double) { return std::sin(t); }
static double f_cos(double t, double) { return std::cos(t); }
static double f_abs(double t, double) { return std::fabs(t); }
static double f_pow(double t, double k) { return std::pow(t, k); }

static Lin linearize(const std::string& op, double lo, double hi, double k = 0) {
	if (op == "EXP") { const double a = std::exp(lo); return {a, std::exp(lo) - a * lo, std::exp(hi) - a * hi, true}; }
	if (op == "LOG" || op == "SQRT") {
		if (hi <= 0) return {0, 0, 0, false};
		lo = std::max(lo, 1e-300);
		if (op == "LOG") { const double a = 1 / hi; return {a, std::log(lo) - a * lo, std::log(hi) - a * hi, true}; }
		const double a = 0.5 / std::sqrt(hi); return {a, std::sqrt(lo) - a * lo, std::sqrt(hi) - a * hi, true};
	}
	if (op == "REC") {
		if (lo <= 0 && hi >= 0) return {0, 0, 0, false};
		const double a = -1 / (lo * hi), t0 = std::copysign(std::sqrt(lo * hi), lo);
		double d1 = INFINITY, d2 = -INFINITY;
		for (double t : {lo, hi, t0}) if (t >= lo && t <= hi) { const double d = 1 / t - a * t; d1 = std::min(d1, d); d2 = std::max(d2, d); }
		return {a, d1, d2, true};
	}
	if (op == "ABS") {
		if (lo >= 0) return {1, 0, 0, true};
		if (hi <= 0) return {-1, 0, 0, true};
		const double a = (std::fabs(hi) - std::fabs(lo)) / (hi - lo);
		return lin_crit(f_abs, 0, lo, hi, a, {0.0});
	}
	if (op == "SIN" || op == "COS") {
		const bool s = op == "SIN";
		double (*f)(double, double) = s ? f_sin : f_cos;
		const double a = (f(hi, 0) - f(lo, 0)) / (hi - lo);
		vector<double> crit;
		if (std::fabs(a) <= 1) {
			// sin: cos t = a -> t = +-acos(a) + 2 pi n ;  cos: -sin t = a -> t = asin(-a) + 2 pi n, pi - asin(-a) + 2 pi n
			const double b1 = s ? std::acos(a) : std::asin(-a), b2 = s ? -std::acos(a) : M_PI - std::asin(-a);
			const long n0 = (long) std::floor((lo - 4 * M_PI) / (2 * M_PI)), n1 = (long) std::ceil((hi + 4 * M_PI) / (2 * M_PI));
			if (n1 - n0 > 400) {      // intervalo enorme: cota gruesa valida
				const double m1 = std::min(a * lo, a * hi), m2 = std::max(a * lo, a * hi);
				return {a, -1 - m2, 1 - m1, true};
			}
			for (long n = n0; n <= n1; n++) { crit.push_back(b1 + 2 * M_PI * n); crit.push_back(b2 + 2 * M_PI * n); }
		}
		return lin_crit(f, 0, lo, hi, a, crit);
	}
	if (op == "POW") {
		const double a = (std::pow(hi, k) - std::pow(lo, k)) / (hi - lo);
		vector<double> crit;            // k t^(k-1) = a
		const double r = a / k; const int e = (int) k - 1;
		if (e % 2 == 1) { const double t = std::copysign(std::pow(std::fabs(r), 1.0 / e), r); crit.push_back(t); }
		else if (r >= 0) { const double t = std::pow(r, 1.0 / e); crit.push_back(t); crit.push_back(-t); }
		return lin_crit(f_pow, k, lo, hi, a, crit);
	}
	return {0, 0, 0, false};
}

static double fval(const std::string& op, double x, double k = 0) {
	if (op == "EXP") return std::exp(x); if (op == "LOG") return std::log(x); if (op == "SQRT") return std::sqrt(x);
	if (op == "REC") return 1 / x; if (op == "SIN") return std::sin(x); if (op == "COS") return std::cos(x);
	if (op == "ABS") return std::fabs(x); return std::pow(x, k);
}

// ------------------------------------------------------------------ cintas
struct Node { std::string op; int a = -1, b = -1; double c = 0; };
struct Expr { std::string name; vector<Node> nodes; };
struct Inst { std::string name; int n; vector<int> walk, depth; vector<Vec> lo, hi; vector<Expr> exprs; };

static vector<Inst> read_tapes(const char* path) {
	std::ifstream in(path);
	vector<Inst> out;
	std::string line;
	while (std::getline(in, line)) {
		std::istringstream is(line);
		std::string tag; is >> tag;
		if (tag == "I") { Inst I; int nb; is >> I.name >> I.n >> nb; out.push_back(I); }
		else if (tag == "B") {
			Inst& I = out.back(); int w, d; is >> w >> d;
			Vec lo(I.n), hi(I.n);
			for (int i = 0; i < I.n; i++) is >> lo[i] >> hi[i];
			I.walk.push_back(w); I.depth.push_back(d); I.lo.push_back(lo); I.hi.push_back(hi);
		} else if (tag == "E") {
			Expr E; int nn; is >> E.name >> nn;
			for (int k = 0; k < nn; k++) {
				std::getline(in, line); std::istringstream ns(line);
				Node N; ns >> N.op;
				if (N.op == "V") ns >> N.a;
				else if (N.op == "C") ns >> N.c;
				else if (N.op == "ADD" || N.op == "SUB" || N.op == "MUL") ns >> N.a >> N.b;
				else if (N.op == "SCL") ns >> N.a >> N.c;
				else if (N.op == "POW") ns >> N.a >> N.c;
				else ns >> N.a;
				E.nodes.push_back(N);
			}
			out.back().exprs.push_back(E);
		}
	}
	return out;
}

// ------------------------------------------------------------------ fAF2 (M00, M0)
struct AF { double c; Vec a; double e; };

static AF af_mul(const AF& x, const AF& y, bool sh) {
	const size_t n = x.a.size();
	AF r{x.c * y.c, Vec(n), 0};
	for (size_t i = 0; i < n; i++) r.a[i] = y.c * x.a[i] + x.c * y.a[i];
	if (!sh) {   // formula (5): centro + 1/2 sum a_i b_i, error (Sx+ex)(Sy+ey) - 1/2 sum |a_i b_i| + |x0| ey + |y0| ex
		double Sx = 0, Sy = 0, sab = 0, sabs = 0;
		for (size_t i = 0; i < n; i++) { Sx += std::fabs(x.a[i]); Sy += std::fabs(y.a[i]); sab += x.a[i] * y.a[i]; sabs += std::fabs(x.a[i] * y.a[i]); }
		r.c += 0.5 * sab;
		r.e = std::fabs(x.c) * y.e + std::fabs(y.c) * x.e + (Sx + x.e) * (Sy + y.e) - 0.5 * sabs;
	} else {     // SH: rango exacto de (a^T eps + ex xi)(b^T eps + ey eta)
		prodk::Pair p; double lo, hi;
		vector<double> u(n + 2, 0.0), v(n + 2, 0.0);
		for (size_t i = 0; i < n; i++) { u[i] = x.a[i]; v[i] = y.a[i]; }
		u[n] = x.e; v[n + 1] = y.e;
		prodk::sh_range(u.data(), v.data(), (int) (n + 2), lo, hi);
		r.c += 0.5 * (lo + hi);
		r.e = std::fabs(y.c) * x.e + std::fabs(x.c) * y.e + 0.5 * (hi - lo);
	}
	return r;
}

static bool eval_af(const Expr& E, const Vec& lo, const Vec& hi, bool sh, double& wlo, double& whi) {
	const size_t n = lo.size();
	vector<AF> v(E.nodes.size());
	for (size_t k = 0; k < E.nodes.size(); k++) {
		const Node& N = E.nodes[k];
		AF r{0, Vec(n, 0.0), 0};
		if (N.op == "V") { r.c = 0.5 * (lo[N.a] + hi[N.a]); r.a[N.a] = 0.5 * (hi[N.a] - lo[N.a]); }
		else if (N.op == "C") r.c = N.c;
		else if (N.op == "ADD" || N.op == "SUB") {
			const double s = N.op == "ADD" ? 1 : -1; const AF &x = v[N.a], &y = v[N.b];
			r.c = x.c + s * y.c; r.e = x.e + y.e;
			for (size_t i = 0; i < n; i++) r.a[i] = x.a[i] + s * y.a[i];
		} else if (N.op == "SCL" || N.op == "NEG") {
			const double s = N.op == "NEG" ? -1 : N.c; const AF& x = v[N.a];
			r.c = s * x.c; r.e = std::fabs(s) * x.e;
			for (size_t i = 0; i < n; i++) r.a[i] = s * x.a[i];
		} else if (N.op == "MUL") r = af_mul(v[N.a], v[N.b], sh);
		else {
			const AF& x = v[N.a];
			double rad = x.e; for (double q : x.a) rad += std::fabs(q);
			const double l = x.c - rad, h = x.c + rad;
			if (h - l <= 1e-15 * std::max(1.0, std::fabs(l))) {
				if ((N.op == "REC" && x.c == 0) || ((N.op == "LOG" || N.op == "SQRT") && x.c <= 0)) return false;
				r.c = fval(N.op, x.c, N.c);
			} else {
				const Lin L = linearize(N.op, l, h, N.c);
				if (!L.ok) return false;
				r.c = L.alpha * x.c + 0.5 * (L.d1 + L.d2); r.e = std::fabs(L.alpha) * x.e + 0.5 * (L.d2 - L.d1);
				for (size_t i = 0; i < n; i++) r.a[i] = L.alpha * x.a[i];
			}
		}
		if (!std::isfinite(r.c) || !std::isfinite(r.e)) return false;
		v[k] = std::move(r);
	}
	const AF& f = v.back();
	double rad = f.e; for (double q : f.a) rad += std::fabs(q);
	wlo = f.c - rad; whi = f.c + rad;
	return std::isfinite(wlo) && std::isfinite(whi);
}

// ------------------------------------------------------------------ formas de grado 2 (M3c)
typedef vector<std::pair<int, double>> SV;                 // vector disperso ordenado por indice
typedef vector<std::pair<long long, double>> SQ;           // Q dispersa: clave i*2^32 + j (i <= j), coef. del monomio
static inline long long qk(int i, int j) { if (i > j) std::swap(i, j); return ((long long) i << 32) | (unsigned) j; }

template <class T, class K> static T merge(const T& x, double sx, const T& y, double sy) {
	T r; r.reserve(x.size() + y.size());
	size_t i = 0, j = 0;
	while (i < x.size() || j < y.size()) {
		if (j == y.size() || (i < x.size() && x[i].first < y[j].first)) { r.push_back({x[i].first, sx * x[i].second}); i++; }
		else if (i == x.size() || y[j].first < x[i].first) { r.push_back({y[j].first, sy * y[j].second}); j++; }
		else { r.push_back({x[i].first, sx * x[i].second + sy * y[j].second}); i++; j++; }
	}
	return r;
}

struct QF { double c = 0; SV lin; SQ Q; };

static void quad_interval(const QF& X, double& lo, double& hi) {
	lo = hi = 0;
	for (auto& t : X.Q) {
		const int i = (int) (t.first >> 32), j = (int) (t.first & 0xffffffff);
		if (i == j) { lo += std::min(0.0, t.second); hi += std::max(0.0, t.second); }
		else { lo -= std::fabs(t.second); hi += std::fabs(t.second); }
	}
}
static double lin_abs(const QF& X) { double s = 0; for (auto& t : X.lin) s += std::fabs(t.second); return s; }

static void qf_range(const QF& X, double& lo, double& hi) {
	lo = hi = X.c;
	std::map<int, double> diag;
	double off = 0;
	for (auto& t : X.Q) {
		const int i = (int) (t.first >> 32), j = (int) (t.first & 0xffffffff);
		if (i == j) diag[i] += t.second; else off += std::fabs(t.second);
	}
	size_t k = 0;
	auto one = [&](double a, double q) {
		double vmin = std::min(q + a, q - a), vmax = std::max(q + a, q - a);
		if (q != 0 && std::fabs(a / (2 * q)) < 1) { const double s = -a * a / (4 * q); vmin = std::min(vmin, s); vmax = std::max(vmax, s); }
		lo += vmin; hi += vmax;
	};
	for (auto& t : X.lin) {
		auto it = diag.find(t.first);
		if (it != diag.end()) { one(t.second, it->second); diag.erase(it); } else one(t.second, 0.0);
		k++;
	}
	for (auto& d : diag) one(0.0, d.second);
	lo -= off; hi += off;
}

static QF qf_scale(const QF& X, double s) {
	QF r; r.c = s * X.c; r.lin = X.lin; r.Q = X.Q;
	for (auto& t : r.lin) t.second *= s;
	for (auto& t : r.Q) t.second *= s;
	return r;
}

struct EvalQ {
	int nsym;
	explicit EvalQ(int n) : nsym(n) {}
	void add_err(QF& F, double e) { if (e > 0) F.lin.push_back({nsym++, e}); }   // indices crecientes: sigue ordenado

	QF mul(const QF& X, const QF& Y, bool same) {
		double qx0, qx1, qy0, qy1;
		quad_interval(X, qx0, qx1); quad_interval(Y, qy0, qy1);
		const double mx = 0.5 * (qx0 + qx1), my = 0.5 * (qy0 + qy1), rx = 0.5 * (qx1 - qx0), ry = 0.5 * (qy1 - qy0);
		const double Cx = X.c + mx, Cy = Y.c + my;
		QF r;
		r.lin = merge<SV, int>(X.lin, Cy, Y.lin, Cx);
		SQ q = merge<SQ, long long>(X.Q, Cy, Y.Q, Cx);
		SQ outer; outer.reserve(X.lin.size() * Y.lin.size());
		for (auto& a : X.lin) for (auto& b : Y.lin) outer.push_back({qk(a.first, b.first), a.second * b.second});
		std::sort(outer.begin(), outer.end(), [](const std::pair<long long, double>& p, const std::pair<long long, double>& s) { return p.first < s.first; });
		SQ oc; for (auto& t : outer) { if (!oc.empty() && oc.back().first == t.first) oc.back().second += t.second; else oc.push_back(t); }
		r.Q = merge<SQ, long long>(q, 1.0, oc, 1.0);
		const double sa = lin_abs(X), sb = lin_abs(Y);
		double lo, hi;
		if (same) { lo = -2 * sa * rx; hi = 2 * sa * rx + rx * rx; }
		else { const double t = sa * ry + sb * rx + rx * ry; lo = -t; hi = t; }
		r.c = Cx * Cy - Cx * my - Cy * mx + 0.5 * (lo + hi);
		add_err(r, 0.5 * (hi - lo));
		return r;
	}

	bool eval(const Expr& E, const Vec& blo, const Vec& bhi, QF& out) {
		vector<QF> v(E.nodes.size());
		for (size_t k = 0; k < E.nodes.size(); k++) {
			const Node& N = E.nodes[k];
			QF r;
			if (N.op == "V") { r.c = 0.5 * (blo[N.a] + bhi[N.a]); const double rd = 0.5 * (bhi[N.a] - blo[N.a]); if (rd > 0) r.lin.push_back({N.a, rd}); }
			else if (N.op == "C") r.c = N.c;
			else if (N.op == "ADD" || N.op == "SUB") {
				const double s = N.op == "ADD" ? 1 : -1; const QF &x = v[N.a], &y = v[N.b];
				r.c = x.c + s * y.c; r.lin = merge<SV, int>(x.lin, 1.0, y.lin, s); r.Q = merge<SQ, long long>(x.Q, 1.0, y.Q, s);
			} else if (N.op == "SCL" || N.op == "NEG") r = qf_scale(v[N.a], N.op == "NEG" ? -1.0 : N.c);
			else if (N.op == "MUL") r = mul(v[N.a], v[N.b], N.a == N.b);
			else {
				const QF& x = v[N.a];
				double l, h; qf_range(x, l, h);
				if (h - l <= 1e-15 * std::max(1.0, std::fabs(l))) {
					if ((N.op == "REC" && x.c == 0) || ((N.op == "LOG" || N.op == "SQRT") && x.c <= 0)) return false;
					r.c = fval(N.op, 0.5 * (l + h), N.c);
				} else {
					const Lin L = linearize(N.op, l, h, N.c);
					if (!L.ok) return false;
					r = qf_scale(x, L.alpha); r.c += 0.5 * (L.d1 + L.d2);
					add_err(r, 0.5 * (L.d2 - L.d1));
				}
			}
			if (!std::isfinite(r.c)) return false;
			v[k] = std::move(r);
		}
		out = v.back();
		return true;
	}
};

// homogeneizada H (m+1)x(m+1) por filas: z^T Q z + a^T z = [z;t]^T H [z;t] con t = 1 (mismo rango con t en [-1,1])
static int homog(const QF& F, Vec& H) {
	std::map<int, int> idx;
	for (auto& t : F.lin) idx[t.first] = 0;
	for (auto& t : F.Q) { idx[(int) (t.first >> 32)] = 0; idx[(int) (t.first & 0xffffffff)] = 0; }
	int m = 0; for (auto& p : idx) p.second = m++;
	const int s = m + 1;
	H.assign(s * s, 0.0);
	for (auto& t : F.Q) {
		const int a = idx[(int) (t.first >> 32)], b = idx[(int) (t.first & 0xffffffff)];
		if (a == b) H[a * s + a] += t.second; else { H[a * s + b] += 0.5 * t.second; H[b * s + a] += 0.5 * t.second; }
	}
	for (auto& t : F.lin) { const int a = idx[t.first]; H[a * s + m] += 0.5 * t.second; H[m * s + a] += 0.5 * t.second; }
	return s;
}

// k-nomial sobre H (como bound_full de ../raak-medicion-synq/knomial_driver.cpp, con redondeo hacia +inf)
static double knomial_width(const Vec& H0, int n, int k) {
	Vec Q = H0, w(n, 2.0);
	vector<ibex::KNomialBlock> blocks;
	ibex::knomial_rewrite(Q, w, n, k, blocks);
	const int old = std::fegetround(); std::fesetround(FE_UPWARD);
	double u = 0, negl = 0;
	for (auto& b : blocks) { double s = 0; for (int i = 0; i < n; i++) s += std::fabs(b.u[i]); if (b.c > 0) u += s * s; else negl += s * s; }
	for (int i = 0; i < n; i++) for (int j = i; j < n; j++) {
		const double q = 0.5 * (H0[i * n + j] + H0[j * n + i]);
		double tu = 0, nt = 0;
		for (auto& b : blocks) if (b.u[i] != 0 && b.u[j] != 0) { const double si = b.c * b.u[i]; tu += si * b.u[j]; nt += (-si) * b.u[j]; }
		const double rup = q + nt, rlo = -(tu - q);
		if (i == j) { if (rup > 0) u += rup; if (rlo < 0) negl += -rlo; }
		else { const double m2 = 2.0 * std::max(std::fabs(rlo), std::fabs(rup)); u += m2; negl += m2; }
	}
	std::fesetround(old);
	return u + negl;
}

// ------------------------------------------------------------------ medicion
template <class F> static double tmed(int reps, F f) {
	vector<double> t;
	for (int r = 0; r < reps; r++) {
		const auto t0 = std::chrono::steady_clock::now(); f();
		t.push_back(std::chrono::duration<double, std::micro>(std::chrono::steady_clock::now() - t0).count());
	}
	std::sort(t.begin(), t.end()); return t[t.size() / 2];
}

int main(int argc, char** argv) {
	ibex::knomial_bound_mode() = 1;
	const int REPS = argc > 2 ? std::atoi(argv[2]) : 3;
	const vector<Inst> I = read_tapes(argv[1]);
	std::printf("instancia,expr,caminata,prof,m,M00,M0,M3c_tt,M3c_kn6,M3c_pc,M3c_full,t_M00,t_M0,t_M3c_eval,t_tt,t_kn6,t_pc,t_full\n");
	for (const Inst& In : I)
		for (size_t b = 0; b < In.lo.size(); b++)
			for (const Expr& E : In.exprs) {
				double l00, h00, l0, h0;
				if (!eval_af(E, In.lo[b], In.hi[b], false, l00, h00) || !eval_af(E, In.lo[b], In.hi[b], true, l0, h0)) continue;
				QF F; EvalQ ev(In.n);
				if (!ev.eval(E, In.lo[b], In.hi[b], F)) continue;
				const double t00 = tmed(REPS, [&] { double a, c; eval_af(E, In.lo[b], In.hi[b], false, a, c); });
				const double t0 = tmed(REPS, [&] { double a, c; eval_af(E, In.lo[b], In.hi[b], true, a, c); });
				const double t3 = tmed(REPS, [&] { QF G; EvalQ e2(In.n); e2.eval(E, In.lo[b], In.hi[b], G); });
				double lt, ht; qf_range(F, lt, ht);
				const double wtt = ht - lt, ttt = tmed(REPS, [&] { double a, c; qf_range(F, a, c); });
				double wkn = wtt, wpc = wtt, wfu = wtt, tkn = 0, tpc = 0, tfu = 0;
				int m = 0;
				if (!F.Q.empty()) {
					Vec H; m = homog(F, H);
					tkn = tmed(REPS, [&] { Vec H2; const int s = homog(F, H2); wkn = std::min(wtt, knomial_width(H2, s, 6)); });
					tpc = tmed(REPS, [&] { Vec H2; const int s = homog(F, H2); auto r = prodk::prodk_interval(H2, s, true, 0, 12); wpc = std::min(wtt, r.U - r.L); });
					tfu = tmed(REPS, [&] { Vec H2; const int s = homog(F, H2); auto r = prodk::prodk_interval(H2, s, true, 40, 12); wfu = std::min(wtt, r.U - r.L); });
				}
				std::printf("%s,%s,%d,%d,%d,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n",
				            In.name.c_str(), E.name.c_str(), In.walk[b], In.depth[b], m, h00 - l00, h0 - l0, wtt, wkn, wpc, wfu,
				            t00, t0, t3, ttt, tkn, tpc, tfu);
			}
}
