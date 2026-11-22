//============================================================================
// Formas afines de grado 2 para Ibex (prototipo de integracion, rama raak)
//
//   F = c + a^T z + z^T Q z,   z = (eps de las variables, simbolos de error nuevos) en [-1,1]
//
// Producto: se conserva todo el grado 2; cada parte cuadratica q se separa como q = m + qhat (m = centro de su
// rango) y el resto de grado 3-4 se acota con qhat y pasa a un simbolo de error nuevo. Funciones: Chebyshev
// (pendiente de la secante y residuo analitico en los puntos criticos) para exp, log, sqrt, 1/x, sin, cos, abs,
// x^k, como el modo por omision de ibex-affine; las demas operaciones de Ibex se evaluan por intervalos.
// Como en AffineEval, en cada nodo se lleva un intervalo (operacion por intervalos & rango de la forma) que es el
// dominio de las linealizaciones y la evaluacion 'ev' del linearizador.
//
// Para la relajacion lineal: la parte lineal en las variables se conserva y el resto
// r(z) = sum_{simbolos de error} a_e e + z^T Q z se acota con un encierro final:
//   TT    barato (cada simbolo exacto con su lineal y su diagonal + resto termino a termino)
//   KN    k-nomial (k=3 por omision, variable Q2_KN_K; ibex_KNomial.h, parches 0012-0015)
//   PC    bloques producto, pares + soporte completo (m <= 12)   (prodk.h)
//   FULL  bloques producto, pares + triples + completo            (prodk.h)
// y se intersecta con TT.
//
// RIGOR: la aritmetica es en dobles. Para la relajacion se agrega un margen relativo QUAD2_MARGIN (1e-13) a la
// cota del resto y a la constante; es un margen de seguridad, no una cota rigurosa del redondeo (pendiente: el
// esquema de knomial_quad_bound, residuo en intervalos).
//============================================================================
#ifndef __IBEX_QUAD2_H__
#define __IBEX_QUAD2_H__

#include "ibex_Function.h"
#include "ibex_Expr.h"
#include "ibex_KNomial.h"
#include "prodk.h"
#include <algorithm>
#include <cstdlib>
#include <cmath>
#include <map>
#include <string>
#include <vector>

namespace ibex {
namespace quad2 {

enum FinalBound { FB_TT = 0, FB_KN = 1, FB_PC = 2, FB_FULL = 3 };
// Margenes de seguridad (no rigurosos) frente al redondeo en dobles: relativos a la magnitud acumulada del calculo
// (suma del tamano de los terminos de todos los nodos), que ya crece con el numero de operaciones. El redondeo es
// del orden de 1e-16 por operacion; se usa 1e-13 (mil veces mas).
static const double QUAD2_MARGIN = 1e-13;

typedef std::vector<double> Vec;
typedef std::vector<std::pair<int, double> > SV;
typedef std::vector<std::pair<long long, double> > SQ;
inline long long qk(int i, int j) { if (i > j) std::swap(i, j); return ((long long) i << 32) | (unsigned) j; }
inline int qi(long long k) { return (int) (k >> 32); }
inline int qj(long long k) { return (int) (k & 0xffffffff); }

template <class T> inline void merge_into(T& r, const T& x, double sx, const T& y, double sy) {
	r.clear(); r.reserve(x.size() + y.size());
	size_t i = 0, j = 0;
	while (i < x.size() || j < y.size()) {
		if (j == y.size() || (i < x.size() && x[i].first < y[j].first)) { r.push_back(std::make_pair(x[i].first, sx * x[i].second)); i++; }
		else if (i == x.size() || y[j].first < x[i].first) { r.push_back(std::make_pair(y[j].first, sy * y[j].second)); j++; }
		else { r.push_back(std::make_pair(x[i].first, sx * x[i].second + sy * y[j].second)); i++; j++; }
	}
}
template <class T> inline T merge(const T& x, double sx, const T& y, double sy) { T r; merge_into(r, x, sx, y, sy); return r; }

struct QF { double c; SV lin; SQ Q; bool ok; QF() : c(0), ok(true) {} };

inline void quad_interval(const QF& X, double& lo, double& hi) {
	lo = hi = 0;
	for (size_t k = 0; k < X.Q.size(); k++) {
		const double v = X.Q[k].second;
		if (qi(X.Q[k].first) == qj(X.Q[k].first)) { lo += std::min(0.0, v); hi += std::max(0.0, v); }
		else { lo -= std::fabs(v); hi += std::fabs(v); }
	}
}
inline double lin_abs(const QF& X) { double s = 0; for (size_t k = 0; k < X.lin.size(); k++) s += std::fabs(X.lin[k].second); return s; }

inline void range_one(double a, double q, double& lo, double& hi) {
	double vmin = std::min(q + a, q - a), vmax = std::max(q + a, q - a);
	if (q != 0 && std::fabs(a / (2 * q)) < 1) { const double s = -a * a / (4 * q); vmin = std::min(vmin, s); vmax = std::max(vmax, s); }
	lo += vmin; hi += vmax;
}

// rango (cota) de c + a^T z + z^T Q z: cada simbolo exacto con su lineal y su diagonal, resto termino a termino.
// lin esta ordenada por indice y Q por clave (i, j), asi que las entradas diagonales (i, i) aparecen en orden de i:
// basta un recorrido doble (sin estructuras auxiliares).
inline void qf_range(const QF& X, double& lo, double& hi, bool with_c = true) {
	lo = hi = with_c ? X.c : 0.0;
	double off = 0;
	size_t a = 0;                                   // puntero en lin
	for (size_t k = 0; k < X.Q.size(); k++) {
		const int i = qi(X.Q[k].first);
		if (i != qj(X.Q[k].first)) { off += std::fabs(X.Q[k].second); continue; }
		while (a < X.lin.size() && X.lin[a].first < i) { range_one(X.lin[a].second, 0.0, lo, hi); a++; }
		if (a < X.lin.size() && X.lin[a].first == i) { range_one(X.lin[a].second, X.Q[k].second, lo, hi); a++; }
		else range_one(0.0, X.Q[k].second, lo, hi);
	}
	for (; a < X.lin.size(); a++) range_one(X.lin[a].second, 0.0, lo, hi);
	lo -= off; hi += off;
}

inline void qf_scale_into(QF& r, const QF& X, double s) {
	r.c = s * X.c;
	r.lin.resize(X.lin.size()); for (size_t k = 0; k < X.lin.size(); k++) r.lin[k] = std::make_pair(X.lin[k].first, s * X.lin[k].second);
	r.Q.resize(X.Q.size()); for (size_t k = 0; k < X.Q.size(); k++) r.Q[k] = std::make_pair(X.Q[k].first, s * X.Q[k].second);
}
inline QF qf_scale(const QF& X, double s) { QF r; qf_scale_into(r, X, s); return r; }

// ------------------------------------------------------------------ linealizacion de funciones
struct Lin { double alpha, d1, d2; bool ok; };
inline Lin make_lin(double a, double d1, double d2) { Lin L; L.alpha = a; L.d1 = d1; L.d2 = d2; L.ok = true; return L; }

inline double fval(int op, double x, double k);
enum Op { V, C, CI, ADD, SUB, MUL, SCL, NEG, REC, POW, EXP, LOG, SQRT, SIN, COS, ABS, IFUN };

inline Lin crit_lin(int op, double k, double lo, double hi, double a, const std::vector<double>& crit) {
	double d1 = INFINITY, d2 = -INFINITY;
	const double ts[2] = {lo, hi};
	for (int i = 0; i < 2; i++) { const double d = fval(op, ts[i], k) - a * ts[i]; d1 = std::min(d1, d); d2 = std::max(d2, d); }
	for (size_t i = 0; i < crit.size(); i++) if (crit[i] > lo && crit[i] < hi) {
		const double d = fval(op, crit[i], k) - a * crit[i]; d1 = std::min(d1, d); d2 = std::max(d2, d); }
	return make_lin(a, d1, d2);
}

inline double fval(int op, double x, double k) {
	switch (op) {
	case EXP: return std::exp(x); case LOG: return std::log(x); case SQRT: return std::sqrt(x);
	case REC: return 1 / x; case SIN: return std::sin(x); case COS: return std::cos(x);
	case ABS: return std::fabs(x); default: return std::pow(x, k);
	}
}

inline Lin linearize(int op, double lo, double hi, double k) {
	Lin bad; bad.ok = false; bad.alpha = bad.d1 = bad.d2 = 0;
	// exp, log, sqrt: Chebyshev (pendiente de la secante, residuo en los extremos y en el punto critico), como el
	// modo por omision de ibex-affine (AF_Lin_Chebyshev)
	if (op == EXP) {
		const double a = (std::exp(hi) - std::exp(lo)) / (hi - lo);
		return crit_lin(EXP, 0, lo, hi, a, std::vector<double>(1, a > 0 ? std::log(a) : lo));
	}
	if (op == LOG || op == SQRT) {
		if (hi <= 0) return bad;
		lo = std::max(lo, 1e-300);
		const double a = (fval(op, hi, 0) - fval(op, lo, 0)) / (hi - lo);
		if (!(a > 0) || !std::isfinite(a)) return bad;
		const double u = op == LOG ? 1 / a : 1 / (4 * a * a);       // f'(u) = a
		return crit_lin(op, 0, lo, hi, a, std::vector<double>(1, u));
	}
	if (op == REC) {
		if (lo <= 0 && hi >= 0) return bad;
		const double a = -1 / (lo * hi);
		std::vector<double> crit(1, std::copysign(std::sqrt(lo * hi), lo));
		return crit_lin(REC, 0, lo, hi, a, crit);
	}
	if (op == ABS) {
		if (lo >= 0) return make_lin(1, 0, 0);
		if (hi <= 0) return make_lin(-1, 0, 0);
		return crit_lin(ABS, 0, lo, hi, (std::fabs(hi) - std::fabs(lo)) / (hi - lo), std::vector<double>(1, 0.0));
	}
	if (op == SIN || op == COS) {
		const double a = (fval(op, hi, 0) - fval(op, lo, 0)) / (hi - lo);
		std::vector<double> crit;
		if (std::fabs(a) <= 1) {
			const bool s = op == SIN;
			const double b1 = s ? std::acos(a) : std::asin(-a), b2 = s ? -std::acos(a) : M_PI - std::asin(-a);
			const double n0 = std::floor((lo - 4 * M_PI) / (2 * M_PI)), n1 = std::ceil((hi + 4 * M_PI) / (2 * M_PI));
			if (n1 - n0 > 400) {
				const double m1 = std::min(a * lo, a * hi), m2 = std::max(a * lo, a * hi);
				return make_lin(a, -1 - m2, 1 - m1);
			}
			for (double n = n0; n <= n1; n++) { crit.push_back(b1 + 2 * M_PI * n); crit.push_back(b2 + 2 * M_PI * n); }
		}
		return crit_lin(op, 0, lo, hi, a, crit);
	}
	if (op == POW) {
		const double a = (std::pow(hi, k) - std::pow(lo, k)) / (hi - lo);
		std::vector<double> crit;
		const double r = a / k; const int e = (int) k - 1;
		if (e % 2 == 1) crit.push_back(std::copysign(std::pow(std::fabs(r), 1.0 / e), r));
		else if (r >= 0) { const double t = std::pow(r, 1.0 / e); crit.push_back(t); crit.push_back(-t); }
		return crit_lin(POW, k, lo, hi, a, crit);
	}
	return bad;
}

// ------------------------------------------------------------------ aproximacion cuadratica de funciones
// Modo de las funciones unarias: 0 = lineal (Chebyshev), 1 = cuadratica cuando mejora (variable Q2_UNARY=lin|quad).
inline int q2_unary_mode() {
	// 1 = cuadratica con cota de Taylor rigurosa; 2 = cuadratica con el error medido por muestreo denso (NO rigurosa,
	// solo para medir el techo de la idea)
	static int m = !getenv("Q2_UNARY") ? 0 : std::string(getenv("Q2_UNARY")) == "quad" ? 1
	             : std::string(getenv("Q2_UNARY")) == "quadmuestreo" ? 2 : 0;
	return m;
}
// Cota superior de |f^(3)| en dom (aritmetica de intervalos de Ibex, redondeo hacia fuera). +inf si no aplica.
inline double third_deriv_bound(int op, const Interval& dom, double k) {
	Interval t;
	switch (op) {
	case EXP: t = exp(dom); break;
	case LOG: if (dom.lb() <= 0) return INFINITY; t = 2.0 / pow(dom, 3); break;
	case SQRT: if (dom.lb() <= 0) return INFINITY; t = 0.375 / (pow(dom, 2) * sqrt(dom)); break;
	case REC: if (dom.contains(0)) return INFINITY; t = 6.0 / pow(dom, 4); break;
	case SIN: t = cos(dom); break;                    // f^(3) = -cos
	case COS: t = sin(dom); break;                    // f^(3) = sin
	case POW: { const int kk = (int) k; if (kk < 3) return 0.0; t = Interval((double) kk * (kk - 1) * (kk - 2)) * pow(dom, kk - 3); break; }
	default: return INFINITY;                         // abs (no derivable), otras
	}
	return t.is_empty() ? INFINITY : t.mag();
}

// ------------------------------------------------------------------ cinta desde una ibex::Function
struct Node { int op, a, b; double c, d; const ExprNode* e; };   // d: cota superior (CI); e: nodo de Ibex (IFUN)

struct Tape {
	std::vector<Node> nodes;
	bool supported;
	Tape() : supported(true) {}

	static Tape build(const Function& f) {
		Tape t;
		if (!f.expr().dim.is_scalar()) { t.supported = false; return t; }
		std::map<const ExprNode*, int> id;
		// primer indice de variable de cada argumento (Function::symbol_index es protegido)
		std::vector<int> off(f.nb_arg() + 1, 0);
		for (int i = 0; i < f.nb_arg(); i++) off[i + 1] = off[i] + f.args()[i].dim.size();
		for (int k = f.nodes.size() - 1; k >= 0; k--) {
			const ExprNode& e = f.nodes[k];
			Node n; n.op = -1; n.a = n.b = -1; n.c = n.d = 0; n.e = &e;
			if (const ExprSymbol* s = dynamic_cast<const ExprSymbol*>(&e)) {
				if (!s->dim.is_scalar()) continue;            // vector: se usa via ExprIndex
				n.op = V; n.a = off[s->key];
			} else if (const ExprIndex* ix = dynamic_cast<const ExprIndex*>(&e)) {
				const ExprSymbol* s = dynamic_cast<const ExprSymbol*>(&ix->expr);
				if (!s || !ix->dim.is_scalar() || !s->dim.is_vector()) { t.supported = false; return t; }
				n.op = V; n.a = off[s->key] + ix->index.first_row();
			} else if (const ExprConstant* c = dynamic_cast<const ExprConstant*>(&e)) {
				if (!c->dim.is_scalar()) { t.supported = false; return t; }
				const Interval v = c->get_value();
				if (v.is_degenerated()) { n.op = C; n.c = v.mid(); } else { n.op = CI; n.c = v.lb(); n.d = v.ub(); }
			} else if (const ExprBinaryOp* b = dynamic_cast<const ExprBinaryOp*>(&e)) {
				if (!b->dim.is_scalar()) { t.supported = false; return t; }
				n.a = id.at(&b->left); n.b = id.at(&b->right);
				if (dynamic_cast<const ExprAdd*>(&e)) n.op = ADD;
				else if (dynamic_cast<const ExprSub*>(&e)) n.op = SUB;
				else if (dynamic_cast<const ExprMul*>(&e)) n.op = MUL;
				else if (dynamic_cast<const ExprDiv*>(&e)) {
					Node r; r.op = REC; r.a = n.b; r.b = -1; r.c = r.d = 0; r.e = NULL;
					t.nodes.push_back(r); n.op = MUL; n.b = (int) t.nodes.size() - 1;
				} else { t.supported = false; return t; }      // max, min, atan2: no soportadas
			} else if (const ExprUnaryOp* u = dynamic_cast<const ExprUnaryOp*>(&e)) {
				if (!u->dim.is_scalar()) { t.supported = false; return t; }
				n.a = id.at(&u->expr);
				if (dynamic_cast<const ExprMinus*>(&e)) n.op = NEG;
				else if (dynamic_cast<const ExprSqr*>(&e)) { n.op = MUL; n.b = n.a; }
				else if (const ExprPower* p = dynamic_cast<const ExprPower*>(&e)) {
					if (p->expon == 2) { n.op = MUL; n.b = n.a; }
					else if (p->expon == 1) { id[&e] = n.a; continue; }
					else if (p->expon > 2) { n.op = POW; n.c = p->expon; }
					else n.op = IFUN;
				}
				else if (dynamic_cast<const ExprExp*>(&e)) n.op = EXP;
				else if (dynamic_cast<const ExprLog*>(&e)) n.op = LOG;
				else if (dynamic_cast<const ExprSqrt*>(&e)) n.op = SQRT;
				else if (dynamic_cast<const ExprSin*>(&e)) n.op = SIN;
				else if (dynamic_cast<const ExprCos*>(&e)) n.op = COS;
				else if (dynamic_cast<const ExprAbs*>(&e)) n.op = ABS;
				else n.op = IFUN;                           // tan, atan, cosh, ...: por intervalos
			} else { t.supported = false; return t; }       // vectores, apply, chi
			t.nodes.push_back(n);
			id[&e] = (int) t.nodes.size() - 1;
		}
		t.root = id.at(&f.expr());
		return t;
	}
	int root;
};

// evaluacion por intervalos de una operacion unaria de Ibex sobre [lo, hi] (para IFUN)
inline Interval eval_unary_interval(const ExprNode* e, const Interval& x) {
	if (const ExprPower* p = dynamic_cast<const ExprPower*>(e)) return pow(x, p->expon);
	if (dynamic_cast<const ExprTan*>(e)) return tan(x);
	if (dynamic_cast<const ExprAtan*>(e)) return atan(x);
	if (dynamic_cast<const ExprAcos*>(e)) return acos(x);
	if (dynamic_cast<const ExprAsin*>(e)) return asin(x);
	if (dynamic_cast<const ExprCosh*>(e)) return cosh(x);
	if (dynamic_cast<const ExprSinh*>(e)) return sinh(x);
	if (dynamic_cast<const ExprTanh*>(e)) return tanh(x);
	if (dynamic_cast<const ExprAcosh*>(e)) return acosh(x);
	if (dynamic_cast<const ExprAsinh*>(e)) return asinh(x);
	if (dynamic_cast<const ExprAtanh*>(e)) return atanh(x);
	if (dynamic_cast<const ExprSign*>(e)) return sign(x);
	if (dynamic_cast<const ExprFloor*>(e)) return floor(x);
	if (dynamic_cast<const ExprCeil*>(e)) return ceil(x);
	return Interval::all_reals();
}

// ------------------------------------------------------------------ evaluador
struct Eval {
	int nvar, nsym;
	// memoria que se reutiliza entre evaluaciones (sin reservas en el bucle una vez estabilizada la capacidad)
	std::vector<QF> v; std::vector<Interval> d; std::vector<double> M;
	SQ outer, oc, qtmp;
	explicit Eval(int nvar) : nvar(nvar), nsym(nvar) {}
	void add_err(QF& F, double e) { if (e > 0) F.lin.push_back(std::make_pair(nsym++, e)); }

	// r = X * Y (r no debe ser X ni Y)
	void mul_into(QF& r, const QF& X, const QF& Y, bool same) {
		double qx0, qx1, qy0, qy1;
		quad_interval(X, qx0, qx1); quad_interval(Y, qy0, qy1);
		const double mx = 0.5 * (qx0 + qx1), my = 0.5 * (qy0 + qy1), rx = 0.5 * (qx1 - qx0), ry = 0.5 * (qy1 - qy0);
		const double Cx = X.c + mx, Cy = Y.c + my;
		merge_into(r.lin, X.lin, Cy, Y.lin, Cx);
		merge_into(qtmp, X.Q, Cy, Y.Q, Cx);
		outer.clear();
		for (size_t i = 0; i < X.lin.size(); i++) for (size_t j = 0; j < Y.lin.size(); j++)
			outer.push_back(std::make_pair(qk(X.lin[i].first, Y.lin[j].first), X.lin[i].second * Y.lin[j].second));
		std::sort(outer.begin(), outer.end());
		oc.clear();
		for (size_t k = 0; k < outer.size(); k++) {
			if (!oc.empty() && oc.back().first == outer[k].first) oc.back().second += outer[k].second; else oc.push_back(outer[k]);
		}
		merge_into(r.Q, qtmp, 1.0, oc, 1.0);
		const double sa = lin_abs(X), sb = lin_abs(Y);
		double lo, hi;
		if (same) { lo = -2 * sa * rx; hi = 2 * sa * rx + rx * rx; }
		else { const double t = sa * ry + sb * rx + rx * ry; lo = -t; hi = t; }
		r.c = Cx * Cy - Cx * my - Cy * mx + 0.5 * (lo + hi);
		add_err(r, 0.5 * (hi - lo));
	}
	QF mul(const QF& X, const QF& Y, bool same) { QF r; mul_into(r, X, Y, same); return r; }

	// Evaluacion de la cinta en la caja. Como AffineEval de ibex-affine, en cada nodo se lleva tambien un intervalo
	// d = (operacion por intervalos sobre los intervalos de los hijos) & (rango de la forma, ensanchado por un margen
	// QUAD2_MARGIN relativo a la magnitud acumulada del calculo): ese intervalo es el dominio de las linealizaciones
	// de funciones y su valor en la raiz es la evaluacion 'ev' que usa el linearizador. Devuelve false si la
	// evaluacion no es posible en esta caja (dominio, valores no finitos, interseccion vacia).
	// mag: magnitud acumulada del calculo (suma, sobre todos los nodos del subarbol, del tamano de sus terminos);
	// acota la escala de los errores de redondeo, tambien cuando hay cancelacion entre operandos grandes.
	// El resultado queda en result() (sin copia); out recibe una copia si se pide.
	// Aproximacion cuadratica f(X) ~ f(m) + alpha (X-m) + gamma (X-m)^2 +- delta, interpolando f en los tres nodos
	// de Chebyshev de dom = [m-rr, m+rr] (m y m +- s, s = rr sqrt(3)/2). Cota del error de interpolacion:
	// |f - p| <= max|f^(3)| rr^3 / 24  (max |(x-m)((x-m)^2 - s^2)| = rr^3/4 en el dominio). (X-m)^2 se calcula con
	// mul_into (grado 2 conservado, resto de grado 3-4 a un simbolo). Se usa solo si delta + |gamma| (resto de
	// (X-m)^2) < dlin (error de la linealizacion de Chebyshev); si no, devuelve false y se usa la lineal.
	QF ybuf, ysq, ay, gs;
	bool unary_quad(QF& r, const QF& x, const Node& N, const Interval& dom, double dlin) {
		const double m = dom.mid(), rr = dom.rad();
		if (!(rr > 0) || !std::isfinite(m)) return false;
		const double sn = rr * 0.8660254037844386;
		const double f0 = fval(N.op, m, N.c), fp = fval(N.op, m + sn, N.c), fm = fval(N.op, m - sn, N.c);
		if (!(std::isfinite(f0) && std::isfinite(fp) && std::isfinite(fm))) return false;
		const double alpha = (fp - fm) / (2 * sn), gamma = (fp + fm - 2 * f0) / (2 * sn * sn);
		double dint;
		if (q2_unary_mode() == 2) {          // muestreo denso del residuo (no riguroso)
			dint = 0;
			for (int i = 0; i <= 2000; i++) {
				const double t = dom.lb() + dom.diam() * i / 2000.0, y = t - m;
				const double fv = fval(N.op, t, N.c);
				if (!std::isfinite(fv)) return false;
				dint = std::max(dint, std::fabs(fv - (f0 + alpha * y + gamma * y * y)));
			}
			dint *= 1.01;
		} else {
			const double M3 = third_deriv_bound(N.op, dom, N.c);
			if (!std::isfinite(M3)) return false;
			dint = M3 * rr * rr * rr / 24.0;
		}
		if (!(dint < dlin)) return false;
		ybuf.c = x.c - m; ybuf.lin = x.lin; ybuf.Q = x.Q;
		const int nsym0 = nsym;
		mul_into(ysq, ybuf, ybuf, true);
		const double eS = (nsym > nsym0) ? ysq.lin.back().second : 0.0;   // error de truncar (X-m)^2
		const double dtot = dint + std::fabs(gamma) * eS;
		if (!(dtot < dlin)) { nsym = nsym0; return false; }
		qf_scale_into(ay, ybuf, alpha);
		qf_scale_into(gs, ysq, gamma);
		r.c = f0 + ay.c + gs.c;
		merge_into(r.lin, ay.lin, 1.0, gs.lin, 1.0);
		merge_into(r.Q, ay.Q, 1.0, gs.Q, 1.0);
		add_err(r, dint);
		quad_unary_count()++;
		return true;
	}
	static long& quad_unary_count() { static long c = 0; return c; }

	const QF& result() const { return v[root_]; }
	bool eval(const Tape& T, const IntervalVector& box, QF& out, Interval& ev) { double mag; return eval(T, box, out, ev, mag); }
	bool eval(const Tape& T, const IntervalVector& box, QF& out, Interval& ev, double& mag_out) {
		if (!eval(T, box, ev, mag_out)) return false;
		out = v[root_];
		return true;
	}
	bool eval(const Tape& T, const IntervalVector& box, Interval& ev, double& mag_out) {
		nsym = nvar;
		const size_t nn = T.nodes.size();
		if (v.size() < nn) { v.resize(nn); d.resize(nn); M.resize(nn); }
		for (size_t k = 0; k < nn; k++) {
			const Node& N = T.nodes[k];
			QF& r = v[k];
			r.c = 0; r.lin.clear(); r.Q.clear();
			Interval di;
			switch (N.op) {
			case V: {
				r.c = box[N.a].mid(); const double rd = box[N.a].rad();
				if (!(std::isfinite(r.c) && std::isfinite(rd))) return false;
				if (rd > 0) r.lin.push_back(std::make_pair(N.a, rd));
				di = box[N.a];
				break; }
			case C: r.c = N.c; di = Interval(N.c); break;
			case CI: r.c = 0.5 * (N.c + N.d); add_err(r, 0.5 * (N.d - N.c)); di = Interval(N.c, N.d); break;
			case ADD: case SUB: {
				const double sg = N.op == ADD ? 1 : -1;
				r.c = v[N.a].c + sg * v[N.b].c; merge_into(r.lin, v[N.a].lin, 1.0, v[N.b].lin, sg); merge_into(r.Q, v[N.a].Q, 1.0, v[N.b].Q, sg);
				di = N.op == ADD ? d[N.a] + d[N.b] : d[N.a] - d[N.b];
				break; }
			case NEG: qf_scale_into(r, v[N.a], -1.0); di = -d[N.a]; break;
			case MUL: mul_into(r, v[N.a], v[N.b], N.a == N.b); di = N.a == N.b ? sqr(d[N.a]) : d[N.a] * d[N.b]; break;
			default: {
				const QF& x = v[N.a];
				const Interval& dx = d[N.a];
				const double l = dx.lb(), h = dx.ub();
				if (!(std::isfinite(l) && std::isfinite(h))) return false;
				switch (N.op) {
				case REC: di = 1.0 / dx; break;
				case POW: di = pow(dx, (int) N.c); break;
				case EXP: di = exp(dx); break;
				case LOG: di = log(dx); break;
				case SQRT: di = sqrt(dx); break;
				case SIN: di = sin(dx); break;
				case COS: di = cos(dx); break;
				case ABS: di = abs(dx); break;
				default: di = eval_unary_interval(N.e, dx);
				}
				if (di.is_empty() || di.is_unbounded()) return false;
				if (N.op == IFUN) {
					r.c = di.mid(); add_err(r, di.rad());
				} else if (h - l <= 1e-15 * std::max(1.0, std::fabs(l))) {
					if ((N.op == REC && x.c == 0) || ((N.op == LOG || N.op == SQRT) && x.c <= 0)) return false;
					r.c = fval(N.op, 0.5 * (l + h), N.c);
				} else {
					const Lin L = linearize(N.op, l, h, N.c);
					if (!L.ok) return false;
					const double dlin = 0.5 * (L.d2 - L.d1);
					bool done = false;
					if (q2_unary_mode() >= 1 && N.op != ABS) done = unary_quad(r, x, N, dx, dlin);
					if (!done) {
						qf_scale_into(r, x, L.alpha); r.c += 0.5 * (L.d1 + L.d2);
						add_err(r, dlin);
					}
				}
			} }
			if (!std::isfinite(r.c)) return false;
			double fl, fh; qf_range(r, fl, fh);
			if (!(std::isfinite(fl) && std::isfinite(fh))) return false;
			// margen relativo a la magnitud acumulada (no al rango): con cancelacion (p. ej. (x-y)^2 + ... - 25 con
			// valor cerca de 0) el error de redondeo es del orden de los operandos, no del resultado
			double mag = std::fabs(r.c);
			for (size_t q = 0; q < r.lin.size(); q++) mag += std::fabs(r.lin[q].second);
			for (size_t q = 0; q < r.Q.size(); q++) mag += std::fabs(r.Q[q].second);
			if (N.op != V && N.op != C && N.op != CI) {
				mag += M[N.a];
				if (N.b >= 0 && N.op != REC) mag += M[N.b];
				mag += std::fabs(d[N.a].mag()) + (N.b >= 0 ? std::fabs(d[N.b].mag()) : 0.0);
			}
			M[k] = mag;
			const double mg = QUAD2_MARGIN * mag + 1e-300;
			di &= Interval(fl - mg, fh + mg);
			if (di.is_empty()) return false;
			d[k] = di;
		}
		root_ = T.root;
		ev = d[T.root];
		mag_out = M[T.root];
		return true;
	}
private:
	int root_ = 0;
};

// ------------------------------------------------------------------ encierro final del resto
// homogeneizada de c? no: de r(z) = lin^T z + z^T Q z  (sin constante); devuelve el tamano m+1
inline int homog(const SV& lin, const SQ& Q, Vec& H) {
	static thread_local std::vector<int> ids;
	ids.clear();
	for (size_t k = 0; k < lin.size(); k++) ids.push_back(lin[k].first);
	for (size_t k = 0; k < Q.size(); k++) { ids.push_back(qi(Q[k].first)); ids.push_back(qj(Q[k].first)); }
	std::sort(ids.begin(), ids.end()); ids.erase(std::unique(ids.begin(), ids.end()), ids.end());
	const int m = (int) ids.size(), s = m + 1;
	auto ix = [&](int id) { return (int) (std::lower_bound(ids.begin(), ids.end(), id) - ids.begin()); };
	H.assign(s * s, 0.0);
	for (size_t k = 0; k < Q.size(); k++) {
		const int a = ix(qi(Q[k].first)), b = ix(qj(Q[k].first));
		if (a == b) H[a * s + a] += Q[k].second; else { H[a * s + b] += 0.5 * Q[k].second; H[b * s + a] += 0.5 * Q[k].second; }
	}
	for (size_t k = 0; k < lin.size(); k++) { const int a = ix(lin[k].first); H[a * s + m] += 0.5 * lin[k].second; H[m * s + a] += 0.5 * lin[k].second; }
	return s;
}

// k-nomial: cota [L, U] de z^T H z en [-1,1]^n (como bound_full del driver de medicion, redondeo hacia +inf)
// k del k-nomial de q2kn: por omision 3 (en las mediciones en la raiz k >= 4 casi no cambia la cota y la busqueda
// de cliques crece mucho); se puede cambiar con la variable de entorno Q2_KN_K.
inline int q2_kn_k() { static int k = getenv("Q2_KN_K") ? std::atoi(getenv("Q2_KN_K")) : 3; return k; }

inline void knomial_interval(const Vec& H0, int n, int k, double& L, double& U) {
	static thread_local Vec Q, w;
	static thread_local std::vector<KNomialBlock> blocks;
	Q = H0; w.assign(n, 2.0); blocks.clear();
	knomial_rewrite(Q, w, n, k, blocks);
	const int old = std::fegetround(); std::fesetround(FE_UPWARD);
	double u = 0, negl = 0;
	for (size_t b = 0; b < blocks.size(); b++) {
		double s = 0; for (int i = 0; i < n; i++) s += std::fabs(blocks[b].u[i]);
		if (blocks[b].c > 0) u += s * s; else negl += s * s;
	}
	for (int i = 0; i < n; i++) for (int j = i; j < n; j++) {
		const double q = 0.5 * (H0[i * n + j] + H0[j * n + i]);
		double tu = 0, nt = 0;
		for (size_t b = 0; b < blocks.size(); b++) if (blocks[b].u[i] != 0 && blocks[b].u[j] != 0) {
			const double si = blocks[b].c * blocks[b].u[i]; tu += si * blocks[b].u[j]; nt += (-si) * blocks[b].u[j]; }
		const double rup = q + nt, rlo = -(tu - q);
		if (i == j) { if (rup > 0) u += rup; if (rlo < 0) negl += -rlo; }
		else { const double m2 = 2.0 * std::max(std::fabs(rlo), std::fabs(rup)); u += m2; negl += m2; }
	}
	std::fesetround(old);
	L = -negl; U = u;
}

// cota [L, U] del resto r(z) = sum_{simbolos e >= nvar} a_e z_e + z^T Q z (la parte lineal en las variables queda fuera)
inline void remainder_bound(const QF& F, int nvar, FinalBound fb, double& L, double& U) {
	static thread_local QF R; static thread_local Vec H;
	R.c = 0; R.Q = F.Q; R.lin.clear();
	for (size_t k = 0; k < F.lin.size(); k++) if (F.lin[k].first >= nvar) R.lin.push_back(F.lin[k]);
	qf_range(R, L, U, false);
	if (fb == FB_TT || R.Q.empty()) return;
	const int s = homog(R.lin, R.Q, H);
	double l2, u2;
	if (fb == FB_KN) knomial_interval(H, s, q2_kn_k(), l2, u2);
	else {
		const prodk::Result res = prodk::prodk_interval(H, s, true, fb == FB_FULL ? 40 : 0, 12);
		l2 = res.L; u2 = res.U;
	}
	L = std::max(L, l2); U = std::min(U, u2);
}

} // namespace quad2
} // namespace ibex

#endif
