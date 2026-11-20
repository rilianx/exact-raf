//============================================================================
// LinearizerQuad2: relajacion lineal con formas afines de grado 2 (ibex_Quad2.h)
//
// Misma estructura que LinearizerAffineT<T>::linearize (ibex-affine, parches raak): para cada restriccion
// g(x) op 0 se evalua la forma de grado 2 en la caja; la parte lineal en las variables da la fila del LP y el
// resto (parte cuadratica + simbolos de error) se acota con el encierro final elegido; esa cota hace de error.
// Las comprobaciones con la evaluacion por intervalos ev (infactibilidad, restriccion inactiva, igualdad
// estrecha) son las mismas. Si la restriccion no tiene cinta (expresion no soportada) o la evaluacion falla en
// la caja, la restriccion no aporta fila (como una forma fAF2 no valida).
//============================================================================
#ifndef __IBEX_LINEARIZER_QUAD2_H__
#define __IBEX_LINEARIZER_QUAD2_H__

#include "ibex_Linearizer.h"
#include "ibex_System.h"
#include "ibex_LPSolver.h"
#include "ibex_Quad2.h"
#include <cmath>
#include <vector>

namespace ibex {

struct Quad2Stats { long evals, failed, unsupported, rows; Quad2Stats() : evals(0), failed(0), unsupported(0), rows(0) {} };
inline Quad2Stats& quad2_stats() { static Quad2Stats s; return s; }

class LinearizerQuad2 : public Linearizer {
public:
	LinearizerQuad2(const System& sys1, quad2::FinalBound fb) : Linearizer(sys1.nb_var), sys(sys1), fb(fb) {
		for (int i = 0; i < sys.nb_ctr; i++) {
			tapes.push_back(quad2::Tape::build(sys.ctrs[i].f));
			if (!tapes.back().supported) quad2_stats().unsupported++;
		}
	}

	virtual int linearize(const IntervalVector& box, LPSolver& lp_solver) {
		const int n = sys.nb_var;
		Vector row(n);
		int cont = 0;
		for (int ctr = 0; ctr < sys.nb_ctr; ctr++) {
			const CmpOp op = sys.ctrs[ctr].op;
			if (!tapes[ctr].supported) continue;
			// evaluacion de grado 2; ev = intervalo natural & intervalos por nodo intersectados con el rango de la forma
			// (como el ev de AffineEval en LinearizerAffineT)
			quad2::QF F; quad2::Eval E(n);
			Interval evq; double mag;
			quad2_stats().evals++;
			if (!E.eval(tapes[ctr], box, F, evq, mag)) { quad2_stats().failed++; continue; }
			const Interval ev = sys.ctrs[ctr].f.eval(box) & evq;
			if (ev.is_empty()) continue;
			// mismas salidas por infactibilidad que LinearizerAffineT
			if (op == LEQ && 0.0 < ev.lb()) return -1;
			if (op == GEQ && ev.ub() < 0.0) return -1;
			if (op == EQ && !ev.contains(0.0)) return -1;
			if (op == LT && ev.lb() == 0.0) return -1;
			if (op == GT && ev.ub() == 0.0) return -1;
			if (op == LT || op == GT) continue;
			if (op == LEQ && !(0.0 < ev.ub())) continue;
			if (op == GEQ && !(ev.lb() < 0.0)) continue;
			if (op == EQ && !(ev.diam() > 2 * lp_solver.tolerance())) continue;

			// fila: coeficientes a_j / r_j de las variables; constante c - sum (a_j / r_j) m_j
			std::vector<double> a(n, 0.0);
			for (size_t k = 0; k < F.lin.size(); k++) if (F.lin[k].first < n) a[F.lin[k].first] = F.lin[k].second;
			bool abort = false;
			Interval center(0.0); double err = 0.0, scale = std::fabs(F.c);
			for (int i = 0; i < n && !abort; i++) {
				const double r = box[i].rad();
				if (r == 0) { if (a[i] != 0) abort = true; row[i] = 0; }
				else {
					row[i] = a[i] / r;
					center += row[i] * box[i].mid();
					err += std::fabs(row[i]) * std::pow(2, -50);
					scale += std::fabs(a[i]);
				}
			}
			if (abort) continue;
			double L, U;
			quad2::remainder_bound(F, n, fb, L, U);
			if (!(std::isfinite(L) && std::isfinite(U))) { quad2_stats().failed++; continue; }
			// margen de seguridad (no riguroso): relativo a la magnitud de la forma
			const double margin = quad2::QUAD2_MARGIN * (scale + mag + std::fabs(L) + std::fabs(U)) + 1e-300;
			// g(x) en [c - center + row.x + L, c - center + row.x + U]
			const Interval cst = Interval(F.c) - center;
			try {
				if (op == LEQ || op == EQ) {        // row.x <= -(c - center) - L
					lp_solver.add_constraint(row, LEQ, (-(cst) - Interval(L) + Interval(err + margin)).ub());
					cont++; quad2_stats().rows++;
				}
				if (op == GEQ || op == EQ) {        // row.x >= -(c - center) - U
					lp_solver.add_constraint(row, GEQ, (-(cst) - Interval(U) - Interval(err + margin)).lb());
					cont++; quad2_stats().rows++;
				}
			} catch (LPException&) { }
		}
		return cont;
	}

private:
	const System& sys;
	quad2::FinalBound fb;
	std::vector<quad2::Tape> tapes;
};

} // namespace ibex

#endif
