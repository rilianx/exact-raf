//============================================================================
// LinearizerCount: envoltorio que cuenta las llamadas a linearize y las filas agregadas (para todas las
// aritmeticas por igual). Cada llamada corresponde a una construccion del politopo en CtcPolytopeHull, que luego
// resuelve hasta 2n LPs; asi se distingue "mas LPs" de "LPs mas caros".
//============================================================================
#ifndef __IBEX_LINEARIZER_COUNT_H__
#define __IBEX_LINEARIZER_COUNT_H__
#include "ibex_Linearizer.h"
#include <chrono>
namespace ibex {
struct LinCountStats { long calls, rows, infeasible; double seconds; LinCountStats() : calls(0), rows(0), infeasible(0), seconds(0) {} };
inline LinCountStats& lin_count_stats() { static LinCountStats s; return s; }
class LinearizerCount : public Linearizer {
public:
	explicit LinearizerCount(Linearizer& inner) : Linearizer(inner.nb_var()), inner(inner) {}
	virtual int linearize(const IntervalVector& box, LPSolver& lp) {
		const auto t0 = std::chrono::steady_clock::now();
		const int r = inner.linearize(box, lp);
		lin_count_stats().seconds += std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
		lin_count_stats().calls++;
		if (r < 0) lin_count_stats().infeasible++; else lin_count_stats().rows += r;
		return r;
	}
private:
	Linearizer& inner;
};
}
#endif
