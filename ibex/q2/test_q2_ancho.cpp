// Ancho del encierro (intervalo por nodo en la raiz, 'ev' del linearizador) de cada restriccion de un sistema, en
// cajas aleatorias dentro de la caja inicial (escalas 1, 1e-1, 1e-2, 1e-3; semilla fija). Sirve para comparar
// modos de Quad2 (p. ej. Q2_UNARY=lin|quad) en las mismas cajas. Salida: una linea por (caja, restriccion):
//   escala ctr ancho
// Uso: test_q2_ancho <sistema.bch> [cajas por escala=50]
#include "ibex.h"
#include "ibex_Quad2.h"
#include <cstdio>
#include <random>
using namespace ibex;
int main(int argc, char** argv) {
	System sys(argv[1]);
	const int nb = argc > 2 ? std::atoi(argv[2]) : 50;
	const int n = sys.nb_var;
	std::vector<quad2::Tape> tapes;
	for (int c = 0; c < sys.nb_ctr; c++) tapes.push_back(quad2::Tape::build(sys.ctrs[c].f));
	quad2::Eval E(n);
	std::mt19937 rng(7);
	std::uniform_real_distribution<double> U01(0, 1);
	const double scales[4] = {1, 1e-1, 1e-2, 1e-3};
	for (int si = 0; si < 4; si++)
		for (int b = 0; b < nb; b++) {
			IntervalVector box(n);
			for (int i = 0; i < n; i++) {
				Interval B = sys.box[i];
				if (B.is_unbounded()) B = Interval(-1e8, 1e8) & B;
				const double w = B.diam() * scales[si];
				const double lo = B.lb() + (B.diam() - w) * U01(rng);
				box[i] = Interval(lo, lo + w) & B;
			}
			for (int c = 0; c < sys.nb_ctr; c++) {
				if (!tapes[c].supported) continue;
				Interval ev; double mag;
				if (!E.eval(tapes[c], box, ev, mag)) continue;
				std::printf("%g %d %.17g\n", scales[si], c, ev.diam());
			}
		}
	std::fprintf(stderr, "unarias_cuadraticas=%ld\n", quad2::Eval::quad_unary_count());
}
