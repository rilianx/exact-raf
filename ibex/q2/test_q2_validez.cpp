// Prueba de validez de LinearizerQuad2: para puntos x* (centros de cajas solucion) y cajas aleatorias que los
// contienen, comprueba en cada restriccion que
//   (1) el intervalo por nodos 'ev' de Quad2 contiene g(x*) (evaluada por intervalos en el punto), y
//   (2) g(x*) cae en [c - center + row.x* + L - margen, c - center + row.x* + U + margen] (la relajacion).
// Uso: test_q2_validez <sistema.bch> <soluciones.txt (formato --boxes de ibexsolve-raak)> [cajas por punto=200] [fb=0..3] [ext=0|1]
// Con ext=1 se prueba el sistema extendido de ibexopt (ExtendedSystem(sys, 1e-8): variable del objetivo y = f(x*)
// e igualdades relajadas), que es el que recibe el linearizador en el optimizador.
#include "ibex.h"
#include "ibex_LinearizerQuad2.h"
#include <cstdio>
#include <fstream>
#include <random>
#include <sstream>
using namespace ibex;

int main(int argc, char** argv) {
	System sys0(argv[1]);
	const bool ext = argc > 5 && std::atoi(argv[5]) != 0;
	ExtendedSystem* es = ext ? new ExtendedSystem(sys0, 1e-8) : NULL;
	const System& sys = ext ? (const System&) *es : sys0;
	const int nb = argc > 3 ? std::atoi(argv[3]) : 200;
	const quad2::FinalBound fb = (quad2::FinalBound) (argc > 4 ? std::atoi(argv[4]) : 0);
	const int n = sys.nb_var;
	std::vector<Vector> pts;
	std::ifstream in(argv[2]); std::string line;
	while (std::getline(in, line)) {
		std::istringstream is(line); std::string t; is >> t;
		const int n0 = sys0.nb_var;
		Vector p0(n0); for (int i = 0; i < n0; i++) { double lo, hi; is >> lo >> hi; p0[i] = 0.5 * (lo + hi); }
		Vector p(n);
		if (ext) {   // variables del extendido: las originales con y en goal_var()
			const int g = es->goal_var(); int k = 0;
			for (int i = 0; i < n; i++) p[i] = (i == g) ? sys0.goal->eval(IntervalVector(p0)).mid() : p0[k++];
		} else p = p0;
		pts.push_back(p);
	}
	std::vector<quad2::Tape> tapes;
	for (int c = 0; c < sys.nb_ctr; c++) tapes.push_back(quad2::Tape::build(sys.ctrs[c].f));
	std::mt19937 rng(1);
	std::uniform_real_distribution<double> U01(0, 1);
	long checks = 0, bad_ev = 0, bad_rel = 0;
	for (size_t p = 0; p < pts.size(); p++)
		for (int b = 0; b < nb; b++) {
			// caja aleatoria que contiene x*: dentro de la caja inicial, con radio relativo aleatorio
			IntervalVector box(n);
			const double scale = std::pow(10.0, -6 * U01(rng));
			for (int i = 0; i < n; i++) {
				const double w = sys.box[i].diam() * scale;
				double lo = pts[p][i] - w * U01(rng), hi = pts[p][i] + w * U01(rng);
				if (U01(rng) < 0.1) lo = hi = pts[p][i];        // a veces, variable fija (radio 0)
				box[i] = Interval(lo, hi) & sys.box[i];
				if (box[i].is_empty()) box[i] = Interval(pts[p][i]);
			}
			for (int c = 0; c < sys.nb_ctr; c++) {
				if (!tapes[c].supported) continue;
				quad2::QF F; quad2::Eval E(n); Interval evq; double mag;
				if (!E.eval(tapes[c], box, F, evq, mag)) continue;
				const Interval gx = sys.ctrs[c].f.eval(IntervalVector(pts[p]));
				checks++;
				if (!evq.intersects(gx)) {
					if (bad_ev++ < 5) std::printf("ev NO contiene g(x*): ctr %d punto %zu evq=[%.17g, %.17g] g=[%.17g, %.17g]\n",
					                              c, p, evq.lb(), evq.ub(), gx.lb(), gx.ub());
				}
				std::vector<double> a(n, 0.0);
				for (size_t k = 0; k < F.lin.size(); k++) if (F.lin[k].first < n) a[F.lin[k].first] = F.lin[k].second;
				double lin = 0, scale_f = std::fabs(F.c);
				for (int i = 0; i < n; i++) if (box[i].rad() > 0) { lin += a[i] / box[i].rad() * (pts[p][i] - box[i].mid()); scale_f += std::fabs(a[i]); }
				double L, Uu; quad2::remainder_bound(F, n, fb, L, Uu);
				const double margin = quad2::QUAD2_MARGIN * (scale_f + mag + std::fabs(L) + std::fabs(Uu));
				const double lo = F.c + lin + L - margin, hi = F.c + lin + Uu + margin;
				if (gx.ub() < lo || gx.lb() > hi) {
					if (bad_rel++ < 5) std::printf("relajacion NO contiene g(x*): ctr %d punto %zu [%.17g, %.17g] g=[%.17g, %.17g] c=%.6g L=%.6g U=%.6g\n",
					                               c, p, lo, hi, gx.lb(), gx.ub(), F.c, L, Uu);
				}
			}
		}
	std::printf("comprobaciones=%ld ev_malos=%ld relajacion_mala=%ld\n", checks, bad_ev, bad_rel);
}
