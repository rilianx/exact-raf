// Driver de prodk.h: una forma por linea, "m h11 h12 ... hmm" (simetrica por filas).
// Salida por linea: "L U bloques tiempo_us" (tiempo de prodk_interval, mediana de REPS repeticiones).
// Argumentos: [REPS=5] [LAZY=1] [MAX_TRIPLES_M=40] [FULL_MAX=12]. Ademas imprime el ancho del encierro barato TT.
#include <algorithm>
#include <chrono>
#include <cstdio>
#include <iostream>
#include <sstream>
#include <string>
#include "prodk.h"

int main(int argc, char** argv) {
	const int REPS = argc > 1 ? std::atoi(argv[1]) : 5;
	const bool LAZY = argc > 2 ? std::atoi(argv[2]) != 0 : true;
	const int MT = argc > 3 ? std::atoi(argv[3]) : 40;
	const int FM = argc > 4 ? std::atoi(argv[4]) : 12;
	std::string line;
	while (std::getline(std::cin, line)) {
		if (line.empty()) continue;
		std::istringstream is(line);
		int m; is >> m;
		prodk::Vec H(m * m);
		for (double& x : H) is >> x;
		std::vector<double> ts;
		prodk::Result r{};
		for (int k = 0; k < REPS; k++) {
			auto t0 = std::chrono::steady_clock::now();
			r = prodk::prodk_interval(H, m, LAZY, MT, FM);
			auto t1 = std::chrono::steady_clock::now();
			ts.push_back(std::chrono::duration<double, std::micro>(t1 - t0).count());
		}
		std::sort(ts.begin(), ts.end());
		double dlo = 0, dhi = 0, off = 0;
		for (int i = 0; i < m; i++) for (int j = 0; j < m; j++) {
			const double v = H[i * m + j];
			if (i == j) { dlo += std::min(v, 0.0); dhi += std::max(v, 0.0); } else off += std::fabs(v);
		}
		std::printf("%.17g %.17g %d %.3f %.17g\n", r.L, r.U, r.blocks, ts[ts.size() / 2], (dhi + off) - (dlo - off));
	}
}
