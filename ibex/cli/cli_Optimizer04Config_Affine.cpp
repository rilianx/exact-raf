//============================================================================
// ibexopt-raak (version local de experimentos)
//============================================================================

#include "ibex_Optimizer04Config_Affine.h"
#include "ibex_LinearizerXTaylor.h"
#include "ibex_LinearizerCompo.h"
#include "ibex_LinearizerAffine2.h"
#include "ibex_LinearizerQuad2.h"
#include "ibex_Exception.h"

#include <sstream>

namespace ibex {

Optimizer04Config_Affine::Optimizer04Config_Affine(int argc, char** argv, const std::string& arith_) :
		Optimizer04Config(argc, argv), arith(arith_) {
}

// El main ya ha traducido los alias (raak, raakcompo, --raak_k) a --lr=art|compo mas --arith y
// --knomial_k, y ha fijado ibex::affine_mode() y ibex::knomial_max_k(). Aqui solo se elige el
// linearizador: AF2/RAA usan LinearizerAffine2 (= LinearizerAffineT<AF_fAF2>); AA usa
// LinearizerAffineT<AF_fAFFullI> (Affine3, parche 0007).
Linearizer* Optimizer04Config_Affine::get_linear_relax() {
	if (linearrelaxation == "xn")
		return &rec(new LinearizerXTaylor(*ext_sys));

	Linearizer* affine = NULL;
	if (linearrelaxation == "art" || linearrelaxation == "compo") {
		if (arith == "aa")
			affine = &rec(new LinearizerAffineT<AF_fAFFullI>(*ext_sys));
		else if (arith == "q2" || arith == "q2kn" || arith == "q2pc" || arith == "q2full") {
			// [raak, este repo] formas afines de grado 2 (ibex/q2/ibex_LinearizerQuad2.h)
			const quad2::FinalBound fb = arith == "q2kn" ? quad2::FB_KN : arith == "q2pc" ? quad2::FB_PC
			                           : arith == "q2full" ? quad2::FB_FULL : quad2::FB_TT;
			affine = &rec(new LinearizerQuad2(*ext_sys, fb));
		}
		else
			affine = &rec(new LinearizerAffine2(*ext_sys));
	}

	if (linearrelaxation == "art")
		return affine;

	if (linearrelaxation == "compo")
		return &rec(new LinearizerCompo(
				rec(new LinearizerXTaylor(*ext_sys)),
				*affine));

	std::stringstream ss;
	ss << "[ibexopt-raak] relajacion lineal desconocida: '" << linearrelaxation
	   << "' (opciones: xn | art | compo; alias raak | raakcompo solo en el main)";
	ibex_error(ss.str().c_str());
	return NULL;
}

} // namespace ibex
