#ifndef DFT_DFT_HPP
#define DFT_DFT_HPP

#include "dft/types.hpp"
#include "core/wavefunction.hpp"
#include "array_ref.hpp"

#include <string>

namespace dft {

    XC fock(const Wavefunction &W, const std::string &functional, 
	    const_array_ref<double,3> dr, const_array_ref<double> w,
	    matrix_type &F,
	    const Parameters &parameters = Parameters());

}

#endif // DFT_DFT_HPP
