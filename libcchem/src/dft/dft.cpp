#include "dft/dft.hpp"
#include "dft/grid.hpp"

#include <boost/numeric/ublas/matrix.hpp>

namespace dft {

    XC fock(const Wavefunction &W, const std::string &functional, 
	    const_array_ref<double,3> dr, const_array_ref<double> w,
	    matrix_type &F, const Parameters &parameters) {

	namespace ublas = boost::numeric::ublas;
	using boost::numeric::ublas::triangular_adaptor;

	Grid grid(parameters);
	Functional f(functional);
	XC xc = grid.evaluate(W, dr, w, f);
	grid.fock(triangular_adaptor<matrix_type, ublas::upper>(F));

	return xc;
    }

}
