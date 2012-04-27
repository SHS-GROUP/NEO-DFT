#ifndef DFT_TYPES_HPP
#define DFT_TYPES_HPP

#include <boost/numeric/ublas/matrix.hpp>

namespace dft {

    typedef boost::numeric::ublas::matrix<
	double, boost::numeric::ublas::column_major> matrix_type;

    struct Parameters {
	Parameters() : rcut(1e-15), ccut(1e-15) {}
	double rcut, ccut;
    };
    
    struct XC {
	double Xa, Xg, Ec, rho, aa, bb, ab;
	XC& operator+=(const XC &xc) {
	    Xa += xc.Xa;
	    Xg += xc.Xg;
	    Ec += xc.Ec;
	    rho += xc.rho;
	    aa += xc.aa;
	    bb += xc.bb;
	    ab += xc.ab;
	    return *this;
	}
    };

}

#endif // DFT_TYPES_HPP
