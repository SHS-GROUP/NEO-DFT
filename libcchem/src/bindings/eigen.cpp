#ifndef BINDINGS_EIGEN_HPP_
#define BINDINGS_EIGEN_HPP_

#ifndef BINDINGS_PROTOTYPE_ONLY
#define BINDINGS_PROTOTYPE_ONLY 1
#endif

#include "eigen/jacoby.hpp"
#include "eigen/sort.hpp"
#include "bindings/storage.hpp"

extern "C" {

    int Eigen_jacoby(int N, double *A, double *w, double *V, int ldV)
#if BINDINGS_PROTOTYPE_ONLY
	;
#else
    {
	// std::cout << "Eigen_jacoby" << std::endl;
	using namespace bindings;
	symmetric_fortran_matrix<double> A_(N, A);
	bindings::vector<double> w_(N, w);
	fortran_matrix<double> V_(ldV,N,V);
	try {
	    eigen::jacoby::apply(A_, w_, V_);
	    eigen::sort(w_, V_);
	    return 0;
	}
	catch (std::exception) { return 1; }
    }
#endif

}

#endif /* BINDINGS_EIGEN_HPP */
