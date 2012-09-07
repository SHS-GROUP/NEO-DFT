#ifndef HF_HF_HPP
#define HF_HF_HPP

/**
 * @file 
 * @brief  HF functions
 */

#include "basis/basis.hpp"
#include "matrix/matrix.hpp"
#include "parallel.hpp"

namespace hf {

    typedef boost::numeric::ublas::matrix<
	double, boost::numeric::ublas::column_major> Matrix;

    template<class M1, class M2>
    void Dmax(const Basis &basis, const M1 &D, M2 &Dmax);

    struct Parameters {
	Parameters(double cscale, double xscale,
		   double cutoff = 1e-10)
	    : cscale(cscale), xscale(xscale), cutoff(cutoff) {}
	double cscale, xscale;
	double cutoff;
    };

    /** 
     * Construct Fock JK matrix operator
     * 
     * @param basis basis set
     * @param D density matrix
     * @param F fock matrix
     * @param tolerance integral tolerance threshold
     * @param screen integral screening matrix
     */
    template<class M>
    void fock(const Basis &basis,
	      const matrix::meta_matrix<BlockMatrix> &D,
	      matrix::meta_matrix<BlockMatrix> &F,
	      const Parameters &parameters,
	      Parallel &pe,
	      const M *Kmax = NULL, const M *Dmax = NULL);

}

#include "hf/screen.hpp"

#endif // HF_HF_HPP
