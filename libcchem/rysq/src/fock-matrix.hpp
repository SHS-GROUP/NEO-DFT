#ifndef _RYSQ_FOCK_MATRIX_HPP_
#define _RYSQ_FOCK_MATRIX_HPP_

#include <algorithm>
#include <cassert>

#include <rysq/fock.hpp>

#include "quartet.hpp"
#include "util.h"
#include "cxx/utility/permute.hpp"

namespace rysq {
    namespace fock {

	typedef Fock::fock_set fock_set;
	typedef Fock::density_set density_set;

	class Matrix {
	    density_set D_;
	    fock_set F_;
	    void *memory_;
	    size_t size_;
	    int ni, nj, nk, nl;
	    int shuffle;
	    bool tbraket, tij, tkl;
	public:
	    Matrix(const Quartet<Shell> &quartet) {
		this->ni = quartet[0].size;
		this->nj = quartet[1].size;
		this->nk = quartet[2].size;
		this->nl = quartet[3].size;
		this->size_ = ni*nj + nk*nl + (ni + nj)*(nk + nl) + 6;
		
		assert(posix_memalign(&this->memory_, 16, 2*size_*sizeof(double)) == 0);
		initialize(ni, nj, nk, nl, this->D_, (double*)memory_);
		initialize(ni, nj, nk, nl, this->F_, (double*)memory_ + size_);

		const bool tbra = quartet.transpose_mask & Transpose::BRA;
		const bool tket = quartet.transpose_mask & Transpose::KET;

		this->tbraket = quartet.transpose_mask & Transpose::BRAKET;
		this->tij = (tbra && !tbraket) || (tket && tbraket);
		this->tkl = (tket && !tbraket) || (tbra && tbraket);

		namespace p = cxx::utility::permutation_;
		int i = p::index(0,quartet.shuffle_mask)%2;
		int j = p::index(1,quartet.shuffle_mask)%2;
		int k = p::index(2,quartet.shuffle_mask)%2;
		int l = p::index(3,quartet.shuffle_mask)%2;

		int ik = index(i,k,tbraket);
		int il = index(i,l,tbraket);
		int jk = index(j,k,tbraket);
		int jl = index(j,l,tbraket);
		this->shuffle = SHUFFLE(ik, il, jk, jl);
	    }

	    ~Matrix() { free(this->memory_); }

	    static int index(int bra, int ket, bool transpose) {
		 return (transpose) ? bra + ket*2 : bra*2 + ket;
	    }

	    density_set& D() { return this->D_; }

	    fock_set& F() { return this->F_; }

	    void set(density_set D) {
		if (tbraket) std::swap(D[0], D[1]);
		cxx::utility::permute(D[2], D[3], D[4], D[5], this->shuffle);

 		// copy(ni, nj, const_cast<double*>(this->D_[0]), D[0], tij);
		// copy(nk, nl, const_cast<double*>(this->D_[1]), D[1], tkl);
		// copy(ni, nk, const_cast<double*>(this->D_[2]), D[2], tbraket);
		// copy(ni, nl, const_cast<double*>(this->D_[3]), D[3], tbraket);
		// copy(nj, nk, const_cast<double*>(this->D_[4]), D[4], tbraket);
		// copy(nj, nl, const_cast<double*>(this->D_[5]), D[5], tbraket);

		std::fill(F_[0], F_[0] + size_, 0.0);
	    }

	    void get(double scale, double scale2, fock_set F) {
		if (tbraket) std::swap(F[0], F[1]);
		cxx::utility::permute(F[2], F[3], F[4], F[5], this->shuffle);

   		// add(ni, nj, scale2, (this->F_[0]),  F[0], tij);
 		// add(nk, nl, scale2, (this->F_[1]),  F[1], tkl);
		// add(ni, nk, -scale, (this->F_[2]), F[2], tbraket);
 		// add(ni, nl, -scale, (this->F_[3]), F[3], tbraket);
 		// add(nj, nk, -scale, (this->F_[4]), F[4], tbraket);
		// add(nj, nl, -scale, (this->F_[5]), F[5], tbraket);
	    }

	    template<class M, typename T>
	    static void initialize(int ni, int nj, int nk, int nl, M &M6, T *Ptr) {
		util::round<2> round2;
		T *ptr = Ptr;
		M6[0] = ptr; ptr += round2(ni*nj);
		M6[1] = ptr; ptr += round2(nk*nl);
		M6[2] = ptr; ptr += round2(ni*nk);
		M6[3] = ptr; ptr += round2(ni*nl);
		M6[4] = ptr; ptr += round2(nj*nk);
		M6[5] = ptr; ptr += round2(nj*nl);
	    }

	    template <typename S, typename D>
	    void copy(int m, int n, D *dest, const S *source, bool transpose) {
		if (transpose) {
		    for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
			    dest[i+j*m] = source[j+i*n];
			}
		    }
		}
		else {
		    util::copy(m*n, dest, source);
		}
	    }


	    template <typename S, typename D>
	    void add(int m, int n, S scale, const S *source, D *dest, bool transpose) {
		if (transpose) {
		    for (int j = 0; j < n; ++j) {
			for (int i = 0; i < m; ++i) {
			    dest[j+i*n] += scale*source[i+j*m];
			}
		    }
		}
		else {
		    util::add(m*n, dest, scale, source);
		}
	    }


	};

    }
}

#endif /* _RYSQ_FOCK_MATRIX_HPP_ */
