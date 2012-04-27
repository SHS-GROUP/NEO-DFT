#ifndef RYSQ_CUDA_KERNEL_QUARTET_HPP
#define RYSQ_CUDA_KERNEL_QUARTET_HPP

#include "kernel/vector.hpp"


#include <math.h>
#include <cuda.h>
#include <cuda_runtime.h>
#include <boost/assert.hpp>
#include <boost/config.hpp>

#include "kernel/vector.hpp"
#include "cuda/detail.hpp"
#include "boost/cuda/device.hpp"
#include "boost/cuda/types.hpp"

namespace rysq {
namespace cuda {
namespace kernel {

    typedef boost::cuda::device::Block thread_block;
    using rysq::kernel::vector;

    struct State {
	vector<3> r[2], dr;
	double dr2;
	ushort K;
	// State(const ushort &K) : K(K) {}
	// static size_t size(ushort K) { return K*9; }
	// size_t size()const { return K*9; }
    protected:
	BOOST_GPU_ENABLED
	size_t data(double *data,
		    double* &e, double* &A, double* &A1,
		    double* &rA, double* &rAi) {
	    e = data;
	    A = data + K;
	    A1 = data + 2*K;
	    rA = data + 3*K;
	    rAi = data + 6*K;
	    return 9*K;
	}
	__device__
	void primitives(const double2 *p,
			double *e, double *A, double *A1,
			double *rA, double *rAi,
			const thread_block &block) {
	    // compute bra values
	    for (ushort Kij = block.rank; Kij < this->K; Kij += block.size) {
		double2 aij = p[Kij];
		double a = aij.x + aij.y;
		double a1 = 1.0/a;
		// printf("%f\n", aij.x + aij.y);
		A[Kij] = a;
		A1[Kij] = a1;
		e[Kij] = exp(-aij.x*aij.y*a1*this->dr2);
		for (ushort i = 0; i < 3; ++i) {
		    rA[i+Kij*3] = a1*(aij.x*this->r[0][i] + aij.y*this->r[1][i]);
		    rAi[i+Kij*3] = rA[i+Kij*3] - this->r[0][i];
		}
	    }
	}
    };

    struct Bra : State {
	ushort mi, mj, m;
	double *e, *rA;
	double *A, *rAi;
	double *A1;
	BOOST_GPU_ENABLED ushort mij() const { return mi*mj; }
	BOOST_GPU_ENABLED size_t data(double *data) {
	    return State::data(data, e, A, A1, rA, rAi);
	}
	__device__ void primitives(const double2 *p,
				   const thread_block &block) {
	    State::primitives(p, e, A, A1, rA, rAi, block);
	}
    };


    struct Ket : State {
	ushort nk, nl, n;
	double *e, *rB;
	double *B, *rBi;
	double *B1;
	BOOST_GPU_ENABLED ushort nkl() const { return nk*nl; }
	BOOST_GPU_ENABLED size_t data(double *data) {
	    return State::data(data, e, B, B1, rB, rBi);
	}
	__device__ void primitives(const double2 *p,
				   const thread_block &block) {
	    State::primitives(p, e, B, B1, rB, rBi, block);
	}
    };


    struct Quartet {
	Quartet(const detail::Quartet &quartet) : centers_() {
	    initialize(quartet);
	} 
	Quartet(const detail::Quartet &quartet,
		const detail::Centers centers,
		unsigned short permutation)
	    : centers_(centers)
	{
	    permutation_ = permutation;
	    initialize(quartet);
#if !RYSQ_CUDA_KERNEL_USE_CONSTANT
	    index2d = quartet.index2d().begin().get();
	    data =  quartet.data().begin().get();
#endif
	}
	BOOST_GPU_ENABLED
	size_t size() const { return size_; }
	__device__ void centers(const int *index, State &bra, State &ket,
				const thread_block &block) {
	    // compute bra and ket distances
	    for (ushort i = block.rank; i < 3; i +=  block.size) {
		bra.r[0][i] = centers_[index[0]].elems[i];
		bra.r[1][i] = centers_[index[1]].elems[i];
		ket.r[0][i] = centers_[index[2]].elems[i];
		ket.r[1][i] = centers_[index[3]].elems[i];
		bra.dr[i] = bra.r[0][i] - bra.r[1][i];
		ket.dr[i] = ket.r[0][i] - ket.r[1][i];
	    }
	    if (block.rank == 0) {
		bra.dr2 = bra.dr.dot();
		ket.dr2 = ket.dr.dot();
	    }
	    __syncthreads();
	}
	BOOST_GPU_ENABLED
	void initialize(Bra &bra, Ket &ket) {
	    bra.K = K(0)*K(1);
	    bra.mi = L(0) + 1;
	    bra.mj = L(1) + 1;
	    bra.m = bra.mi + bra.mj - 1;
	    ket.K = K(2)*K(3);
	    ket.nk = L(2) + 1;
	    ket.nl = L(3) + 1;
	    ket.n = ket.nk + ket.nl - 1;
	}
	BOOST_GPU_ENABLED
	ushort K(int i) const { return ((K_ >> i*4) & 0xf); }
	BOOST_GPU_ENABLED
	size_t K() const { return K(0)*K(1)*K(2)*K(3); }
	BOOST_GPU_ENABLED
	ushort L(int i) const { return abs(type(i)); }
	// BOOST_GPU_ENABLED ushort L(int i) const { return ((L_ >> i*4) & 0xf); }
	BOOST_GPU_ENABLED
	rysq::type type(int i) const {
	    return rysq::type(((type_ >> i*4) & 0xf) - 1);
	}
	BOOST_GPU_ENABLED
	bool hybrid(int i) const { return type(i) == rysq::SP; }
	BOOST_GPU_ENABLED
	ushort size(int i) const {
	    return ((L(i)+2)*(L(i)+1))/2 + hybrid(i);
	}
	BOOST_GPU_ENABLED
	ushort permutation()const { return permutation_; }
#if !RYSQ_CUDA_KERNEL_USE_CONSTANT
	const ushort3 *index2d;
	const void *data;
#endif
    private:
	unsigned short permutation_;
	unsigned short size_;
	unsigned short K_, type_; //L_, type_;
	const detail::Centers centers_;
	template<class Q>
	void initialize(const Q &quartet) {
	    size_ = quartet.size();
	    // permutation_ = quartet.permutation();
	    // L_ = 0;
	    type_ = 0;
	    K_ = 0;
	    for (int i = 0; i < 4; ++i) {		        
		// L_ = L_ | (quartet[i].L << i*4);
		type_ = type_|((quartet[i].type+1) << i*4);
		K_ = K_ | (quartet[i].K << i*4);
	    }
	}
    };

}
}
}

#endif /* RYSQ_CUDA_KERNEL_QUARTET_HPP */
