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
#include "cuda/kernel/device.hpp"

namespace rysq {
namespace cuda {
namespace kernel {

    typedef detail::Shell Shell;

    typedef ::cuda::device::Threads Threads;
    typedef ::cuda::device::Warp Warp;
    typedef ::cuda::device::Block Block;

    using rysq::kernel::vector;

    __device__
    inline int divmod(int ij, int K, int &i, int &j) {
	j = ij/K;
	i = ij - j*K;
	return 0;
    }

    struct State {
	vector<3> r[2], dr;
	double dr2;
	int L, K;
	int ptr_;
	BOOST_GPU_ENABLED
	int size() const { return 6*K; }
    protected:
	BOOST_GPU_ENABLED
	State() {}
	BOOST_GPU_ENABLED
	void initialize(const Shell &p, const Shell &q) {
	    this->L = p.L() + q.L();
	    this->K = p.K()*q.K();
	    this->ptr_ = 0;
	}
	BOOST_GPU_ENABLED
	State(const Shell &p, const Shell &q) {
	    initialize(p, q);
	}
	__device__
	State(const Shell &p, const Shell &q, double *data) {
	    initialize(p, q);
	    this->ptr_ = data - this->data();
	}
	__device__ double* data() const {
	    extern __shared__ double shmem[];
	    return shmem + ptr_;
	}
	__device__ double e(int i) const { return data()[0*K+i]; }
	__device__ double X(int i) const { return data()[1*K+i]; }
	//__device__ double X1(int i) const { return 1/X(i); }//data_[1*K+i]; }
	__device__ double X1(int i) const { return data()[2*K+i]; }
	__device__ double* rX(int i) const { return data()+3*K+i*3; }

    public:
	__device__
	void primitive(const Shell::Primitive &p,
		       const Shell::Primitive &q,
		       int ij) {
	    double *e = data()+0*K;
	    double *A = data()+1*K;
	    double *A1 = data()+2*K;
	    double *rA = data()+3*K;
	    double a = p.a + q.a;
	    double a1 = 1/a;
	    e[ij] = exp(-p.a*q.a*a1*this->dr2)/a;
	    e[ij] *= (p.C*q.C);
	    A[ij] = a;
	    A1[ij] = a1;
#pragma unroll
	    for (int k = 0; k < 3; ++k) {
		rA[k+ij*3] = a1*(p.a*this->r[0][k] + q.a*this->r[1][k]);
	    }
	}

	__device__
	void primitives(const Shell &P, const Shell &Q,
			double* (&C)[2], int thread, int block) {
	    int K = P.K()*Q.K();
	    int Ki = P.K();
	    for (int ij = thread; ij < K; ij += block) {
		int i, j;
		divmod(ij, Ki, i, j);
		Shell::Primitive p = P(i);
		Shell::Primitive q = Q(j);
		primitive(p, q, ij);
		if (P == rysq::SP) C[0][ij] = p.Cs;//P(i,0)/p;
		if (Q == rysq::SP) C[1][ij] = q.Cs;//Q(j,0)/q;
	    }
	}

    };

    struct Bra : State {
	BOOST_GPU_ENABLED
	Bra() {}
	template<class Q>
	BOOST_GPU_ENABLED Bra(const Q &quartet)
	    : State(quartet[0], quartet[1]) {}
	template<class Q>
	__device__ Bra(const Q &quartet, double *shmem)
	    : State(quartet[0], quartet[1], shmem) {}
	__device__ double e(int i) const { return State::e(i); }
	__device__ double A(int i) const { return State::X(i); }
	__device__ double A1(int i) const { return State::X1(i); }
	__device__ double* rA(int i) const { return State::rX(i); }
	template<class Q>
	__device__
	void primitives(const Q &quartet, double* (&C)[2],
			int thread, int block) {
	    State::primitives(quartet[0], quartet[1], C, thread, block);
	}
    };


    struct Ket : State {
	BOOST_GPU_ENABLED
	Ket() {}
	template<class Q>
	BOOST_GPU_ENABLED Ket(const Q &quartet)
	    : State(quartet[2], quartet[3]) {}
	template<class Q>
	__device__ Ket(const Q &quartet, double *shmem)
	    : State(quartet[2], quartet[3], shmem) {}
	__device__ double e(int i) const { return State::e(i); }
	__device__ double B(int i) const { return State::X(i); }
	__device__ double B1(int i) const { return State::X1(i); }
	__device__ double* rB(int i) const { return State::rX(i); }
	template<class Q>
	__device__
	void primitives(const Q &quartet, double* (&C)[2],
			int thread, int block) {
	    State::primitives(quartet[2], quartet[3], C, thread, block);
	}
    };


    struct Quartet {

    private:
	// Shell elems[4];
	const detail::Centers centers_;

	struct {
	    const Shell::Primitive* data[4];
	    signed char type[4], K[4];
	} shell_;
	int size_;
	int K_;
	char order[4];

    public:

	explicit
	Quartet(const Shell::Quartet &quartet,
		const detail::Centers &centers = detail::Centers())
	    :  centers_(centers)
	{
	    for (int i = 0; i < 4; ++i) {
		//elems[i] = quartet.elems[i];
		Shell s = quartet.elems[i];
		//elems[i] = Shell(s.data(), s.type(), 0);//quartet.elems[i];
		shell_.data[i] = s.data();
		shell_.type[i] = s.type();
		shell_.K[i] = s.K();
		order[i] = quartet.order[i];
		//std::cout << order[i] << std::endl;
	    }
	    size_ = quartet.size();
	    K_ = quartet.K();
	}

	__device__ __forceinline__
	void centers(const int *index,
		     State &bra, State &ket,
		     int thread) const {
	    // compute bra and ket distances

	    if (thread < 3) {
		int i = thread;
		// printf("%i %i\n", order[0], index[order[0]]);
		// printf("%i %i\n", order[1], index[order[1]]);
		// printf("%i %i\n", order[2], index[order[2]]);
		// printf("%i %i\n", order[3], index[order[3]]);
		// printf("%p\n", &centers_[0]);
		bra.r[0][i] = centers_[ index[order[0]] ].elems[i];
		bra.r[1][i] = centers_[ index[order[1]] ].elems[i];
		ket.r[0][i] = centers_[ index[order[2]] ].elems[i];
		ket.r[1][i] = centers_[ index[order[3]] ].elems[i];
		bra.dr[i] = bra.r[0][i] - bra.r[1][i];
		ket.dr[i] = ket.r[0][i] - ket.r[1][i];
	    }
	    if (thread == 0) {
		bra.dr2 = bra.dr.dot();
		ket.dr2 = ket.dr.dot();
	    }
	}

	BOOST_GPU_ENABLED
	Shell operator[](int i) const {
	    return Shell(shell_.data[i], rysq::type(shell_.type[i]),
			 shell_.K[i]);
	}

	BOOST_GPU_ENABLED
	int K() const { return K_; }

	BOOST_GPU_ENABLED
	int size() const { return size_; }

	template<int N>
	__device__
	void primitives(Bra *bra, double* (&C)[2],
			int thread, int block) const {
	    primitives<N>((*this)[0], (*this)[1], bra, C, thread, block);
	}

	template<int N>
	__device__
	void primitives(Ket *ket, double* (&C)[2],
			int thread, int block) const {
	    primitives<N>((*this)[2], (*this)[3], ket, C, thread, block);
	}

    private:
	template<int N, class S>
	__device__
	static void primitives(const Shell &P, const Shell &Q,
			       S *s, double* (&C)[2],
			       int thread, int block) {
	    int K = P.K()*Q.K();
	    int Ki = P.K();
	    for (int ij = thread; ij < K; ij += block) {
		int i, j;
		divmod(ij, Ki, i, j);
		Shell::Primitive p = P(i);
		Shell::Primitive q = Q(j);
#pragma unroll
		for (int i = 0; i < N; ++i) {
		    s[i].primitive(p, q, ij);
		}
		if (P == rysq::SP) C[0][ij] = p.Cs;//P(i,0)/p;
		if (Q == rysq::SP) C[1][ij] = q.Cs;//Q(j,0)/q;
	    }
	}

    };


}
}
}

#endif /* RYSQ_CUDA_KERNEL_QUARTET_HPP */
