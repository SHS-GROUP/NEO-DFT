#include <algorithm>

#include <boost/preprocessor/tuple/rem.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include <boost/utility/enable_if.hpp>
#include "boost/mpl/bool.hpp"

#include "boost/cuda/device/device.hpp"
#include "boost/cuda/assert.h"

#include <rysq/core.hpp>
#include "cuda/detail.hpp"
#include "cuda/kernel/quartet.hpp"
#include "cuda/kernel/eri.hpp"


namespace rysq {
namespace cuda {
namespace kernel {
namespace fock {

    typedef rysq::cuda::detail::fock_index index;

    template<size_t A0, size_t A1, size_t A2, size_t A3>
    struct mask {
	template<size_t B0, size_t B1, size_t B2, size_t B3>
	struct compare {
	    static const bool value =
		(A0 == B0 && A1 == B1 && A2 == B2 && A3 == B3);
	};
    };


    __device__ volatile int lock;//[6] = { 0, 0, 0, 0, 0, 0 };

    struct mutex {
    	__device__
    	static void lock(int &mutex) {
	    while (atomicCAS((int*)&mutex, 0, 1)) {}
    	}
    	__device__
	static void unlock(int &mutex) {
	    mutex = 0;
	}
    };


    template<int begin, int end>
    struct range {
	template<typename T, int N>
	__host__ __device__
	static T multiply(const T (&v)[N]) {
	    T value = 1;
	    for (int i = begin; i < end; ++i) {
		value *= v[i];
	    }
	    return value;
	}
    };

    template<size_t _1,size_t _2, bool consecutive = (_1 == _2 - 1)>
    struct integral_index_ {
	template<typename T, typename U>
	__device__
	static T eval(const T (&N)[4], const U &index) {
	    // int shift;
	    // if (N[_1] == 1) shift = (1<<17)/1;
	    // else if (N[_1] == 3) shift = (1<<17)/3;
	    // else if (N[_1] == 6) shift = (1<<17)/6;
	    // // else (N[_1] == 10) shift = (1<<17)/10;
	    // else shift = (1<<17)/10;
	    // T j = (index*shift)>>17;
	    // j += ((j + 1)*N[_1]-1 < index);

	    T j = index/N[_1];
	    return ((index - j*N[_1])*range<0,_1>::multiply(N) +
		    j*range<0,_2>::multiply(N));
	}
    };

    template<size_t _1,size_t _2>
    struct integral_index_<_1, _2, true> {
	template<typename T, typename U>
	__device__
	static T eval(const T (&N)[4], const U &index) {
	    return index*range<0,_1>::multiply(N);
	}
    };

    template<size_t _1,size_t _2, typename T, typename U>
    __device__
    T integral_index(const T (&N)[4], const U &index) {
	return integral_index_<_1,_2>::eval(N, index);
    }

    template<class B, class K, class enable = void>
    struct transform;

    struct element {
	__host__ __device__
	element() : value(0), flag(0) {}
	double value;
	bool flag;
    };

    template<int a, int b, int c, int d, class enable>
    struct transform<index::tuple<a,b>, index::tuple<c,d>, enable> {
	template< typename T, typename U>
	__device__
	static element apply(const T (&N)[4], const double *D, const double *Q,
			     const U &rank, double *shmem = NULL) {

	    __syncthreads();
	    if (rank < N[c]*N[d]) shmem[rank] = D[rank];
	    __syncthreads();

	    element f;
	    const T nc = range<0,c>::multiply(N);
	    const T nd = range<0,d>::multiply(N) - N[c]*nc;
	    if (rank < N[a]*N[b]) {
		f.flag = true;
	    	Q += integral_index<a,b>(N, rank);
	    	for (T j = 0; j < N[d]; ++j) {
	    	    for (T i = 0; i < N[c]; ++i) {
	    		f.value += (*shmem)*(*Q);
	    		++shmem;
	    		Q += nc;
	    	    }
	    	    Q += nd;
	    	}
	    }
	    return f;
	}
    };


    template<>
    struct transform<index::tuple<0,3>, index::tuple<1,2> > {
	template<typename T, typename U>
	__device__
	static element apply(const T (&N)[4], const double *D, const double *Q,
			     const U &rank, double *shmem = NULL) {

	    __syncthreads();
	    if (rank < N[1]*N[2]) shmem[rank] = D[rank];
	    __syncthreads();
	    
	    element f;
	    if (rank < N[0]*N[3]) {
		f.flag = true;
		Q += integral_index<0,3>(N, rank);
		for (T i = 0; i < N[1]*N[2]; ++i) {
		    f.value += shmem[i]*(*Q);
		    Q += N[0];
		}
	    }
	    return f;
	}
    };

    template<int a, int b, int c, int d>
    struct transform<index::tuple<a,b>, index::tuple<c,d>, 
		     typename boost::enable_if_c<
			 (mask<a,b,c,d>::template compare<0,1,2,3>::value ||
			  mask<a,b,c,d>::template compare<2,3,0,1>::value)
			 >::type> {
	template<typename T, typename U>
	__device__
	static element apply(const T (&N)[4], const double *D, const double *Q,
			     const U &rank, double *shmem = NULL) {
	    const bool transposed = (mask<a,b,c,d>::template compare<2,3,0,1>::value);

	    __syncthreads();
	    if (rank < N[c]*N[d]) shmem[rank] = D[rank];
	    __syncthreads();

	    element f;
	    if (rank < N[a]*N[b]) {
		f.flag = true;
	    	Q += rank*((transposed) ? N[c]*N[d] : 1);
	    	for (T ab = 0; ab < N[c]*N[d]; ++ab) {
	    	    f.value += shmem[ab]*(*Q);
	    	    Q += ((transposed) ? 1 : N[a]* N[b]);
	    	}
	    }
	    return f;
	}
    };

    template<class S = detail::matrix_set<double> >
    struct Transform {
	typedef S Matrix;
	mutable Matrix matrix;
	float scale_[2];
	int* mutex_[6];
	int offset;

	Transform(const S &set, const double (&scale)[2],
		  int* const (&mutex)[6])
	    : matrix(set)
	{
	    this->scale_[0] = scale[0];
	    this->scale_[1] = scale[1];
	    for (int i = 0; i < 6; ++i) {
		this->mutex_[i] = mutex[i];
	    }
	}

	template<int M, typename T, typename R>
	__device__
	void operator()(int index, const T (&quartet)[4],
			const R (&Q)[M], size_t size,
			const thread_block &block, double *shmem) const {
	    eri::Transform::apply<M>(Q, size, shmem, block);
	    apply(index, quartet, shmem, matrix, scale_,
		  mutex_, block, shmem + size);
	}

	template<typename T>
	__device__
	void operator()(int index, const T (&quartet)[4],
			double *Q, size_t size,
			const thread_block &block, double *shmem) const {
	    eri::Transform::apply(Q, size, Q, block);
	    apply(index, quartet, Q, matrix, scale_, mutex_, block, shmem);
	}

	template<typename T>
	__device__
	static void apply(int index, const T (&quartet)[4], const double *Q,
			  Matrix &matrix, const float (&scale)[2],
			  int* const (&mutex)[6], const thread_block &block,
			  double *shmem) {

	    double sym = 1.0/rysq::Quartet<rysq::Shell>::symmetry(quartet);
	    int rank = block.rank;

	    __syncthreads();
#define UPDATE(i) {							\
		typedef typename index::ij<i>::type ij;			\
		typedef typename index::kl<i>::type kl;			\
		element f =						\
		    transform<ij,kl>::apply(matrix.block,		\
					    matrix.density(kl(quartet)), \
					    Q, rank, shmem);		\
		int index = matrix.index(ij(quartet));			\
		if (rank == 0) mutex::lock(mutex[i][index]);		\
		__syncthreads();					\
		if (f.flag)						\
		    matrix.fock(ij(quartet))[rank] +=			\
			(sym*scale[i > 1])*f.value;			\
		__syncthreads();					\
		if (rank == 0) mutex::unlock(mutex[i][index]);		\
	    }

	    UPDATE(0);
	    UPDATE(1);
	    UPDATE(2);
	    UPDATE(3);
	    UPDATE(4);
	    UPDATE(5);
	}
    };
 
} // namespace fock
}
}
} // namespace rysq

