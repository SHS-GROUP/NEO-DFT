#include <algorithm>

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/sort.hpp>
#include <boost/mpl/at.hpp>

#include "rysq-core.hpp"
#include "cuda/fock.hpp"
#include "cuda/kernel/quartet.hpp"
#include "cuda/kernel/eri.hpp"


namespace rysq {
namespace cuda {
namespace kernel {
namespace fock {

    typedef detail::Mutex Mutex;
    typedef detail::fock_index index;

    __device__ volatile int lock_;

    struct mutex {
    private:
	int *data_;
	bool locked_;
    public:
	__device__	
	mutex(int *data) : data_(data), locked_(false) {} // (int*)&lock_) {}
	    // : data_(lock_) {}
	__device__
	~mutex() { unlock(); }
	__device__
	void lock() {
	    lock(data_);
	    locked_ = true;
	}
	__device__
	void unlock() {
	    if (locked_) unlock(data_);
	    locked_ = false;
	}
    public:
	__device__
	static void lock(int *data) {
	    while (atomicCAS((int*)data, 0, 1)) {}
	}
	__device__
	static void unlock(int *data) {
	    *data = 0;
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


    template<class S = detail::matrix_set<double> >
    struct Transform {
	typedef S Matrix;
	mutable Matrix matrix_;
	mutable Mutex mutex_;
	float scale_[2];

	Transform(const S &set, const double (&scale)[2],
		   Mutex mutex)
	    : matrix_(set), mutex_(mutex)
	{
	    this->scale_[0] = scale[0];
	    this->scale_[1] = scale[1];
	}	

	template<class BRAKET>
	struct max {
	    typedef boost::mpl::vector_c<int,
		BRAKET::A::size,
		BRAKET::B::size,
		BRAKET::C::size,
		BRAKET::D::size> QUARTET;
	    typedef typename boost::mpl::sort<QUARTET>::type sorted;
	    static const int value =
 		(boost::mpl::at_c<sorted, 2>::type::value*
		 boost::mpl::at_c<sorted, 3>::type::value);				 
	};

	template<class BRAKET>
	static int block(const Context &context) {
	    return (max<BRAKET>::value+3)/4;
	}

	template<class BRAKET>
	static size_t shmem(const Context &context, dim3 block) {
	    int warpSize = context.properties().warpSize;
	    int nt = block.x*block.y*block.z;
	    BOOST_ASSERT(nt%warpSize == 0);
	    int shared = max<BRAKET>::value;
	    return shared*std::min<int>(nt/warpSize, 2);
	}

	template<class BRAKET, typename T>
	__device__
	void operator()(BRAKET, int index, const T (&quartet)[4],
			double *Q, double *shmem,
			int rank = ::cuda::device::threads::rank(),
			int block = ::cuda::device::threads::size()) const {

	    typedef boost::mpl::vector_c<int,
		BRAKET::A::size,
		BRAKET::B::size,
		BRAKET::C::size,
		BRAKET::D::size> QUARTET;

	    shmem += (rank/warpSize)*max<BRAKET>::value;

	    int warp = (1 << rank/warpSize);
	    if (block/warpSize <= 1) warp |= 0x2; 
	    rank = rank%warpSize;
	    block = warpSize;

#define APPLY(IJ, KL, SCALE)			       \
	    apply<QUARTET>(IJ, KL, rank, block, SCALE, \
			   Q, matrix_, mutex_, shmem)

	    {
		index::tuple<0,1> ij(quartet);
		index::tuple<2,3> kl(quartet);
		double scale = rysq::SQRT_4PI5*scale_[0];
		scale /= rysq::Quartet<rysq::Shell>::symmetry(quartet);
		if (warp & 0x1) APPLY(ij, kl, scale);
		if (warp & 0x2) APPLY(kl, ij, scale);
	    }

	    {
		index::tuple<0,3> ij(quartet);
		index::tuple<1,2> kl(quartet);
		double scale = rysq::SQRT_4PI5*scale_[1];
		scale /= rysq::Quartet<rysq::Shell>::symmetry(quartet);
		if (warp & 0x1) APPLY(ij, kl, scale);
		if (warp & 0x2) APPLY(kl, ij, scale);
	    }

	    {
		index::tuple<0,2> ij(quartet);
		index::tuple<1,3> kl(quartet);
		double scale = rysq::SQRT_4PI5*scale_[1];
		scale /= rysq::Quartet<rysq::Shell>::symmetry(quartet);
		if (warp & 0x1) APPLY(ij, kl, scale);
		if (warp & 0x2) APPLY(kl, ij, scale);
	    }

#undef APPLY

	}

	template<class BRAKET, class IJ, class KL>
	__device__
	static void apply(const IJ &ij, const KL &kl,
			  int t, int nt,
			  const double &scale, const double *Q,
			  Matrix &matrix, Mutex &mtx,
			  double *shmem) {

#define SYNCHRONIZE() if (nt > warpSize) __syncthreads()

	    namespace mpl = boost::mpl;

	    mutex mutex(&mtx[ij.j][ij.i]);

	    int nkl = (mpl::at_c<BRAKET, KL::first>::type::value*
		       mpl::at_c<BRAKET, KL::second>::type::value);

	    for (int i = t; i < nkl; i += nt) {
		shmem[i] = matrix.density(kl)[i];
	    }
	    SYNCHRONIZE();

	    int nij = (mpl::at_c<BRAKET, IJ::first>::type::value*
		       mpl::at_c<BRAKET, IJ::second>::type::value);

	    do {
		double f = 0;
		if (t < nij) f = apply<BRAKET>(ij, kl, t, Q, shmem);
		if (t == 0) mutex.lock();
		if (t < nt) SYNCHRONIZE();
		if (t < nij) matrix.fock(ij)[t] += scale*f;
		t += nt;
	    } while (t < nij); 

	    SYNCHRONIZE();
	    mutex.unlock();

#undef SYNCHRONIZE

	}

	template<class Seq, int It, int End>
	struct multiply {
	    static const int value =
		(boost::mpl::at_c<Seq,It>::type::value*
		 multiply<Seq, It+1, End>::value);
	};

	template<class Seq, int It>
	struct multiply<Seq, It, It> {
	    static const int value = 1;
	};


	template<class BRAKET, int I, int J, class = void>
	struct integral_index {
	    static const int NI = boost::mpl::at_c<BRAKET,I>::type::value;
	    template<typename T>
	    __device__ static T value(const T &t) {
		T j = t/NI;
		return ((t - j*NI)*multiply<BRAKET,0,I>::value +
			j*multiply<BRAKET,0,J>::value);
	    }
	};

	template<class BRAKET, int I, int J>
	struct integral_index<
	    BRAKET, I, J, typename boost::enable_if_c<(I == J-1)>::type>
	{
	    template<typename T>
	    __device__ static T value(const T &t) {
		return t*multiply<BRAKET,0,I>::value;
	    }
	};


	template<class BRAKET, int I, int J, int K, int L>
	__device__
	static double apply(const index::tuple<I,J>&,
			    const index::tuple<K,L>&,
			    int ij,
			    const double *Q, const double *D) {
	    using boost::mpl::at_c;;
	    int nc = multiply<BRAKET,0,K>::value;
	    int nd = (multiply<BRAKET,0,L>::value -
		      nc*at_c<BRAKET,K>::type::value);
	    double f = 0;
	    Q += integral_index<BRAKET,I,J>::value(ij);
	    for (int l = 0; l < at_c<BRAKET,L>::type::value; ++l) {
		for (int k = 0; k < at_c<BRAKET,K>::type::value; ++k) {
		    f += (*D)*(*Q);
		    ++D;
		    Q += nc;
		}
		Q += nd;
	    }
	    return f;
	}    	

    };
 
} // namespace fock
}
}
} // namespace rysq

