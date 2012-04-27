#ifndef RYSQ_CUDA_KERNEL_QUADRATURE_HPP
#define RYSQ_CUDA_KERNEL_QUADRATURE_HPP

#include <exception>
#include <boost/utility/binary.hpp>
#include <boost/utility/binary.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/assert.hpp>

#include "boost/cuda/exception.hpp"
#include "boost/cuda/types.hpp"
#include "boost/cuda/device/device.hpp"

#include "externals/cxx/utility.hpp"

#include "transpose.hpp"
#include "cuda/detail.hpp"
#include "cuda/kernel/quartet.hpp"
#include "cuda/kernel/kernel.hpp"
#include "cuda/kernel/recurrence.hpp"
#include "cuda/kernel/transfer.hpp"
#include "cuda/kernel/global.hpp"
#include "cuda/kernel/properties.hpp"

#include "types.hpp"
#include "roots/roots.hpp"

#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/for_each.hpp>

#include "boost/utility/profiler.hpp"

// using namespace ::cuda;
// using namespace rysq::cuda;


#define Q0000 BOOST_BINARY(0000)
#define Q1000 BOOST_BINARY(0001)
#define Q1010 BOOST_BINARY(0101)
#define Q1100 BOOST_BINARY(0011)
#define Q1110 BOOST_BINARY(0111)
#define Q1011 BOOST_BINARY(1101)
#define Q1111 BOOST_BINARY(1111)

#define Q0101 BOOST_BINARY(1010)
#define Q0011 BOOST_BINARY(1100)
#define Q0001 BOOST_BINARY(1000)
#define Q0100 BOOST_BINARY(0010)


namespace rysq {
namespace cuda {
namespace kernel {


template<int mask, int N, int M = 0>
struct Parameters;

template<int mask, int N, int M, class Transform, typename T>
struct quadrature<Parameters<mask, N, M>, Transform, T> {

    kernel::Quartet quartet;
    const rysq::Int4 *quartets;
    float cutoff;

    quadrature(kernel::Quartet quartet, const rysq::Int4 *quartets, float cutoff)
	: quartet(quartet),
	  quartets(quartets),
	  cutoff(cutoff/(rysq::SQRT_4PI5*quartet.K()*1e2)) {}
    
    __device__ void operator()(Transform &transform) {

	namespace device = boost::cuda::device;

    __shared__ kernel::Bra bra;
    __shared__ kernel::Ket ket;

    const device::Block block;
    __shared__ int index[4];

    __shared__ T *Ix, *Iy, *Iz, *tmp;
    __shared__ double *C, *Cps[4];

    __shared__ double *t2W;

    // __shared__  int index[4];

    if (block.rank == 0)  {
	quartet.initialize(bra,	ket);

	size_t size2d = bra.mij()*ket.nkl();

	size_t offset = 0;//(sizeof(int)*4)/sizeof(double);

	offset += bra.data(GLOBAL(shmem) + offset);
	offset += ket.data(GLOBAL(shmem) + offset);

	C = GLOBAL(shmem) + offset;
	offset += quartet.K();

	for (int i = 0; i < 4; ++i) {
	    Cps[i] = GLOBAL(shmem) + offset;
	    offset += ((i < 2) ? bra.K : ket.K)*quartet.hybrid(i);
	}

	t2W = GLOBAL(shmem) + offset;
	offset += 2*N*quartet.K();
	T *I_ = (T*)(GLOBAL(shmem) + offset);
	Ix = I_ + 0*size2d;
	Iy = I_ + 1*size2d;
	Iz = I_ + 2*size2d;
	tmp = I_ + size2d*3*N;

    }

    {
	int *pindex = (int*)GLOBAL(shmem);
	if (block.rank < 4) {
	    int i = block.rank;
	    index[i] = quartets[device::grid::shift_rank()].elems[i];
	    pindex[i] = 
		cxx::utility::permute(i, index, quartet.permutation());
	}
	quartet.centers(pindex, bra, ket, block);
    }

    __syncthreads();

    {
	const double2 *cAij = RYSQ_CUDA_KERNEL_QUARTET_DATA(quartet);
    	const double2 *cAkl = cAij + bra.K;
    	bra.primitives(cAij, block);
    	ket.primitives(cAkl, block);
    }

    const double2 *cAB = RYSQ_CUDA_KERNEL_QUARTET_DATA(quartet) + (bra.K + ket.K);
    {
	ushort size = 0;
#pragma unroll
	for (int i = 0; i < 4; ++i) {
	    size += quartet.hybrid(i)*quartet.K(i);
	}
	double *tmp = t2W; // use as temporary
	copy((const double*)(cAB + quartet.K()), tmp, size, block);
	__syncthreads();
	
#define COEFFICIENTS(index)						\
	if (quartet.hybrid(index)) {					\
	    ushort K1 = quartet.K(index + !(index%2));			\
	    for (ushort j =  block.rank; j < K1;  j += block.size) {	\
		char K0 = quartet.K(index - (index%2));			\
		for (char i = 0; i < K0;  ++i){				\
		    Cps[index][i+j*K0] = tmp[((index%2 == 0) ? i : j)];	\
		}							\
	    }								\
	    tmp += quartet.K(index);					\
	}
	
	COEFFICIENTS(0);
	COEFFICIENTS(1);
	COEFFICIENTS(2);
	COEFFICIENTS(3);
#undef COEFFICIENTS
    }
    __syncthreads();

    // known at compile time, correct up to eight roots
    const ushort N2 = (N < 5) ? N + N%2 : 8;

    // partition block into N2 parts
    for (ushort K = block.rank/N2; K < bra.K*ket.K; K += block.size/N2) {
    	ushort Kkl = K/bra.K;
    	ushort Kij = K - Kkl*bra.K;

    	ushort pRank = block.rank%N2;
    	double2 AB = cAB[K];
    	if (pRank == 0) C[K] = AB.x*bra.e[Kij]*ket.e[Kkl];

    	double *t2 = t2W + 2*N*K;
    	double *W = t2 + N;

    	// broken: if(fabs(C[K]) < cutoff) continue; 
    	double AB1 = AB.y;

    	double rho = (bra.A[Kij])*(ket.B[Kkl])*AB1;
    	double X = rho*device::distance2(&bra.rA[Kij*3], &ket.rB[Kkl*3]);
    	rysq::roots<N>(X, t2, W, pRank);
    	if (pRank < N) t2[pRank] *= AB1;
    }

    ushort3 index2d[M];
    T Q[M];
    unsigned int hybrid = 0;
#pragma unroll
    for (int i = 0; i < M; ++i) {
    	ushort j = block.rank + i*block.size;
    	index2d[i] = (j < quartet.size()) ?
	    RYSQ_CUDA_KERNEL_QUARTET_INDEX2D(quartet)[j] : make_ushort3(0,0,0);
    	Q[i] = 0.0;
	ushort size = 1;
#pragma unroll
	for (int k = 0; k < 4; ++k) {
	    if (quartet.hybrid(k)) hybrid |= ((j/size)%4 == 0) << (i*4+k);
	    size *=quartet.size(k);
	}
    }

    __syncthreads();

    size_t size2d = bra.mij()*ket.nkl();
    ushort stride = size2d*3;

    // ket contractions
    for (ushort Kkl = 0; Kkl < ket.K; ++Kkl) {
    	for (ushort Kij = 0; Kij <  bra.K; ++Kij) {
    	    ushort K = Kij + Kkl*bra.K;

	    double Cps_max = 1.0;
#pragma unroll
	    for (int k = 0; k < 4; ++k) {
		if (quartet.hybrid(k))
		    Cps_max = max(Cps_max, fabs(Cps[k][((k < 2)? Kij : Kkl)]));
	    }
	    if (fabs(Cps_max*C[K]) < cutoff) continue;

	    double *t2 = t2W + 2*N*K;
	    double *W = t2 + N;

	    __syncthreads();
	    __shared__ double rAB[3];

	    if (block.rank < 3)
	    	rAB[block.rank] = bra.rA[block.rank+Kij*3] - ket.rB[block.rank+Kkl*3];
	    __syncthreads();

	    T *G = (!(mask & Q0101) or (mask == Q1111)) ? Ix : tmp;
	    const ushort n = (mask & Q0011) ? ket.n : 1;

	    recurrence<N>(bra.m, n,
	    		  bra.A[Kij], ket.B[Kkl], bra.A1[Kij], ket.B1[Kkl],
	    		  rAB, &bra.rAi[Kij*3], &ket.rBi[Kkl*3], t2, W, G,
	    		  block.size, block.rank);
	    __syncthreads();


	    if (!(mask & Q0011)) { // xx|ss
	    	transfer1<N>(bra.mi, bra.mj, 1, bra.dr, G, Ix);
	    }
	    else if (!(mask & Q0001)) { // xx|xs
	    	transfer1<N>(bra.mi, bra.mj, ket.nk, bra.dr, G, Ix);
	    }
	    else if (!(mask & Q0100)) { // xs|xx
	    	transfer2<N>(bra.mi, ket.nk, ket.nl, ket.dr, G, Ix);
	    }
	    else if (mask == Q1111) { // xx|xx
	    	transfer1<N>(bra.mi, bra.mj, ket.n, bra.dr, G, (T*)tmp);
	    	__syncthreads();
	    	transfer2<N>(bra.mi*bra.mj, ket.nk, ket.nl, ket.dr, (T*)tmp, Ix);
	    }
	    __syncthreads();

#pragma unroll
	    for (int i = 0; i < M; ++i) {
		T c = C[K];
		if ((hybrid >> i*4) & 0x1) c *= Cps[0][Kij];
		if ((hybrid >> i*4) & 0x2) c *= Cps[1][Kij];
		if ((hybrid >> i*4) & 0x4) c *= Cps[2][Kkl];
		if ((hybrid >> i*4) & 0x8) c *= Cps[3][Kkl];

		const ushort3 &idx = index2d[i];
		T q = 0.0;
#pragma unroll
		for (ushort a = 0, off = 0; a < N; ++a, off += stride) {
		    q += Ix[idx.x+off]*Iy[idx.y+off]*Iz[idx.z+off];  
		}
		Q[i] += c*q;
	    }

	}
    }    

    // double Q_[M];
    // for (int i = 0; i < M; ++i) { Q_[i] = Q[i]; }
    __syncthreads();
    // int *index_ = (int*)(GLOBAL(shmem));
    // double *shmem = (double*)(index_ + 4);
    // int index[4] = { index_[0], index_[1], index_[2], index_[3] };
    transform.operator()<M>(device::grid::shift_rank(), index,
			    Q, quartet.size(), block, GLOBAL(shmem));
}
};


template<int mask, int N, int M_, class Transform>
struct Kernel<Parameters<mask, N, M_>, Transform> : Eri<Transform> {

    static const int M = M_;
    typedef quadrature<Parameters<mask, N, M>, Transform> Quadrature;
    typedef global::function<Quadrature, Transform> Function;

    static typename Function::attributes_type attributes() {
	return Function::attributes();
    }

    static size_t get_occupancy(const detail::Quartet &quartet) {
      try {
	    size_t threads = std::max<size_t>(multiply(block(quartet)), 32);
	    size_t max_registers = properties().regsPerBlock;
	    size_t max_shared = properties().sharedMemPerBlock;
	    return std::min(max_registers/(threads*60),
			    max_shared/shared(quartet));
	}
      catch (std::range_error) {
	     return 0;
	 }	
    }

    static dim3 block(const detail::Quartet &quartet) {
	using namespace cxx::utility;
	size_t m = (quartet[0].L + quartet[1].L + 1);
	size_t n = (quartet[2].L + quartet[3].L + 1);
	dim3 block;
	block.x = ceiling2(max<int>(m, n, N));
	block.y = 3;
	block.z = max(N, qceiling<int>(quartet.size(), (M*block.x*block.y)));

	if (multiply(block) > attributes().maxThreadsPerBlock) {
	    // std::cout << quartet << " too many registers: "
	    // 	      << multiply(block)*60 << " > "
	    // 	      << properties().regsPerBlock << std::endl;
	    throw std::range_error("too many threads");
	}
	return block;
    }

    static size_t shared(const detail::Quartet &quartet) {
	kernel::Bra bra;
	kernel::Ket ket;
	(kernel::Quartet(quartet)).initialize(bra, ket);

	size_t	shared = bra.data(NULL) + ket.data(NULL);
	shared += quartet.K(); // coefficients

	// hybrid s/p coefficients
	shared += quartet[0].is_hybrid()*bra.K;
	shared += quartet[1].is_hybrid()*bra.K;
	shared += quartet[2].is_hybrid()*ket.K;
	shared += quartet[3].is_hybrid()*ket.K;
	
	shared += 2*N*quartet.K(); // roots and weights
	shared += bra.mij()*ket.nkl()*(3*N); // 2-D integrals

	// temporary transfer memory
	size_t tmp = 0;
	if (mask == Q1000 || mask == Q1010) tmp = 0;
	else if (mask == Q1111) tmp = (bra.mi*bra.mj)*(ket.n)*(3*N);
	else tmp = (bra.m)*(ket.n)*(3*N);
	shared += tmp;

	shared *= sizeof(double);

	if (shared > (properties().sharedMemPerBlock - attributes().sharedSizeBytes)) {
	    // std::cout << quartet << " too much memory: "
	    // 	      << shared << " < "
	    // 	      << properties().sharedMemPerBlock-1024 << std::endl;
	    throw std::range_error("resources exceeded");
	}

	return shared;
    }

    Kernel() : quartet_() {}

    Kernel(const detail::Quartet *quartet)
	:  block_(block(*quartet)), shared_(shared(*quartet)),
	   quartet_(quartet)
    {
	BOOST_VERIFY(quartet->L()/2 +1 == N);
    }

    void operator()(const detail::Centers &centers,
		    const detail::Quartets &quartets,
		    const cuda::Eri::Parameters &p,
		    const boost::cuda::stream &stream,
		    Transform transform,
		    ushort threads, size_t shared) {

	BOOST_VERIFY(quartet_->index2d().begin());

#if RYSQ_CUDA_KERNEL_USE_CONSTANT
	try {
	    //BOOST_PROFILE_TYPE(Eri<Transform>,  "memcpy");
	    copy(quartet->data(), BOOST_PP_STRINGIZE(GLOBAL(quartet_data)), stream);
	    boost::cuda::check_status();
	    copy(quartet_->index2d(), BOOST_PP_STRINGIZE(GLOBAL(quartet_index2d)), stream);
	    boost::cuda::check_status();
	}
	catch (boost::cuda::exception e) {
	    std::cout << __FILE__ << ":" << __LINE__ << ": "
		      << e << std::endl;
	    throw;
	}
#endif

	// BOOST_PROFILE_TYPE(Eri<Transform>);

	ushort permutation = transpose_permutation(transpose_.value);
	kernel::Quartet quartet(*quartet_, centers, permutation);

	using std::max;
	using cxx::utility::qceiling;

	dim3  block = block_;
	block.z = max(qceiling<typeof(block.z)>(threads, block.x*block.y), block.z);
	if (block.x* block.y* block.z < threads)
	    throw std::runtime_error(__FILE__);

	shared = max(this->shared_, shared);
	
	for (int i = 0; i <  quartets.size(); i += 65535) {
	    dim3 grid = dim3(std::min<size_t>(quartets.size() - i, 65535));
//	    std::cout <<  grid << block << shared << std::endl;
	    transform.offset = i;
	    Quadrature kernel(quartet, quartets.data() + i, p.cutoff);
	    global(grid, block, shared, stream)(kernel, transform);

	    boost::cuda::check_status();
	}

    }

private:
    dim3 block_;
    size_t shared_;
    const detail::Quartet *quartet_;
    Transpose transpose_;
};




struct configuration {

    struct data {
	int M;
	size_t occupancy;
	const detail::Quartet &quartet;
	data(const detail::Quartet &quartet)
	    : M(0), occupancy(0), quartet(quartet) {}
    };

    struct find {
	data &data_;
	explicit find(configuration::data &data) : data_(data) {}
	template<class K>
	void operator()(const K &kernel) {
	    int o = kernel.get_occupancy(data_.quartet);
	    if (o > data_.occupancy) {
		// std::cout << K::M << std::endl;
		data_.occupancy = o;
		data_.M = K::M;
	    }
	}
    };

    struct create {
	typedef void* pointer;
	const configuration::data &data;
	pointer &kernel;
	explicit create(const configuration::data &data, pointer &kernel)
	    : data(data), kernel(kernel) {}
	template<class K>
	void operator()(const K &kernel) {
	    if (data.M == K::M) {
		if (this->kernel)
		    throw std::runtime_error("runtime_error");
		this->kernel = new K(&data.quartet);
	    }
	}
	template<class T>
	operator Eri<T>*() {
	    return static_cast<Eri<T>*>(this->kernel);
	}
    };

};

namespace mpl = boost::mpl;

template<int mask, int N, class Transform>
struct Kernel<Parameters<mask, N>, Transform> {
    typedef mpl::vector_c<int,1,3,6> M;

    template<class M>
    struct kernel_type {
	typedef Kernel<Parameters<mask, N, M::value>, Transform> type;
    };

    typedef typename mpl::transform<M, kernel_type<mpl::_1>
    				    >::type kernel_types;

    static Eri<Transform>* new_(const detail::Quartet &quartet,
				const rysq::Transpose &transpose) {
	//std::cout <<  "new_" << std::endl;
	configuration::data data(quartet);
	mpl::for_each<kernel_types>(configuration::find(data));
	void *kernel = NULL;
	mpl::for_each<kernel_types>(configuration::create(data, kernel));
	return static_cast<Eri<Transform>*>(kernel);
    }
};


} // namespace kernel
} // namespace cuda
} // namespace rysq


#endif // RYSQ_CUDA_KERNEL_QUADRATURE_HPP
