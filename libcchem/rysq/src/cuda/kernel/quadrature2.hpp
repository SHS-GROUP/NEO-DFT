#ifndef RYSQ_CUDA_KERNEL_QUADRATURE2_HPP
#define RYSQ_CUDA_KERNEL_QUADRATURE2_HPP

#include <algorithm>

#include <boost/utility/binary.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/preprocessor/seq/for_each_product.hpp>
#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/seq/first_n.hpp>
#include <boost/preprocessor/seq/enum.hpp>

#include "boost/cuda/exception.hpp"
#include "boost/cuda/types.hpp"
#include "boost/cuda/device/device.hpp"
#include "boost/cuda/device/algorithm.hpp"

#include "externals/cxx/utility.hpp"

#include "transpose.hpp"
#include "cuda/kernel/kernel.hpp"
#include "cuda/detail.hpp"
#include "cuda/kernel/quartet.hpp"
#include "cuda/kernel/recurrence.hpp"
#include "cuda/kernel/transfer.hpp"
#include "cuda/kernel/global.hpp"
#include "cuda/kernel/properties.hpp"

#include "kernel/quadrature2-impl.hpp"

#include "roots/roots.hpp"

namespace rysq {
namespace cuda {
namespace kernel {


struct coefficients {

    __device__
    static void load(const kernel::Quartet &quartet,
		     const double *Cps_, double* (&Cps)[4], double *tmp,
		     const kernel::thread_block &block) {
	ushort size = 0;
#pragma unroll
	for (int i = 0; i < 4; ++i) {
	    size += quartet.hybrid(i)*quartet.K(i);
	}
	copy(Cps_, tmp, size, block);
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

    template<size_t I, size_t N>
    __device__
    static void bra(const double &c, const ushort i,
		    double* const (&Csp)[4], double (&C)[N][1]) {
	C[I][0] = c;
    }

    template<size_t I, size_t N>
    __device__
    static void bra(const double &c, const ushort i,
		    double* const (&Csp)[4], double (&C)[N][2]) {
	// in this case Csp[0] and Csp[1] should be same
	C[I][0] = c*Csp[0][i];
	C[I][1] = c;
	// printf("%e %e\n", Csp[0][i], C[1][i]);
    }

    template<size_t I, size_t N>
    __device__
    static void bra(const double &c, const ushort i,
		    double* const (&Csp)[4], double (&C)[N][4]) {
	C[I][0] = c*Csp[0][i]*Csp[1][i];
	C[I][1] = c*Csp[1][i];
	C[I][2] = c*Csp[0][i];
	C[I][3] = c;
    }

    template<size_t N>
    __device__
    static void braket(double c, const ushort i, const ushort k,
    		       double* const (&Csp)[4], double (&C)[1][N]) {
	bra<0>(c, i, Csp, C);
    }

    template<size_t N>
    __device__
    static void braket(double c, const ushort i, const ushort k,
    		       double* const (&Csp)[4], double (&C)[2][N]) {
	bra<0>(c*Csp[2][k], i, Csp, C);
	bra<1>(c, i, Csp, C);
    }

    template<size_t N>
    __device__
    static void braket(double c, const ushort i, const ushort k,
		       double* const (&Csp)[4], double (&C)[4][N]) {
	bra<0>(c*Csp[2][k]*Csp[3][k], i, Csp, C);
	bra<1>(c*Csp[3][k], i, Csp, C);
	bra<2>(c*Csp[2][k], i, Csp, C);
	bra<3>(c, i, Csp, C);
    }
			       
};

template<class Braket, class Transform, typename T>
struct quadrature<Braket, Transform, T,
		  typename boost::enable_if_c<
		      (rysq::kernel::quadrature::impl<Braket>::value) &&
		      (Braket::size <= 20)>::type>
{
    typedef rysq::kernel::quadrature::impl<Braket> kernel_;
    static const bool value = true;
    static const int L = Braket::L;
    static const int N = Braket::L/2 + (Braket::L > 0);

    kernel::Quartet quartet;
    const rysq::Int4 *quartets;
    float cutoff;

    quadrature(kernel::Quartet quartet, const rysq::Int4 *quartets, double cutoff)
	: quartet(quartet),
	  quartets(quartets),
	  cutoff(cutoff/(rysq::SQRT_4PI5*quartet.K()*1e2)) {}

    __device__
    void operator()(Transform &transform) {

	namespace device = boost::cuda::device;
    
	__shared__ kernel::Bra bra;
	__shared__ kernel::Ket ket;

	const device::Block block;

	__shared__ double *Cps[4];
	__shared__ int index[4];
	{

	    if (block.rank == 0)  {
		size_t offset = 0;

		quartet.initialize(bra, ket);
		offset += bra.data(GLOBAL(shmem) + offset);
		offset += ket.data(GLOBAL(shmem) + offset);

		for (int i = 0; i < 4; ++i) {
		    Cps[i] = GLOBAL(shmem) + offset;
		    offset += ((i < 2) ? bra.K : ket.K)*quartet.hybrid(i);
		}
	    }

	    int *pindex = (int*)GLOBAL(shmem);
	    if (block.rank < 4) {
		int i = block.rank;
		index[i] = quartets[device::grid::shift_rank()].elems[i];
		pindex[i] = cxx::utility::permute(i, index,
						  quartet.permutation());
	    }
	    quartet.centers(pindex, bra, ket, block);
	}
	__syncthreads();

	const double2 *cAB = RYSQ_CUDA_KERNEL_QUARTET_DATA(quartet) + (bra.K + ket.K);
	{
	    const double *Cps_ = (const double*)(cAB + quartet.K());
	    coefficients::load(quartet, Cps_, Cps, (GLOBAL(shmem)), block);
	}
	__syncthreads();

	{
	    const double2 *cAij = RYSQ_CUDA_KERNEL_QUARTET_DATA(quartet);
	    const double2 *cAkl = cAij + bra.K;
	    bra.primitives(cAij, block);
	    ket.primitives(cAkl, block);
	}
	__syncthreads();

	double I[Braket::size];
	for (int i = 0; i < int(Braket::size); ++i) {
	    I[i] = 0;
	}

	for (ushort K = block.rank; K < quartet.K(); K += block.size) {
	    ushort Kkl = K/bra.K;
	    ushort Kij = K - Kkl*bra.K;

	    double2 AB = cAB[K];
	    double C0 = AB.x*bra.e[Kij]*ket.e[Kkl];
	    double AB1 = AB.y;

	    vector<N> t2;
	    vector<N+(L==0)> W;

	    double rho = (bra.A[Kij])*(ket.B[Kkl])*AB1;
	    double X = rho*device::distance2(&bra.rA[Kij*3], &ket.rB[Kkl*3]);
	    rysq::roots<N>(X, t2.elems, W.elems);
	    t2 *= AB1;

	    vector<3>::adapter rA(&bra.rA[Kij*3]);
	    vector<3>::adapter rB(&ket.rB[Kkl*3]);
	    vector<3> rAB(rA - rB);

	    vector<3>::adapter rAi(&bra.rAi[Kij*3]);
	    vector<3>::adapter rBk(&ket.rBi[Kkl*3]);

	    double C[Braket::ket::nc][Braket::bra::nc];
	    coefficients::braket(C0, Kij, Kkl, Cps, C);

	    kernel_::eval(bra.A[Kij], ket.B[Kkl],
	    		  rAi, rBk, rAB, bra.dr, ket.dr,
	    		  t2, W, C, I);
	}
	__syncthreads();

	device::block::reduce(I, GLOBAL(shmem));
	
	__syncthreads();	
	if (block.rank == 0)
	    kernel_::reorder(GLOBAL(shmem), GLOBAL(shmem) + Braket::size);
	__syncthreads();

	transform(device::grid::shift_rank(), index,
		  GLOBAL(shmem), Braket::size,
		  block, GLOBAL(shmem) + Braket::size);

    }
};

template<rysq::type A, rysq::type B, rysq::type C, rysq::type D, class Transform>
struct Kernel<meta::braket<A,B,C,D>, Transform,
	      typename boost::enable_if<
		  quadrature<meta::braket<A,B,C,D>, Transform> >::type>
: Eri <Transform> {

    typedef meta::braket<A,B,C,D> Braket;

    //     static size_t get_occupancy(const detail::Quartet &quartet) {
    // 	try {
    // 	    size_t threads = std::max<size_t>(multiply(block(quartet)), 32);
    // 	    return std::min(16384/(threads*60), 16384/shared(quartet));
    // 	}
    // 	catch (std::exception &e) {
    // 	    return 0;
    // 	}	
    //     }

    static dim3 block(const detail::Quartet &quartet) {
	size_t min = 4;//std::max(4, quartet.size());
	size_t max = std::min<int>(64, quartet.K());
	dim3 block(std::max(min,max));
	return block;
    }

    static size_t shared(const detail::Quartet &quartet) {
	kernel::Bra bra;
	kernel::Ket ket;
	(kernel::Quartet(quartet)).initialize(bra, ket);

	size_t primitives = bra.data(NULL) + ket.data(NULL);
	// hybrid s/p coefficients
	primitives += quartet[0].is_hybrid()*bra.K;
	primitives += quartet[1].is_hybrid()*bra.K;
	primitives += quartet[2].is_hybrid()*ket.K;
	primitives += quartet[3].is_hybrid()*ket.K;

	size_t num_threads = block(quartet).x + 1;
	size_t reduce = quartet.size();
	// memory to reorder/reduce integrals
	reduce += std::max<size_t>(quartet.size(), num_threads/2);

	size_t shared = std::max(reduce, primitives);
	//shared = 15000;
	shared *= sizeof(double);

	if (shared > (properties().sharedMemPerBlock-1024))
	    throw std::range_error("resources exceeded");

	return shared;
    }

    Kernel() : quartet_() {}

    Kernel(const detail::Quartet *quartet)
	:  block_(block(*quartet)), shared_(shared(*quartet)),
	   quartet_(quartet)
    {
	// if (quartet.K() > 1) throw std::exception();
	// if (quartet.ket().hybrid())  throw std::exception();
	// std::cout <<  "new kernel " << M << std::endl;
    }

    void operator()(const detail::Centers &centers,
		    const detail::Quartets &quartets,
		    const cuda::Eri::Parameters &p,
		    const boost::cuda::stream &stream,
		    Transform transform,
		    ushort threads, size_t shared) {

#if RYSQ_CUDA_KERNEL_USE_CONSTANT
	try {
	    copy(quartet_.data(), BOOST_PP_STRINGIZE(GLOBAL(quartet_data)));
	}
	catch (boost::cuda::exception e) {
	    std::cout << __FILE__ << ":" << __LINE__ << ": "
		      << e << std::endl;
	    throw;
	}
#endif

	ushort permutation = transpose_permutation(transpose_.value);
	kernel::Quartet quartet(*quartet_, centers, permutation);
	BOOST_VERIFY(quartet.size() == Braket::size);

	using std::max;
	using cxx::utility::qceiling;

	dim3 block = block_;
	block.x = max(ushort(block.x), threads);
	shared = max(shared_, Braket::size + shared);

	for (int i = 0; i <  quartets.size(); i += 65535) {
	    dim3 grid = dim3(std::min<size_t>(quartets.size() - i, 65535));
	    // std::cout <<  grid << block << shared << std::endl;

	    // std::cout << "gpu: " << quartets.size() << " "
	    // 	      << A << B << C << D << std::endl;
	    transform.offset = i;
	    quadrature<Braket, Transform>
		quadrature(quartet, quartets.data() + i, p.cutoff);
	    global(grid, block, shared, stream)(quadrature, transform);

	    boost::cuda::check_status();
	}
	
    }

    static Eri<Transform>* new_(const detail::Quartet &quartet,
				const rysq::Transpose &transpose) {
	return new Kernel(&quartet);
	// // std::cout <<  "new_" << std::endl;
	// try { return new Kernel(quartet); }
	// catch (const std::exception&) { return NULL; }
    }

private:
    dim3 block_;
    size_t shared_;
    const detail::Quartet *quartet_;
    Transpose transpose_;
};


} // namespace kernel
} // namespace cuda
} // namespace rysq


#endif // RYSQ_CUDA_KERNEL_QUADRATURE2_HPP
