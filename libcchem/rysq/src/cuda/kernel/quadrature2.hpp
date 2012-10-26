#ifndef RYSQ_CUDA_KERNEL_QUADRATURE2_HPP
#define RYSQ_CUDA_KERNEL_QUADRATURE2_HPP

#include <algorithm>

#include <boost/typeof/typeof.hpp>
#include <boost/utility/enable_if.hpp>
#include "boost/utility/profiler.hpp"

#include "cuda/kernel/kernel.hpp"
#include "cuda/kernel/device.hpp"
#include "cuda/kernel/quartet.hpp"
#include "cuda/kernel/global.hpp"

#include "kernel/quadrature2-impl.hpp"
#include "roots/roots.hpp"

namespace rysq {
namespace cuda {
namespace kernel {


template<class Braket, class Transform, typename T = double, typename  = void>
struct quadrature2 : boost::mpl::false_ {};
 
template<class Braket, class Transform, typename T>
struct quadrature2<Braket, Transform, T,
		  typename boost::enable_if_c<
		      (rysq::kernel::quadrature::impl<Braket>::value) &&
		      //false &&
		      (Braket::size <= 24) && //false &&
		      (Braket::A::L >= Braket::B::L &&
		       Braket::C::L >= Braket::D::L)
			  >::type>
{
    typedef rysq::kernel::quadrature::impl<Braket> kernel_;
    static const bool value = true;
    static const int L = Braket::L;
    static const int N = Braket::L/2 + 1;
    static const int M = (N == 3) ? 4 : N;

    static const int SP_ = 
	((int(Braket::A::type == rysq::SP) << 0) |
	 (int(Braket::B::type == rysq::SP) << 1) |
	 (int(Braket::C::type == rysq::SP) << 2) |
	 (int(Braket::D::type == rysq::SP) << 3));

    kernel::Quartet quartet;
    detail::Quartets quartets;
    double cutoff;

    quadrature2(kernel::Quartet quartet, detail::Quartets quartets, double cutoff)
	: quartet(quartet),
	  quartets(quartets),
	  cutoff(cutoff/(quartet.K()*rysq::SQRT_4PI5*1e2)) {}

    __device__
    void operator()(Transform &transform) {

    namespace device = ::cuda::device;

    extern __shared__ double shmem_[];
    __shared__ int ptr[2];
    
    __shared__ kernel::Bra bra_;
    __shared__ kernel::Ket ket_;

    __shared__ kernel::Bra bras[2];
    __shared__ kernel::Ket kets[2];

    int thread = device::threads::rank();
    int block = device::threads::size();

    __shared__ double *Cps[2][2];

    __shared__ Int4 index[2];
    __shared__ int R[2];

    if (thread%warpSize == 0) {

	int offset = 0;
	kernel::Bra bra = kernel::Bra(quartet);
	kernel::Ket ket = kernel::Ket(quartet);

	for (int i = 0; i < 4; ++i) {
	    double *c = NULL;
	    if (SP_ & (1 << i)) {
		c = shmem_ + offset;
		offset += (i/2 == 0) ? bra.K : ket.K;
	    }
	    if (threadIdx.z == 0)
		Cps[i/2][i%2] = c;  
	}

	int roots = N*(warpSize/M)*(L ? 2 : 1); // r/w
	int prims = (bra.size() + ket.size());
	int ints = (warpSize/2)*Braket::size;

	int w = thread/warpSize;
	offset += w*max(roots+prims, ints);
	bra.ptr_ = offset+roots;
	ket.ptr_ = offset+roots+bra.size();
	ptr[w] = offset;

    	bras[w] = bra;
    	kets[w] = ket;

	int r = device::grid::rank();
	r = w + r*2;
	if (r < quartets.size()) {
	    R[w] = device::shuffle<64>(r, quartets.size());
	    index[w] = quartets[R[w]];
	}
	else {
	    R[w] = -1;
	    index[w] = quartets[0];
	}
    }

    int t = thread%warpSize;
    int w = thread/warpSize;

    quartet.centers(index[w].elems, bras[w], kets[w], t);
    __syncthreads();

    // if (w == 0) bras[0].primitives(quartet, Cps[0], t, warpSize);
    // if (w == 1) kets[0].primitives(quartet, Cps[1], t, warpSize);
    if (w == 0) quartet.primitives<2>(bras, Cps[0], t, warpSize);
    if (w == 1) quartet.primitives<2>(kets, Cps[1], t, warpSize);
    __syncthreads();

    if (R[threadIdx.z] >= 0) {
	kernel::Bra &bra = bras[w];
	kernel::Ket &ket = kets[w];
	double *shmem = shmem_+ptr[w];
	(*this)(bra, ket, Cps, shmem, threadIdx);
	transform(Braket(), R[w], index[w].elems,
		  shmem, shmem + Braket::size,
		  t, warpSize);
    }

    // for (int i = 0; i < 2; ++i) {
    // 	if (R[i] < 0) continue;
    // 	__syncthreads();
    // 	double *shmem = shmem_+ptr[i];
    // 	transform(Braket(), R[i], index[i].elems,
    // 		  shmem, shmem + Braket::size);
    // }

    }

    struct vector_difference {
	const double *u;
	const double *v;
	__device__
	vector_difference(const double *u, const double *v)
	    : u(u), v(v) {}
	__device__ double operator[](int i) const {
	    return u[i] - v[i];
	}
    };

__device__
void operator()(const kernel::Bra &bra, const kernel::Ket &ket,
		double* (&Cps)[2][2], double *shmem,
		dim3 thread) {

    namespace device = ::cuda::device;

    double I[Braket::size] = {};

    for (int K = thread.y; K < quartet.K(); K += warpSize/M) {
	if (thread.x >= N) break; // idle thread

	int kl = K/bra.K;
	int ij = K - kl*bra.K;

	const double &A = bra.A(ij);
	const double &B = ket.B(kl);

	//double AB = 1/(A+B);
	{
	    double Cmax = 1.0;
	    double C = rsqrt(A+B)*(bra.e(ij)*ket.e(kl));
	    if (SP_ & 0x1) Cmax = max(Cmax, fabs(Cps[0][0][ij]));
	    if (SP_ & 0x2) Cmax = max(Cmax, fabs(Cps[0][1][ij]));
	    if (SP_ & 0x4) Cmax = max(Cmax, fabs(Cps[1][0][kl]));
	    if (SP_ & 0x8) Cmax = max(Cmax, fabs(Cps[1][1][kl]));
	    if (fabs(Cmax*C) < cutoff) continue;
	}

	double *W = shmem + thread.y*N;
	double *t2 = W + (warpSize/M)*N;

	{
	    double AB = 1/(A+B);
	    double rho = A*B*AB;
	    double X = rho*device::distance2(bra.rA(ij), ket.rB(kl));
	    rysq::roots<(L ? N : 0)>(X, t2, W, thread.x);
	    if (L > 0) t2[thread.x] *= AB;
	}   

	double c = (bra.e(ij)*ket.e(kl))*rsqrt(A+B);
	double C[Braket::ket::nc][Braket::bra::nc];

#pragma unroll
	for (int l = 0; l < Braket::D::nc; ++l) {
#pragma unroll
	for (int k = 0; k < Braket::C::nc; ++k) {
	    double Ckl = c;
	    if (SP_ & 0x4 && k == 0) Ckl *= Cps[1][0][kl];
	    if (SP_ & 0x8 && l == 0) Ckl *= Cps[1][1][kl];
#pragma unroll
	    for (int j = 0; j < Braket::B::nc; ++j) {
#pragma unroll
	    for (int i = 0; i < Braket::A::nc; ++i) {
		double Cij = 1;
		if (SP_ & 0x1 && i == 0) Cij *= Cps[0][0][ij];
		if (SP_ & 0x2 && j == 0) Cij *= Cps[0][1][ij];
		C[k+l*Braket::C::nc][i+j*Braket::A::nc] = Cij*Ckl;
	    }
	    }
	}
	}

	// vector_difference rAi(bra.rA(ij), bra.r[0]);
	// vector_difference rBk(ket.rB(ij), ket.r[0]);
	// vector_difference rAB(bra.rA(ij), ket.rB(kl));

	double rAB[3], rAi[3], rBk[3];
	for (int i = 0; i < 3; ++i) {
	    rAi[i] = bra.rA(ij)[i] - bra.r[0][i];
	    rBk[i] = ket.rB(kl)[i] - ket.r[0][i];
	    rAB[i] = bra.rA(ij)[i] - ket.rB(kl)[i];
	}
		           
	kernel_::eval<1>(A, B, rAi, rBk, rAB, bra.dr, ket.dr,
			 t2+thread.x, W+thread.x, C, I);

    }

    int t = thread.x + thread.y*M;

    // sum warp
    {
    	double *Q = shmem + (t/2)*Braket::size;
	if (t%2 == 0)
    	for (int i = 0; i < Braket::size; ++i) {
	    Q[i] = I[i];
	    // double q = I[i];
    	    // if (t%2 == 0) Q[i] = q;//I[i];
	    // if (t%2 != 0) Q[i] += q;//
    	}
	if (t%2 != 0)
    	for (int i = 0; i < Braket::size; ++i) {
    	    Q[i] += I[i];
    	}
    }

    for (int i = t; i < Braket::size; i += warpSize) {
    	for (int j = 1; j < warpSize/2; ++j) {
    	    shmem[i] += shmem[i + j*Braket::size];
    	}
    }

    for (int i = t; i < Braket::size; i += warpSize) {
    	shmem[i+Braket::size] = shmem[i];
    }

    int *index = (int*)shmem;
    if (t%warpSize == 0) {	
	kernel_::index(index);
    }

    for (int i = t; i < Braket::size; i += warpSize) {
	int j = index[i];
    	shmem[j] = shmem[i+Braket::size];
    }

}

static cudaFuncAttributes attributes() {
    return global::function<quadrature2, Transform>::attributes();
}

};


template<class BRAKET, class Transform>
std::ostream& operator<<(std::ostream &os,
			 const quadrature2<BRAKET, Transform> &Q) {
    BOOST_AUTO(const &attributes, Q.attributes());
    os << "Kernel: "
       << "regs: " << attributes.numRegs << ", "
       << "smem: " << attributes.sharedSizeBytes << ", " 
       << "(" << Q.quartet[0] << Q.quartet[1]
       << Q.quartet[2] << Q.quartet[3] << ") "
       << " K: "
       << Q.quartet[0].K() << ":" << Q.quartet[1].K() << ":"
       << Q.quartet[2].K() << ":" << Q.quartet[3].K();
	//<< " " << Q << " shared: " << shared
    return os;
}


template<rysq::type A, rysq::type B, rysq::type C, rysq::type D,
	 class Transform, class = void>
struct Kernel2 : Eri <Transform> {
    Kernel2(const Shell::Quartet &quartet, const Context &Context) {
	throw kernel::kernel_not_implemented();
    }
};

template<rysq::type A, rysq::type B, rysq::type C, rysq::type D, class Transform>
struct Kernel2<A, B, C, D, Transform,
	      typename boost::enable_if<
		  quadrature2<meta::braket<A,B,C,D>, Transform> >::type>
: Eri <Transform> {

    typedef meta::braket<A,B,C,D> Braket;
    typedef quadrature2<Braket, Transform> Quadrature;
    static const int N = Quadrature::N;
    static const int M = Quadrature::M;
    static const int warpSize = 32;

    static Dim3 block(const Shell::Quartet &quartet) {
	Dim3 block(M, warpSize/M, 2);
	return block;
    }

    static size_t shared(const Shell::Quartet &quartet,
			 const Context &context,
			 Dim3 block) {
	using std::max;

	kernel::Bra bra(quartet);
	kernel::Ket ket(quartet);

	size_t roots = N*(block.size()/Quadrature::M);
	roots *= (Quadrature::L > 0 ? 2 : 1);

	size_t primitives = 2*(bra.size() + ket.size());

	size_t reduce = Braket::size*(block.size()/2);

	size_t shared = max(reduce, roots+primitives);

	// // hybrid s/p coefficients
	shared += bra.K*(int(A == rysq::SP) + int(B == rysq::SP));
	shared += ket.K*(int(C == rysq::SP) + int(D == rysq::SP));

	shared *= sizeof(double);

	return shared;
    }

    Kernel2(const Shell::Quartet &quartet, const Context &context)
	: quartet_(quartet)
    {
	block_ = block(quartet);
	shared_ = shared(quartet, context, block_);
	//BOOST_PROFILE_LINE;
	// if (quartet.K() > 1) throw std::exception();
	// if (quartet.ket().hybrid())  throw std::exception();
	// std::cout <<  "new kernel " << M << std::endl;
    }

    void operator()(const detail::Centers &centers,
		    const detail::Quartets &quartets,
		    const cuda::Eri::Parameters &p,
		    cudaStream_t stream,
		    Transform transform) {

	BOOST_PROFILE_LINE;

	kernel::Quartet quartet(quartet_, centers);
	BOOST_VERIFY(quartet.size() == Braket::size);

	{
	    Dim3 grid = Dim3((quartets.size()+1)/2);
	    Quadrature quadrature(quartet, quartets, p.cutoff);
	    // std::cout << "launch: " << block_.x << " "
	    // 	      << "dynamic shmem: " << shared_ << " "
	    // 	      << quadrature << std::endl;
	    CUDA_VERIFY
		(global(grid, block_, shared_, stream)(quadrature, transform));
	    // CUDA_VERIFY(cudaStreamSynchronize(stream));
	}
	
    }

private:
    Dim3 block_;
    size_t shared_;
    Shell::Quartet quartet_;
};


} // namespace kernel
} // namespace cuda
} // namespace rysq


#endif // RYSQ_CUDA_KERNEL_QUADRATURE2_HPP
