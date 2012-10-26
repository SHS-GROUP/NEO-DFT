#ifndef RYSQ_CUDA_KERNEL_QUADRATURE_HPP
#define RYSQ_CUDA_KERNEL_QUADRATURE_HPP

#include <exception>
#include <cuda_runtime.h>

#include <boost/typeof/typeof.hpp>
#include "boost/utility/profiler.hpp"

#include "assert.hpp"
#include "roots/roots.hpp"
#include "meta.hpp"

#include "cuda.hpp"
#include "cuda/kernel/device.hpp"
#include "cuda/detail.hpp"
#include "cuda/kernel/kernel.hpp"
#include "cuda/kernel/quartet.hpp"
#include "cuda/kernel/recurrence.hpp"
#include "cuda/kernel/transfer.hpp"
#include "cuda/kernel/global.hpp"


namespace rysq {
namespace cuda {
namespace kernel {


template<int K>
struct power_of_two {
    template<int V, bool>
    struct evaluate {
	enum { value = (evaluate<2*V, (2*V < K)>::value) };
    };
    template<int V>
    struct evaluate<V, false> {
	enum { value = V };
    };
    enum { value = (evaluate<1, (1 < K)>::value) };
};


template<class BRAKET, int M_, class Transform, typename T = double>
struct quadrature1 {

    typedef kernel::Bra Bra;
    typedef kernel::Ket Ket;

    typedef typename BRAKET::Bra BRA;
    typedef typename BRAKET::Ket KET;
    static const int N = (BRAKET::L)/2 + 1;

    static const int M = (M_ > 0 ? M_ : -M_);
    static const bool use_warp = (M_ < 0);

    // known at compile time, correct up to eight roots
    static const int P = power_of_two<N>::value;

    static const int TRANSFER =
	(((BRAKET::A::L && BRAKET::B::L) << 0) |
	 ((BRAKET::C::L && BRAKET::D::L) << 1));

    static const int SP = 
	(((BRAKET::A::type == rysq::SP) << 0) |
	 ((BRAKET::B::type == rysq::SP) << 1) |
	 ((BRAKET::C::type == rysq::SP) << 2) |
	 ((BRAKET::D::type == rysq::SP) << 3));

    struct Ints2d {
	int size;
	static const ushort stride =
	    ((BRA::A::L+1)*(BRA::B::L+1)*
	     (KET::A::L+1)*(KET::B::L+1)*3);
	__device__ T* x(int w) { return x_ + w*size; }
	__device__ T* y(int w) { return y_ + w*size; }
	__device__ T* z(int w) { return z_ + w*size; }
	__device__ T* g(int w) { return g_ + w*size; }
	__device__ T* tmp(int w) { return tmp_ + w*size; }
	__device__ Ints2d() {}
	BOOST_GPU_ENABLED
	Ints2d(const Bra &bra, const Ket &ket, T *data) {
	    const int mi = BRAKET::A::L+1;
	    const int mj = BRAKET::B::L+1;
	    const int nk = BRAKET::C::L+1;
	    const int nl = BRAKET::D::L+1;
	    const int n = nk+nl+1;
	    size = mi*mj*nk*nl;
	    x_ = data + 0*size;
	    y_ = data + 1*size;
	    z_ = data + 2*size;
	    size *= 3*N;
	    if (TRANSFER == 0x0) {
		g_ = x_;
		tmp_ = NULL;
	    }
	    else if (TRANSFER == 0x3) {
		g_ = x_;
		tmp_ = data + size;
		size += mi*mj*n*3*N;
	    }
	    else {
		g_ = data + size;
		tmp_ = NULL;
		size += mi*mj*n*3*N;
	    }
	}
    private:
	T *g_, *x_, *y_, *z_, *tmp_;
    };

    struct Team {
	static const int size =
	    power_of_two<(BRAKET::Bra::L+(BRAKET::Ket::L > 0))*3*N>::value;
	//int count, rank;
	static const bool warp = quadrature1::use_warp;
	__device__ operator int() const {
	    return warp ? threadIdx.y : 0;
	}
	__device__ int count() const {
	    return warp ? blockDim.y : 1;
	}
	__device__ void synchronize() const {
	    if (!warp) __syncthreads();
	}
    };


    kernel::Quartet quartet;
    const rysq::Int4 *quartets;
    const ushort4 *index;
    float cutoff;

    quadrature1(kernel::Quartet quartet, const rysq::Int4 *quartets,
		const ushort4 *index, double cutoff)
	: quartet(quartet),
	  quartets(quartets),
	  index(index),
	  cutoff(cutoff/(rysq::SQRT_4PI5*quartet.K()*1e2))
    {
	RYSQ_ASSERT(index);
    }
    
    __device__ void operator()(Transform &transform) {

	namespace device = ::cuda::device;

    __shared__ Bra bra;
    __shared__ Ket ket;
    __shared__ Ints2d ints2d;

    extern __shared__ double shmem_[];

    int block = device::threads::size();

    struct {
	int rank;
	int3 recurrence1, recurrence2;
	int3 transfer1;
	struct transfer2 : int3 { int1 dim; } transfer2;
	__device__
	operator int() { return rank; }
    } thread;
    thread.rank = device::threads::rank();

    __shared__ Int4 index;
    __shared__ int R;

    __shared__ double *t2W;
    __shared__ double* Cps[2][2];

    if (thread == 0)  {

	size_t offset = 0;
	bra = Bra(quartet, shmem_+offset);
	offset += bra.size();
	ket = Ket(quartet, shmem_+offset);
	offset += ket.size();

	t2W = shmem_ + offset;
	offset += 2*N*quartet.K();	

	for (int i = 0; i < 4; ++i) {
	    double *c = NULL;
	    if (SP & (1 << i)) {
		c = shmem_ + offset;
		offset += (i/2 == 0) ? bra.K : ket.K;
	    }
	    Cps[i/2][i%2] = c;	    
	}

	ints2d = Ints2d(bra, ket, shmem_+offset);

	R = device::grid::shift_rank();
	index = quartets[R];
    }

    quartet.centers(index.elems, bra, ket, thread);
    __syncthreads();

    {
	if (thread/warpSize == 0)
	    bra.primitives(quartet, Cps[0], thread%warpSize, warpSize);
	if (thread/warpSize == 1 || block <= warpSize)
	    ket.primitives(quartet, Cps[1], thread%warpSize, warpSize);
    }
    __syncthreads();


    // partition block into P parts
    for (int K = thread/P; K < quartet.K(); K += block/P) {
    	int kl = K/bra.K;
    	int ij = K - kl*bra.K;
    	int p = thread%P;
    	double *t2 = t2W + 2*N*K;
    	double *W = t2 + N;
    	const double &A = bra.A(ij);
    	const double &B = ket.B(kl);

	//double AB = 1/(A+B);
	double C = rsqrt(A+B)*(bra.e(ij)*ket.e(kl));

	{
	    double Cmax = 1.0;
	    if (SP & 0x1) Cmax = max(Cmax, fabs(Cps[0][0][ij]));
	    if (SP & 0x2) Cmax = max(Cmax, fabs(Cps[0][1][ij]));
	    if (SP & 0x4) Cmax = max(Cmax, fabs(Cps[1][0][kl]));
	    if (SP & 0x8) Cmax = max(Cmax, fabs(Cps[1][1][kl]));
	    if (p < N) W[p] = 0;
	    if (fabs(Cmax*C) < cutoff) continue;
	}

    	if (p < N) {
	    double AB = 1/(A+B);
	    double rho = A*B*AB;
	    double X = rho*device::distance2(bra.rA(ij), ket.rB(kl));
	    rysq::roots<N>(X, t2, W, p);
	    t2[p] *= AB;
	    W[p] *= C;
	}
    }

    __syncthreads();

    Team team;

    if (team.warp) {
	thread.rank = thread%Team::size;
	block = Team::size;
    }

    {
	int t = thread;

	thread.recurrence1.x = 0;
	thread.recurrence1.y = t%3;
	thread.recurrence1.z = t/3;

	if (team.warp) {
	    thread.recurrence2.x = (t%(BRA::L+1));
	    thread.recurrence2.y = (t/(BRA::L+1))%3;
	    thread.recurrence2.z = (t/(BRA::L+1))/3;
	    //thread.transfer1 = thread.recurrence2;
	    thread.transfer1.x = (t%(BRA::L));
	    thread.transfer1.y = (t/(BRA::L))%3;
	    thread.transfer1.z = (t/(BRA::L))/3;
	}
	else {
	    const int m = BRA::L+1;
	    thread.recurrence2.x = (t%m);
	    thread.recurrence2.y = (t%power_of_two<3*m>::value)/m;
	    thread.recurrence2.z = (t/power_of_two<3*m>::value);
	    thread.transfer1.x = (t%BRA::L);
	    thread.transfer1.y = (t%power_of_two<3*BRA::L>::value)/BRA::L;
	    thread.transfer1.z = (t/power_of_two<3*BRA::L>::value);
	}

	if (TRANSFER & 0x2) {
	    const int n = block/power_of_two<3*N>::value;
	    thread.transfer2.dim.x = n;
	    thread.transfer2.x = (t%n);
	    thread.transfer2.y = (t/n)%3;
	    thread.transfer2.z = (t/n)/3;
	    if (thread.transfer2.x >= n) {
		thread.transfer2.z = N;
	    }
	}

    }

    T Q[M];
    for (int i = 0; i < M; ++i) {
    	Q[i] = 0.0;
    }

    // ket contractions
    for (int K = team; K < quartet.K(); K += team.count()) {

	double *t2 = t2W + 2*N*K;
	double *W = t2 + N;

	if (W[0] == 0) continue;

    	int kl = K/bra.K;
    	int ij = K - kl*bra.K;
    
	const double &A = bra.A(ij);
	const double &B = ket.B(kl);

	T *G = ints2d.g(team);

	team.synchronize();

	recurrence1(N, BRA::L+1, KET::L+1,
		    //bra.m, ket.n,
		    bra.A1(ij), B, bra.r[0], ket.r[0],
		    bra.rA(ij), ket.rB(kl),
		    t2, W, G,
		    thread.recurrence1);
	team.synchronize();

	recurrence2(N, BRA::L+1, ket.L+1, //KET::L+1,
		    // ket.n, //KET::L+1,
		    // bra.m, ket.n,
		    A, ket.B1(kl), bra.r[0], ket.r[0],
		    bra.rA(ij), ket.rB(kl),
		    t2, W, G,
		    thread.recurrence2);
	team.synchronize();

	if (TRANSFER & 0x1) {
	    team.synchronize();
	    T *I_ = (TRANSFER == 0x3) ? ints2d.tmp(team) : ints2d.x(team);
	    transfer1(N, BRA::A::L+1, BRA::B::L+1, KET::L+1,
		      //bra.mi, bra.mj, ket.n,
		      bra.dr, G, I_, thread.transfer1);
	}

	if (TRANSFER & 0x2) {
	    T *Ix = ints2d.x(team);
	    G = (TRANSFER == 0x3) ? ints2d.tmp(team) : ints2d.g(team);
	    team.synchronize();
	    transfer2(N, (BRA::A::L+1)*(BRA::B::L+1),
		      KET::A::L+1, KET::B::L+1,
		      ket.dr, G, Ix,
		      (int3)thread.transfer2, thread.transfer2.dim);
	}

	team.synchronize();

	T *Ix = ints2d.x(team);
	T *Iy = ints2d.y(team);
	T *Iz = ints2d.z(team);

#pragma unroll
	for (int i = 0; i < M; ++i) {
	    ushort j = thread + i*block;
	    ushort4 idx = { 0, 0, 0, 0 };
	    if (j < quartet.size()) idx = this->index[j];
	    T q = 0.0;
#pragma unroll
	    for (int a = 0; a < N; ++a) {
		q += Ix[idx.x]*Iy[idx.y]*Iz[idx.z];  
		idx.x += ints2d.stride;
		idx.y += ints2d.stride;
		idx.z += ints2d.stride;
	    }
	    //q *= C;
	    const ushort &sp = idx.w;
	    if ((SP & 0x1) && (sp & 0x1)) q *= Cps[0][0][ij];
	    if ((SP & 0x2) && (sp & 0x2)) q *= Cps[0][1][ij];
	    if ((SP & 0x4) && (sp & 0x4)) q *= Cps[1][0][kl];
	    if ((SP & 0x8) && (sp & 0x8)) q *= Cps[1][1][kl];
	    Q[i] += q;
	}

    }
    __syncthreads();

    if (team == 0) {
    	for (int i = 0; i < M; ++i) {
    	    ushort j = thread + i*block;
    	    if (j < quartet.size()) shmem_[j] = Q[i];
    	}
    }

    for (int k = 1; k < team.count(); ++k) {
	__syncthreads();
	if (team != k) continue;
#pragma unroll
	for (int i = 0; i < M; ++i) {
	    ushort j = thread + i*block;
	    if (j < quartet.size()) shmem_[j] += Q[i];
	}
    }

    __syncthreads();
    
    transform(BRAKET(), R, index.elems, shmem_, shmem_ + quartet.size());

}

static cudaFuncAttributes attributes() {
    return global::function<quadrature1, Transform>::attributes();
}

};


template<class BRAKET, int M, class Transform, typename T>
std::ostream& operator<<(std::ostream &os,
			 const quadrature1<BRAKET, M, Transform, T> &Q) {
    BOOST_AUTO(const &attributes, Q.attributes());
    os << "Kernel: "
       << "regs: " << attributes.numRegs << ", "
       << "smem: " << attributes.sharedSizeBytes << ", " 
       << "ints/thread: " << Q.M << ", "
       << "use warp: " << Q.use_warp << ", "
       << "(" << Q.quartet[0] << Q.quartet[1]
       << Q.quartet[2] << Q.quartet[3] << ") "
       << Q.SP << " "
       << " K: "
       << Q.quartet[0].K() << ":" << Q.quartet[1].K() << ":"
       << Q.quartet[2].K() << ":" << Q.quartet[3].K();
	//<< " " << Q << " shared: " << shared
    return os;
}


template<rysq::type A, rysq::type B, rysq::type C, rysq::type D>
struct Braket1 : meta::braket<A,B,C,D> {
    static const int size = meta::braket<A,B,C,D>::size;
    static const int L = meta::braket<A,B,C,D>::L;
    static const int N = L/2+1;
    typedef typename meta::braket<A,B,C,D>::bra Bra;
    typedef typename meta::braket<A,B,C,D>::ket Ket;
};

template<rysq::type A, rysq::type B, rysq::type C, rysq::type D,
	 class Transform, class T = void>
struct Kernel1 : kernel::Eri<Transform> {
    Kernel1(const Shell::Quartet quartet, const Context &context) {
	throw kernel::kernel_not_implemented();
    }
};

template<rysq::type A, rysq::type B, rysq::type C, rysq::type D>
struct enable1 {
    typedef Braket1<A,B,C,D> braket;
    static const bool value =
	((A != 0) &&
	 (Braket1<A,B,C,D>::N <= 5) &&
	 (braket::A::L >= braket::B::L) && 
	 (braket::C::L >= braket::D::L) && 
	 (braket::Bra::L >= braket::Ket::L)//  &&
	 //  (Braket1<A,B,C,D>::size == 24) ||
	 //  (Braket1<A,B,C,D>::size == 36) ||
	 //  (Braket1<A,B,C,D>::size == 64) ||
	 //  (Braket1<A,B,C,D>::size == 96) ||
	 //  (Braket1<A,B,C,D>::size == 144) ||
	 //  (Braket1<A,B,C,D>::size == 216) ||
	 //  (Braket1<A,B,C,D>::size == 256) ||
	 //  (Braket1<A,B,C,D>::size == 384)||
	 //  (Braket1<A,B,C,D>::size == 576)||
	 //  (Braket1<A,B,C,D>::size == 864)||
	 //  (Braket1<A,B,C,D>::size == 1296)
	 //  )
	 );
};

template<rysq::type A, rysq::type B, rysq::type C, rysq::type D,
	 class Transform>
struct Kernel1<A, B, C, D, Transform,
	       typename boost::enable_if< enable1<A,B,C,D> >::type
> : kernel::Eri<Transform> {

    typedef Braket1<A,B,C,D> BRAKET;

    static const int warpSize = 32;
    static const int N = BRAKET::N;

    static const int Mask = (((A!=0) << 0) | ((B!=0) << 1) |
			     ((C!=0) << 2) | ((D!=0) << 3));

    static const int TEAM = quadrature1<BRAKET, 1, Transform>::Team::size;

    typedef quadrature1<BRAKET,
			-int((BRAKET::size + TEAM - 1)/(TEAM)),
			Transform> Quadrature1;
    typedef quadrature1<BRAKET,
			(BRAKET::size + 2*warpSize - 1)/(2*warpSize),
			Transform> Quadrature2;
    typedef quadrature1<BRAKET,
			(BRAKET::size + 3*warpSize - 1)/(3*warpSize),
			Transform> Quadrature3;
    typedef quadrature1<BRAKET,
			(BRAKET::size + 4*warpSize - 1)/(4*warpSize),
			Transform> Quadrature4;
    typedef quadrature1<BRAKET,
    			(BRAKET::size + 8*warpSize - 1)/(8*warpSize),
    			Transform> Quadrature8;

    static bool use_warp(const Shell::Quartet &quartet) {
	return (quartet.K() > 1 && (TEAM <= warpSize));
    }

    template<class F>
    static Dim3 block(const Context &context,
		      const Shell::Quartet &quartet,
		      bool team = false) {
	using std::max;

	int transform = Transform::template block<BRAKET>(context);

	if (team) {
	    int nt = align(transform, 2*warpSize);
	    Dim3 block(1,1,1);
	    block.x = TEAM;
	    block.y = nt/TEAM;
	    return block;
	}

	const int m = BRAKET::Bra::L + 1;
	const int n = BRAKET::Ket::L + 1;

	Dim3 block(1,1,1);
	block.x = max<int>(power_of_two<3*m>::value, 3*n)*N;
	block.x = max<int>(block.x, (quartet.size()+F::M-1)/F::M);
	block.x = max<int>(block.x, transform);
	block.x = max<int>(block.x, 2*warpSize);
	block.x = align(block.x, warpSize);
	return block;
    }

    template<class F>
    static size_t shared(const Context &context,
			 const Shell::Quartet &quartet,
			 const Dim3 &block,
			 bool team = false) {
	kernel::Bra bra(quartet);
	kernel::Ket ket(quartet);
	size_t shared = bra.size() + ket.size();
	shared += 2*N*quartet.K(); // roots and weights
	shared += bra.K*(int(A == rysq::SP) + int(B == rysq::SP));
	shared += ket.K*(int(C == rysq::SP) + int(D == rysq::SP));
	size_t ints2 = (typename F::Ints2d(bra, ket, NULL).size);
	if (team) ints2 *= (block.size())/TEAM;
	shared += ints2;
	size_t transform = Transform::template shmem<BRAKET>(context, block);
	shared = std::max(shared, quartet.size()+transform);
	shared *= sizeof(double);
	return shared;
    }

    template<class F>
    static int occupancy(const Context &context,
			 const Dim3 &block, int shared) {
	BOOST_AUTO(const &props, context.properties());
	int threads = block.size();
	int registers = F::attributes().numRegs;
	int shmem = F::attributes().sharedSizeBytes + shared;
	int occupancy = std::min<int>(props.regsPerBlock/(threads*registers),
				      props.sharedMemPerBlock/shmem);
	// printf("occupancy = %i, threads = %i, registers = %i, shmem = %i, regsPerBlock = %i, sharedMemPerBlock = %i\n", occupancy, threads, registers, shmem, 
	//        props.regsPerBlock, props.sharedMemPerBlock);
	return occupancy;
    }

    Kernel1(const Shell::Quartet &quartet, const Context &context)
	: quartet_(Shell::Quartet(quartet))
    {
	//if (quartet[1] != 0) throw kernel::kernel_not_implemented();

	this->index_ = (const ushort4*)context.index
	    (quartet[0], quartet[1], quartet[2], quartet[3]);

	if (use_warp(quartet)) {
	    block_ = block<Quadrature2>(context, quartet, true);
	    shared_ = shared<Quadrature2>(context, quartet, block_, true);
	    this->Q = 1;
	    int occ = this->occupancy<Quadrature2>(context, block_, shared_);
	    if (occ == 0) throw kernel::kernel_not_implemented();
	    return;
	}

	this->Q = 0;
	this->block_ = Dim3(0,0,0);
	this->shared_ = 0;
	int occupancy = 0;

#define CONFIGURE(TYPE)								\
	{									\
	    typedef Quadrature ## TYPE Quadrature;				\
	    Dim3 block = this->block<Quadrature>(context, quartet);		\
	    int shared = this->shared<Quadrature>(context, quartet, block);	\
	    int occ = this->occupancy<Quadrature>(context, block, shared);	\
	    if (occ >= occupancy) {						\
		occupancy = occ;						\
		block_ = block;							\
		shared_ = shared;						\
		this->Q = TYPE;							\
	    }									\
	}

	CONFIGURE(2);
	CONFIGURE(3);
	CONFIGURE(4);
	CONFIGURE(8);

#undef CONFIGURE

	if (occupancy == 0) throw kernel::kernel_not_implemented();

    }


    void operator()(const detail::Centers &centers,
		    const detail::Quartets &quartets,
		    const cuda::Eri::Parameters &p,
		    cudaStream_t stream,
		    Transform transform) {
	if (this->Q == 1)
	    launch<Quadrature1>(centers, quartets, p, stream, transform);
	if (this->Q == 2)
	    launch<Quadrature2>(centers, quartets, p, stream, transform);
	if (this->Q == 3)
	    launch<Quadrature3>(centers, quartets, p, stream, transform);
	if (this->Q == 4)
	    launch<Quadrature4>(centers, quartets, p, stream, transform);
	if (this->Q == 8)
	    launch<Quadrature8>(centers, quartets, p, stream, transform);
    }

private:

    template<class F>
    void launch(const detail::Centers &centers,
		const detail::Quartets &quartets,
		const cuda::Eri::Parameters &p,
		cudaStream_t stream,
		Transform transform) {
	BOOST_PROFILE_LINE;
	kernel::Quartet quartet(quartet_, centers);	
	Dim3 grid = Dim3(quartets.size());
	F kernel(quartet, quartets.data(), index_, p.cutoff);
	print(kernel);
	// printf("grid(%i,%i,%i), block(%i,%i,%i), shared = %i\n",
	//        grid.x, grid.y, grid.z,
	//        block_.x, block_.y, block_.z,
	//        shared_);
	CUDA_VERIFY(global(grid, block_, shared_, stream)(kernel, transform));
    }

    template<class F>
    void print(const F &f) {
	// std::cout << "launch: ("
	// 	  << block_.x << ":" << block_.y  << ":" << block_.z
	// 	  << ") dynamic shmem: " << shared_ << " "
	// 	  << f << std::endl;
    }

private:
    Dim3 block_;
    size_t shared_;
    Shell::Quartet quartet_;
    const ushort4 *index_;
    int Q;
};



} // namespace kernel
} // namespace cuda
} // namespace rysq


#endif // RYSQ_CUDA_KERNEL_QUADRATURE_HPP
