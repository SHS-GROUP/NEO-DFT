#ifndef RYSQ_CUDA_KERNEL_GLOBAL_HPP
#define RYSQ_CUDA_KERNEL_GLOBAL_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include "boost/cuda/exception.hpp"
#include <boost/utility/enable_if.hpp>
#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/stringize.hpp>
#include <boost/mpl/assert.hpp>

#ifndef KERNEL_GLOBAL_PREFIX
#error "KERNEL_GLOBAL_PREFIX undefined"
#else
#define GLOBAL(variable) BOOST_PP_CAT(RYSQ_CUDA_KERNEL_, \
				      BOOST_PP_CAT(KERNEL_GLOBAL_PREFIX, variable))
#endif

extern __shared__ double GLOBAL(shmem)[];

#if RYSQ_CUDA_KERNEL_USE_CONSTANT
__constant__ double2 GLOBAL(quartet_data)[36*36 + 4*36];
__constant__ ushort3 GLOBAL(quartet_index2d)[2048];
#define RYSQ_CUDA_KERNEL_QUARTET_DATA(Q) (GLOBAL(quartet_data))
#define RYSQ_CUDA_KERNEL_QUARTET_INDEX2D(Q) (GLOBAL(quartet_index2d))
#else
#define RYSQ_CUDA_KERNEL_QUARTET_DATA(Q) ((const double2*)((Q).data))
#define RYSQ_CUDA_KERNEL_QUARTET_INDEX2D(Q) ((Q).index2d)
#endif

#include <typeinfo>

namespace  rysq {
namespace cuda {
namespace kernel {

    template<class F, class A0, class A1>
    __global__ void
    //__launch_bounds__(100)
	global__(F f, A0 a0, A1 a1) { f(a0, a1); }

    template<class F, class A0>
    __global__ void
    //__launch_bounds__(100)
	global__(F f, A0 a0) { f(a0); }

    template<class F>
    __global__ void
    //__launch_bounds__(100)
	global__(F f) { f(); }

    struct global {
    private:
	dim3 grid_, block_;
	size_t shared_;
	cudaStream_t stream_;
    public:

	template<class F, class A0 = void>
	struct function {
	    typedef cudaFuncAttributes attributes_type;
	    static attributes_type attributes() {
		typedef void (*type)(F, A0);
		attributes_type attributes;
		cudaFuncGetAttributes(&attributes, type(&global__<F,A0>));
		return attributes;
	    }
	};

	template<class F>
	struct function<F>;

	global(dim3 grid, dim3 block, size_t shared, const boost::cuda::stream &stream)
	    : grid_(grid), block_(block), shared_(shared),
	      stream_(stream.data<cudaStream_t>()) {}

	template<class F>
	// typename boost::enable_if_c<(sizeof(F) <= 256)>::type
	void operator()(F &f) const {
	    global__<F><<<grid_, block_, shared_, stream_>>>(f);
	}

	template<class F, class A>
	typename boost::enable_if_c<(sizeof(F) + sizeof(A) <= 256)>::type
	operator()(F &f, A &a) const {
	    global__<F,A><<<grid_, block_, shared_, stream_>>>(f, a);
	}

	template<class F, class A0, class A1>
	typename boost::enable_if_c<(sizeof(F) + sizeof(A0) + sizeof(A1)
				     <= 256)>::type
	operator()(F &f, A0 &a0, A1 &a1) const {
	    global__<F,A0,A1><<<grid_, block_, shared_, stream_>>>(f, a0, a1);
	}

    };

}
}
}

#endif // RYSQ_CUDA_KERNEL_GLOBAL_HPP
