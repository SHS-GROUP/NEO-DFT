#ifndef RYSQ_CUDA_KERNEL_GLOBAL_HPP
#define RYSQ_CUDA_KERNEL_GLOBAL_HPP

#include <cuda.h>
#include <cuda_runtime.h>
#include <boost/preprocessor/cat.hpp>
#include <boost/static_assert.hpp>

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
    // __launch_bounds__(64,4)
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
	    static cudaError_t set_cache_config(enum cudaFuncCache config) {
		typedef void (*type)(F, A0);
		return cudaFuncSetCacheConfig(type(&global__<F,A0>), config);
	    }
	};

	template<class F>
	struct function<F>;

	template<class F, class A>
	static cudaError_t set_cache_config(enum cudaFuncCache cache_config) {
	    return cudaFuncSetCacheConfig
		((void(*)(F, A))(&global__<F,A>), cache_config);
	}

	global(dim3 grid, dim3 block, size_t shared, cudaStream_t stream)
	    : grid_(grid), block_(block), shared_(shared), stream_(stream) {}

	template<class F>
	// typename boost::enable_if_c<(sizeof(F) <= 256)>::type
	cudaError_t operator()(F &f) const {
	    BOOST_STATIC_ASSERT(sizeof(F) <= 1024);
	    global__<F><<<grid_, block_, shared_, stream_>>>(f);
	    return cudaGetLastError();
	}

	template<class F, class A>
	cudaError_t operator()(F &f, A &a) const {
	    //BOOST_STATIC_ASSERT(sizeof(F) + sizeof(A) <= 1024);
	    //set_cache_config<F,A>(cudaFuncCachePreferL1);
	    //set_cache_config<F,A>(cudaFuncCachePreferShared);
	    global__<F,A><<<grid_, block_, shared_, stream_>>>(f, a);
	    return cudaGetLastError();
	}

	template<class F, class A0, class A1>
	cudaError_t operator()(F &f, A0 &a0, A1 &a1) const {
	    BOOST_STATIC_ASSERT(sizeof(F) + sizeof(A0) + sizeof(A1) <= 1024);
	    global__<F,A0,A1><<<grid_, block_, shared_, stream_>>>(f, a0, a1);
	    return cudaGetLastError();
	}

    };

}
}
}

#endif // RYSQ_CUDA_KERNEL_GLOBAL_HPP
