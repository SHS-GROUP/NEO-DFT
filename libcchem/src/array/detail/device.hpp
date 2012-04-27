#ifndef ARRAY_DETAIL_DEVICE_HPP 
#define ARRAY_DETAIL_DEVICE_HPP

#include "array/detail/utility.hpp"
#include <boost/config.hpp>
#include <boost/typeof/typeof.hpp>

#ifdef __CUDACC__

#include <cuda.h>

namespace array {
namespace detail {

    namespace device {

	typedef ::dim3 dim3;
	typedef ::cudaStream_t stream_t;

	struct grid {
	    static __device__ int x() { return blockIdx.x; }
	    static __device__ int y() { return blockIdx.y; }
	    static __device__ int z() { return 0; }
	    static __device__ int nx() { return gridDim.x; }
	    static __device__ int ny() { return gridDim.y; }
	    static __device__ int nz() { return 1; }
	};

	struct block {
	    static __device__ int x() { return threadIdx.x; }
	    static __device__ int y() { return threadIdx.y; }
	    static __device__ int nx() { return blockDim.x; }
	    static __device__ int ny() { return blockDim.y; }
	};

	template<class F, class A0, class A1>
	__global__
	void launch(F f, A0 a0, A1 a1) {
	    f(a0, a1, device::grid(), device::block());
	}
	
#define CHECK_(F) {							\
	    cudaError_t err = (F);					\
	    if (err != cudaSuccess)					\
		throw std::runtime_error(std::string(__FILE__) + ": " +	\
					 cudaGetErrorString(err));	\
	}

	struct kernel {
	    kernel(dim3 grid = dim3(1), dim3 block = dim3(1),
		   size_t shared = 0, stream_t stream = 0) 
		: grid_(grid), block_(block),
		  shared_(shared), stream_(stream) {}
	    template<class F, class A0, class A1>
	    void operator()(const F &f, const A0 &a0, const A1 &a1) const {
		launch<<< grid_, block_, shared_, stream_ >>>(f, a0, a1);
		CHECK_((cudaGetLastError()));
		//CHECK_((cudaThreadSynchronize()));
	    }
	private:
	    dim3 grid_, block_;
	    size_t shared_;
	    stream_t stream_;
	};

#undef CHECK_

    }

}
}

#endif // __CUDACC__

#endif // ARRAY_DETAIL_DEVICE_HPP
