#include <cuda.h>
#include <cuda_runtime.h>
#include <assert.h>

//#include "boost/cuda/detail/implementation.hpp"
#include "boost/cuda/runtime.hpp"
#include "boost/cuda/assert.h"

#include <rysq/cuda.hpp>

#include <boost/type_traits/remove_pointer.hpp>

using namespace rysq::cuda;

#define cuda_throw_error( ) {						\
	cudaError_t status = cudaGetLastError();			\
	if (status != cudaSuccess)					\
	    throw std::runtime_error(cudaGetErrorString(status));	\
    }

std::vector<int> rysq::cuda::get_devices() {
    int count;
    cudaGetDeviceCount(&count);
    cuda_throw_error( );    

    std::vector<int> devices;    
    for (int i = 0; i < count; ++i) {
	struct cudaDeviceProp prop;
	cudaGetDeviceProperties(&prop, i);
	if ((prop.major == 1 && prop.minor >= 3) || 
	    (prop.major > 1)) {
	    // std::cout << "device " << i << std::endl;
	    devices.push_back(i);
	}
    }
    return devices;
}

void rysq::cuda::initialize(int device) {
    namespace cuda = boost::cuda;
    cuda::initialize(device);
    // cuda::flags::set(cuda::flags::blocking_sync);
    // cuda::cache::set(cuda::cache::l1);
}

void rysq::cuda::finish() {
    boost::cuda::reset();
}

rysq::cuda::thread::thread(int device) {
#if CUDA_VERSION < 4000
#warning CUDA_VERSION < 4.0:  no shared context 
  rysq::cuda::initialize(device);
#else
    boost::cuda::initialize(device);
#endif
}

//rysq::cuda::thread::~thread() is in quartet.cpp
