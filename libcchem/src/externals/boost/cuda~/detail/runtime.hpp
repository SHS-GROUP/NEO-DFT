#ifndef BOOST_CUDA_DETAIL_RUNTIME_HPP
#define BOOST_CUDA_DETAIL_RUNTIME_HPP

#include <vector>

#include "boost/cuda/forward.hpp"
#include "boost/cuda/device_ptr.hpp"
#include "boost/cuda/stream.hpp"

namespace boost {
namespace cuda {

namespace detail {

    template<class>
    struct runtime_implementation {
	static std::vector<int> devices(float capability = 0);
	static bool is_active();
	static void initialize(int device);
	static void set_flags(int value);
	static void set_cache_config(int value);
	static void thread_synchronize();
	static void reset();
	static void* malloc_device(size_t size);
	static void free_device(void *ptr);
	static void* malloc_mapped(size_t size);
	static void free_host(void *ptr);
	static void* mapped_to_device(void *mapped);
	static const void* mapped_to_device(const void *mapped);
	static void copy(device_ptr<const void> from, void *to, size_t size,
			 stream stream);
	static void copy(const void *from, device_ptr<void> to, size_t size,
			 stream stream);
	static void copy(device_ptr<const void> from, const std::string &symbol,
			 size_t size, const stream &stream);
	static void fill(device_ptr<void> ptr, size_t size, char value);
    private:
	static void copy(const void *from, void *to, size_t size,
			 int kind, stream stream);
    };

    typedef runtime_implementation<void> runtime;


}
}
}

#endif // BOOST_CUDA_DETAIL_RUNTIME_HPP


#ifdef BOOST_CUDA_DETAIL_RUNTIME_IMPLEMENTATION

#warning "boost::cuda::detail::runtime instantiation"

#include <vector>
#include "boost/cuda/device_ptr.hpp"
#include "boost/cuda/stream.hpp"

#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>
#include "boost/cuda/exception.hpp"

namespace boost {
namespace cuda {
namespace detail {

    inline cudaDeviceProp device_properties(int device) {
	cudaDeviceProp prop;
	BOOST_CUDA_CHECK(cudaGetDeviceProperties(&prop, device));
	return prop;
    }

    template<class _>
    std::vector<int> runtime_implementation<_>::devices(float capability) {
	int N;
	BOOST_CUDA_CHECK(cudaGetDeviceCount(&N));
	std::vector<int> devices;
	for (int i = 0; i < N; ++i) {
	    cudaDeviceProp prop = device_properties(i);
	    float c = prop.major + 0.1*prop.minor;
	    if (capability > (prop.major + 0.1*prop.minor)) continue;
	    devices.push_back(i);
	}
	return devices;
    }

    template<class _>
    void runtime_implementation<_>::initialize(int device) {
      //std::cout << "::initialize " << device << std::endl;
	BOOST_CUDA_CHECK(cudaSetDevice(device));
    }

    // template<class _>
    // bool runtime_implementation<_>::is_active() {
    // 	CUdevice device;
    // 	cudaError_enum error = cuCtxGetDevice(&device);
    // 	return (error == CUDA_SUCCESS);
    // }

    template<class _>
    void runtime_implementation<_>::set_flags(int value) {
	if (!value) return;
	int f = 0;
	if (value & flags::schedule_auto) f |= cudaDeviceScheduleAuto;
	if (value & flags::schedule_spin) f |= cudaDeviceScheduleSpin;
	if (value & flags::schedule_yield) f |= cudaDeviceScheduleYield;
	if (value & flags::blocking_sync) f |= cudaDeviceBlockingSync;
	if (value & flags::map_host) f |= cudaDeviceMapHost;
	//std::cout << "flag " << f_ << std::endl;
	BOOST_CUDA_CHECK(cudaSetDeviceFlags(f));
    }

    template<class _>
    void runtime_implementation<_>::set_cache_config(int value) {
	if (!value) return;
	int f = 0;
	if (value & cache::none) f |= cudaFuncCachePreferNone;
	if (value & cache::shared) f |= cudaFuncCachePreferShared;
	if (value & cache::l1) f |= cudaFuncCachePreferL1;
	BOOST_CUDA_CHECK(cudaSetDeviceFlags(f));
    }
    
    template<class _>
    void runtime_implementation<_>::thread_synchronize() {
	BOOST_CUDA_CHECK(cudaThreadSynchronize());
    }
    
    template<class _>
    void runtime_implementation<_>::reset() {
#if CUDA_VERSION >= 4000 
	BOOST_CUDA_CHECK(cudaDeviceReset());
#else
	BOOST_CUDA_CHECK(cudaThreadExit());
#endif
    }
    
    template<class _>
    void* runtime_implementation<_>::malloc_device(size_t size) {
	void *ptr;
	BOOST_CUDA_CHECK(cudaMalloc(&ptr, size));
	return ptr;
    }
    template<class _>
    void runtime_implementation<_>::free_device(void *ptr) {
	cudaFree(ptr);
	check_status();
    }
    
    template<class _>
    void* runtime_implementation<_>::malloc_mapped(size_t size) {
	void *ptr = 0;
	if (size)
	    BOOST_CUDA_CHECK(cudaHostAlloc(&ptr, size, cudaHostAllocMapped));
	// std::cout << "malloc_mapped " << ptr <<  " " << size << std::endl;
	return ptr;
    }

    template<class _>
    void runtime_implementation<_>::free_host(void *ptr) {
	BOOST_CUDA_CHECK(cudaFreeHost(ptr));
    }
    
    template<class _>
    void* runtime_implementation<_>::mapped_to_device(void *mapped) {
	void *device = 0;
	// std::cout <<"mapped_to_device " << mapped <<  std::endl;
	if (mapped)
	    BOOST_CUDA_CHECK(cudaHostGetDevicePointer(&device, mapped, 0));
	return device;
    }
    
    template<class _>
    const void* runtime_implementation<_>::mapped_to_device(const void *mapped) {
	return mapped_to_device(const_cast<void*>(mapped));
    }
    
    template<class _>
    void runtime_implementation<_>::copy(device_ptr<const void> from, void *to, size_t size,
						stream stream) {
	copy(from.get(), to, size, cudaMemcpyDeviceToHost, stream);
	check_status();
    }
    
    template<class _>
    void runtime_implementation<_>::copy(const void *from, device_ptr<void> to, size_t size,
						stream stream) {
	copy(from, to.get(), size, cudaMemcpyHostToDevice, stream);
	check_status();
    }
    
    template<class _>
    void runtime_implementation<_>::copy(device_ptr<const void> from,
						const std::string &symbol,
						size_t size, const stream &stream) {
	if (stream.synchronous()) {
	    BOOST_CUDA_CHECK
		(cudaMemcpyToSymbol(symbol.c_str(), from.get(), size,
				    0, cudaMemcpyDeviceToDevice));
	}
	else {
	    BOOST_CUDA_CHECK
		(cudaMemcpyToSymbolAsync(symbol.c_str(), from.get(), size,
					 0, cudaMemcpyDeviceToDevice,
					 stream.data<cudaStream_t>()));
	}
	check_status();
    }
    
    template<class _>
    void runtime_implementation<_>::fill(device_ptr<void> ptr, size_t size, char value) {
	if (value != 0)
	    throw boost::cuda::exception("fill value != 0");
	BOOST_CUDA_CHECK(cudaMemset(ptr.get(), value, size));
    }

    template<class _>
    void runtime_implementation<_>::copy(const void *from, void *to, size_t size,
						int kind, stream stream) {
	if (stream.synchronous()) {
	    BOOST_CUDA_CHECK
		(cudaMemcpy(to, from, size, cudaMemcpyKind(kind)));
	}
	else {
	    BOOST_CUDA_CHECK
		(cudaMemcpyAsync(to, from, size, cudaMemcpyKind(kind),
				 stream.data<cudaStream_t>()));
	}
	check_status();
    }

    
    template class runtime_implementation<void>;
    typedef runtime_implementation<void> runtime;

}
}
}

#endif // BOOST_CUDA_DETAIL_RUNTIME_IMPLEMENTATION

