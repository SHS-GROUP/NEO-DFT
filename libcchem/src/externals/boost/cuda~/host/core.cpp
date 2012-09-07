#include "cuda/host/core.hpp"

#include <cuda.h>
#include <cuda_runtime.h>
#include <iostream>

namespace cuda {
    namespace host {

	int initialize(size_t device) {
	    cudaSetDevice(device);
	    return 0;
	}

	namespace detail {

	    void* malloc(size_t size) {
		void *ptr;
		cudaMalloc(&ptr, size);
		return ptr;
	    }

	    void free(void *ptr) {
		cudaFree(ptr);
	    }

	    void copy(device_ptr<const void> from, void *to, size_t size) {
		cudaMemcpy(to, from.data(), size, cudaMemcpyDeviceToHost);
	    }

	    void copy(const void *from, device_ptr<void> to, size_t size) {
		cudaMemcpy(to.data(), from, size, cudaMemcpyHostToDevice);
	    }

	    void copy(device_ptr<const void> from,
		      const std::string &symbol, size_t size) {
		// std::cout << size << " " << from.data() << "\n";
		cudaMemcpyToSymbol(symbol.c_str(), from.data(), size,
				   0, cudaMemcpyDeviceToDevice);
	    }

	    void copy(device_ptr<const void> from,
		      const std::string &symbol, size_t size,
		      stream stream) {
		// std::cout << size << " " << from.data() << "\n";
		//std::cout << stream.data() << std::endl;
		cudaMemcpyToSymbolAsync(symbol.c_str(), from.data(), size,
					0, cudaMemcpyDeviceToDevice,
					static_cast<cudaStream_t>(stream.data()));
	    }

	    void set(device_ptr<void> ptr, char value, size_t size) {
		cudaMemset(ptr.data(), value, size);
	    }

	}

    }
}
