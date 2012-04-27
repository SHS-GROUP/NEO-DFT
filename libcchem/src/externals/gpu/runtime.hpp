#ifndef GPU_RUNTIME_HPP
#define GPU_RUNTIME_HPP

#include <cuda_runtime.h>
#include <stdexcept>

namespace gpu {

    inline void throw_(cudaError_t error) {
	if (error != cudaSuccess) {
	    throw std::runtime_error(cudaGetErrorString(error));
	}
    }

    inline void initialize(int device = 0) {
	gpu::throw_(cudaSetDevice(device));
    }

    inline void shutdown() {
	gpu::throw_(cudaThreadExit());
    }

    struct thread {
	thread() {}
	thread(int device) { gpu::initialize(device); }
	~thread() { } //gpu::shutdown(); }
    };

}

#endif // GPU_GPU_HPP
