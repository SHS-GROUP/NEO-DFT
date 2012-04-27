#ifndef GPU_FORWARD_HPP
#define GPU_FORWARD_HPP

#include <cuda_runtime.h>

namespace gpu {

    template<cudaMemcpyKind K>
    struct copy_tag {
	static const cudaMemcpyKind value = K;
    };

    typedef copy_tag<cudaMemcpyHostToHost> host_to_host;
    typedef copy_tag<cudaMemcpyHostToDevice> host_to_device;
    typedef copy_tag<cudaMemcpyDeviceToHost> device_to_host;
    typedef copy_tag<cudaMemcpyDeviceToDevice> device_to_device;

    struct device;

}

#endif // GPU_FORWARD_HPP
