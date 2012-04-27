#ifndef CUDA_PROPERTIES_HPP
#define CUDA_PROPERTIES_HPP

#include <cuda_runtime.h>

inline cudaDeviceProp properties() {
    int device;
    cudaGetDevice(&device);
    cudaDeviceProp properties;
    cudaGetDeviceProperties(&properties, device);
    return properties;
}



#endif // CUDA_PROPERTIES_HPP
