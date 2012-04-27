#ifndef RYSQ_CUDA_KERNEL_ERI_HPP
#define RYSQ_CUDA_KERNEL_ERI_HPP

#include "boost/cuda/device/device.hpp"
#include "cuda/kernel/kernel.hpp"

namespace rysq {
namespace cuda {
namespace kernel {
namespace eri {

struct Transform {
    typedef double* const* List;
    List eri;
    size_t offset;

    __host__ __device__ Transform() : eri(NULL) {}
    explicit Transform(List eri): eri(eri), offset(0) {}

    template<typename T>
    __device__
    void operator()(int index, const T (&quartet)[4],
		    const double *Q, size_t size,
    		    const thread_block &block, double *shmem) const {
	apply(Q, size, eri[offset + index], block);
    }

    template<int M, typename T, typename R>
    __device__    
    void operator()(int index, const T (&quartet)[4],
		    const R (&Q)[M], size_t size,
		    const thread_block &block, double *shmem) const {
	apply<M>(Q, size, eri[offset + index], block);
    }

    __device__
    static void apply(const double *Q, size_t size, double *eri,
		      const thread_block &block) {
	__syncthreads();
	for (size_t i = block.rank; i < size; i += block.size) {
	    eri[i] = rysq::SQRT_4PI5*Q[i];
	}
    }

    template<int M, typename R>
    __device__    
    static void apply(const R (&Q)[M], size_t size, double *eri,
		      const thread_block &block) {
	__syncthreads();
	eri += block.rank;
#pragma unroll
	for (int i = 0; i < M; ++i) {
	    if (block.rank + i*block.size < size) *eri = rysq::SQRT_4PI5*Q[i];
	    eri += block.size;
	}
    }
    // template<>
    // __device__    
    // static void apply<1>(const double (&Q)[1], size_t size, double *eri,
    // 			 const thread_block &block) {
    // 	if (block.rank < size) eri[block.rank] = rysq::SQRT_4PI5*Q[0];
    // }
};

} // namespace eri
} // namespace kernel
} // namespace cuda
} // namespace rysq

#endif // RYSQ_CUDA_KERNEL_ERI_HPP
