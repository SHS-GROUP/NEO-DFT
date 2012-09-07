#ifndef RYSQ_CUDA_KERNEL_ERI_HPP
#define RYSQ_CUDA_KERNEL_ERI_HPP

#include "cuda/kernel/kernel.hpp"
#include "cuda/kernel/device.hpp"

namespace rysq {
namespace cuda {
namespace kernel {
namespace eri {

struct Transform {
    typedef detail::Eri::List List;
    List eri;

    // __host__ __device__
    Transform() : eri() {}
    explicit Transform(List eri): eri(eri) {}

    template<class BRAKET>
    static int block(const Context &context) {
	return 1;
    }

    template<class BRAKET>
    static size_t shmem(const Context &context, dim3 block) {
	return 0;
    }

    template<typename T, class BRAKET>
    __device__
    void operator()(BRAKET, const int &index, const T (&quartet)[4],
		    const double *Q, double *shmem,
		    int thread = ::cuda::device::threads::rank(),
		    int num_threads = ::cuda::device::threads::size()) {
	apply(Q, BRAKET::size, this->eri[index], thread, num_threads);
    }

    __device__
    static void apply(const double *Q, int size, double *eri,
		      int thread, int num_threads) {
	//__syncthreads();
	for (int i = thread; i < size; i += num_threads) {
	    eri[i] = rysq::SQRT_4PI5*Q[i];
	}
    }

};

} // namespace eri
} // namespace kernel
} // namespace cuda
} // namespace rysq

#endif // RYSQ_CUDA_KERNEL_ERI_HPP
