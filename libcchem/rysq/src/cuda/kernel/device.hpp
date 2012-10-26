#ifndef CUDA_DEVICE_HPP
#define CUDA_DEVICE_HPP

#include <cuda.h>
#include <cuda_runtime.h>


#ifndef __CUDACC__
#warning "emulating __CUDACC__"
#define __constant__
#define __device__
namespace {

    dim3 threadIdx, blockIdx, gridDim, blockDim;
    void __syncthreads();

    template<typename T>
    const T& min(const T &a, const T &b);

}
#endif


namespace cuda {
namespace device {

    template<typename T>
    __device__ __forceinline__
    void swap(T &a, T &b) {
	T q = a;
	a = b;
	b = q;
    }

    template<typename T>
    __device__
    void copy(T *dest, const T source, int N, int stride = 1, int first = 0) {
	for (int i = first; i < N; i += stride) {
	    dest[i] = source;
	}
    }

    template<typename T>
    __device__
    void copy(T *dest, const T *source, int N, int stride = 1, int first = 0) {
	for (int i = first; i < N; i += stride) {
	    dest[i] = source[i];
	}
    }

    template<typename T>
    __device__
    void scale(T *dest, T scale, const T *source, int N, int stride = 1, int first = 0) {
	for (int i = first; i < N; i += stride) {
	    dest[i] = scale*source[i];
	}
    }

    template<typename T>
    __device__
    void scale(T scale, T *dest, int N, int stride = 1, int first = 0) {
	for (int i = first; i < N; i += stride) {
	    dest[i] *= scale;
	}
    }

    template<typename T>
    __device__
    T add(int N, const T *v, int stride = 1) {
	T q = 0;
	for (int i = 0; i < N; i += stride) {
	    q += v[i];
	}    
	return q;
    }

    template<typename T>
    __device__
    T distance2(const T *x, const T *y) {
	T q = 0;
	for (int i = 0; i < 3; ++i) {
	    T dxy = x[i] - y[i];
	    q += dxy*dxy;
	}    
	return q;
    }


    struct grid {
	static __device__
	int rank() { return x() + y()*gridDim.x; }
	static __device__
	int size() { return gridDim.x*gridDim.y; }
	static __device__
	int x() { return blockIdx.x; }
	static __device__
	int y() { return blockIdx.y; }
	static __device__ int shift_rank() {
	    int i = grid::rank();
	    int size = grid::size();
	    const int N = 32;
	    int n = N*(size/N);
	    if (i < n) return (i%N)*(n/N) + i/N;
	    return i;
	}
    };

    template<int N>
    inline __device__ int shuffle(int i, int size) {
	int n = N*(size/N);
	if (i < n) return (i%N)*(n/N) + i/N;
	return i;
    }


    struct Block {
	ushort rank, size;
	__device__
	Block(ushort rank, ushort size) : rank(rank), size(size) { }
    };


    namespace threads {

	/** 
	 * Return BLOCK rank, i.e. thread index
	 * 
	 * 
	 * @return thread index
	 */
	inline __device__
	ushort rank() {
	    return (threadIdx.x + blockDim.x*(threadIdx.y + blockDim.y*threadIdx.z));
	}

	/** 
	 * Return BLOCK size, i.e. total number of threads
	 * 
	 * @return BLOCK size
	 */
	inline __device__
	ushort size() {
	    return (blockDim.x*blockDim.y*blockDim.z);
	}

	inline __device__
	short num_warps() {
	    ushort size = threads::size();
	    return (size+warpSize-1)/warpSize;
	}

	inline __device__
	short warp() { return threads::rank()/warpSize; }

	inline __device__
	short warp_size() {
	    return min(threads::size() - threads::warp()*warpSize, warpSize);
	}	    
    }


    struct Threads : Block {
	__device__
	Threads() : Block(threads::rank(), threads::size()) { }
	static __device__
	short warp() { return threads::warp(); }
	static __device__
	short num_warps() { return threads::num_warps(); }
    };


    struct Warp : Block {
	__device__
	Warp() : Block(threads::rank()%warpSize, threads::warp_size()) {}
    };


    template<typename T>
    __device__
    T inner(int N, const T *v1, const T *v2, int stride = 1) {
	T q = 0;
	for (int i = 0; i < N; i += stride) {
	    q += v1[i]*v2[i];
	}    
	return q;
    }

    namespace reduce {

    	// template<size_t I>
    	// struct int_ { static const int value = I; };

    	// template<size_t I, class F, class T>
    	// __device__
	// void apply(const F &f, const T (&A)[I], T *R, T *tmp, int_<I>) {
    	// }

    	// template<size_t I, class F, class T, size_t N>
    	// __device__
	// void apply(const F &f, const T (&A)[N], T *R, T *tmp, int_<I>) {
    	//     T r = f(A[I], tmp);
	//     if (f.rank == 0) R[I] += r;
    	//     apply(f, A, R, tmp, int_<I+1>());
    	// }

	struct Warp {
	    const int size, rank;
	    void *tmp;
	    __device__
	    Warp(void *tmp)
		: size(warpSize),
		  rank(threads::rank()%warpSize),
		  tmp(tmp) {}
	    template<class T>
	    __device__
	    void operator()(const T &r, T *R) const {
		T *shmem = (T*)tmp;
		// sum within the warp using tree 
		if (rank % 2 == 0) {
		    shmem[rank/2] = r;
		}
		if (rank % 2 != 0) {
		    shmem[rank/2] += r;
		}
		if (rank % 4 == 0) {
		    shmem[rank/2] += shmem[rank/2+1];
		}
		if (rank % 8 == 0) {
		    shmem[rank/2] += shmem[rank/2+2];
		}
		if (rank % 16 == 0) {
		    shmem[rank/2] += shmem[rank/2+4];
		}
		if (rank % 32 == 0) {
		    shmem[rank/2] += shmem[rank/2+8];
		}
		if (rank == 0) *R += *shmem;
	    }
	};

	// template<typename T>
	// __device__
	// T warp(const T &t, size_t size, T *tmp) {
	//     return Warp(size)(t, tmp);
	// }

    	// template<size_t N, typename T>
    	// __device__
	// void warp(const T (&A)[N], T *R, T *tmp) {
    	//     apply(Warp(), A, R, tmp, int_<0>());
    	// }

    	// /** 
    	//  * Reduce local variable in BLOCK
    	//  * 
    	//  * @param r local value
    	//  * @param shmem shared memory
    	//  * 
    	//  * @return reduced value for rank 0, zero value otherwise
    	//  */
	// struct Threads {
	//     size_t rank, size;
	//     __device__
	//     Threads() {
	// 	device::Threads threads;
	// 	rank = threads.rank;
	// 	size = threads.size;
	//     }
	//     template<class T>
	//     __device__
	//     T operator()(const T &r, T *shmem) const {
	// 	// first reduce without sync
	// 	reduce::warp(r, &shmem[threads::warp()*(warpSize/2)]);
	// 	__syncthreads();
	// 	// rank 0 adds across warps
	// 	T q = T(0);
	// 	if (threads::rank() == 0) {
	// 	    int n = threads::num_warps();
	// 	    q = add(n*16, shmem, 16);
	// 	}
	// 	__syncthreads();
	// 	return q;
	//     }	    
	// };

	// template<size_t N, typename T>
	// __device__
	// void threads(const T (&A)[N], T *shmem) {
	//     apply(Threads(), A, shmem, int_<0>());
	// }


    }


    template<typename T>
    __device__
    void copy(int size, const T *from, T *to, const Block &block) {
	for (int i = block.rank; i < size; i += block.size) {
	    to[i] = from[i];
	}
    }

    template<typename T>
    __device__
    void fill(int size, T *to, const T &value, const Block &block) {
	for (int i = block.rank; i < size; i += block.size) {
	    to[i] = value;
	}
    }

} // namespace device
}

#endif /* CUDA_DEVICE_HPP */

