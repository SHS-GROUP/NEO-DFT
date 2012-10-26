#ifndef _CUDAXX_DEVICE_DEVICE_HPP_
#define _CUDAXX_DEVICE_DEVICE_HPP_

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

// #ifndef pow2
// #define pow2(x) ((x)*(x))
// #endif

namespace boost {
namespace cuda {
namespace device {

    template <typename T> class ptr3 {
    public:
	T *x, *y, *z;
	// 	ptr3() : x(NULL), y(NULL), z(NULL) {}
	// 	ptr3(T *x, T *y, T *z) : x(x), y(y), z(z) {}
	// 	ptr3(const ptr3<T> &lv) : x(lv.x), y(lv.y), z(lv.z) {}
	// 	template< typename P> __device__
	// 	ptr3<T>& operator+=(const P &p) {
	// 	    this->x += p;
	// 	    this->y += p;
	// 	    this->z += p;
	// 	    return *this;
	// 	}
	// 	__device__
	// 	ptr3<T> operator+(const ushort3 &v) {
	// 	    return ptr3<T>(this->x + v.x, this->y + v.y, this->z + v.z);
	// 	}
	// 	template< typename P> __device__
	// 	ptr3<T>& operator=(P p) {
	// 	    this->x = (T*)p;
	// 	    this->y = (T*)p;
	// 	    this->z = (T*)p;
	// 	    return *this;
	//	}
    };

    template<typename T> __device__ void swap(T &a, T &b) {
	T q = a;
	a = b;
	b = q;
    }

    template<typename T> __device__
    void copy(T *dest, const T source, int N, int stride = 1, int first = 0) {
	for (int i = first; i < N; i += stride) {
	    dest[i] = source;
	}
    }

    template<typename T> __device__
    void copy(T *dest, const T *source, int N, int stride = 1, int first = 0) {
	for (int i = first; i < N; i += stride) {
	    dest[i] = source[i];
	}
    }

    template<typename T> __device__
    void scale(T *dest, T scale, const T *source, int N, int stride = 1, int first = 0) {
	for (int i = first; i < N; i += stride) {
	    dest[i] = scale*source[i];
	}
    }

    template<typename T> __device__
    void scale(T scale, T *dest, int N, int stride = 1, int first = 0) {
	for (int i = first; i < N; i += stride) {
	    dest[i] *= scale;
	}
    }

    template<typename T> __device__
    T add(int N, const T *v, int stride = 1) {
	T q = 0;
	for (int i = 0; i < N; i += stride) {
	    q += v[i];
	}    
	return q;
    }

    template<typename T> __device__ T distance2(const T *x, const T *y) {
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
	static __device__ 
	int shift_rank() {
	    int rank = grid::rank();
	    int size = grid::size();
	    int i = (rank/4) + (rank%4)*(size/4+(size%4>0));
	    if (rank%4 == 3 && (size%4 == 1 || size%4 == 2)) {
		if (size+1 == i) i /= 2;
		else i -= 1;
	    }
	    return i;
	}
    };

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


    // namespace warp {

    // 	/** 
    // 	 * Reduce local variables in a warp
    // 	 * the reduced value is stored in shmem[warp*16]
    // 	 * 
    // 	 * @param r local variable
    // 	 * @param[out] shmem shared memory
    // 	 */
    // 	template<class T>
    // 	__device__ void reduce(const T &r, T *shmem) {
    // 	    ushort rank = Warp().rank();
    // 	    ushort size = Warp().size();

    // 	    __syncthreads();
    // 	    // sum within the warp using tree 
    // 	    if(rank % 2 == 0) shmem[rank/2] = r;
    // 	    __syncthreads();
    // 	    if(rank % 2 != 0) shmem[rank/2] += r;
    // 	    __syncthreads();
    // 	    if(rank % 4 == 0 && rank+2 < size) {
    // 		shmem[rank/2] += shmem[rank/2+1];
    // 	    }
    // 	    __syncthreads();
    // 	    if(rank % 8 == 0 && rank+4 < size) {
    // 		shmem[rank/2] += shmem[rank/2+2];
    // 	    }
    // 	    __syncthreads();
    // 	    if(rank % 16 == 0 && rank+8 < size) {
    // 		shmem[rank/2] += shmem[rank/2+4];
    // 	    }
    // 	    __syncthreads();
    // 	    if(rank % 32 == 0 && rank+16 < size) {
    // 		shmem[rank/2] += shmem[rank/2+8];
    // 	    }
    // 	}

    // }

    // namespace block {

    // 	/** 
    // 	 * Reduce local variable in BLOCK
    // 	 * 
    // 	 * @param r local value
    // 	 * @param shmem shared memory
    // 	 * 
    // 	 * @return reduced value for rank 0, zero value otherwise
    // 	 */
    // 	template<class T> __device__ T reduce(const T &r, T *shmem){
    // 	    // first reduce without sync
    // 	    warp::reduce(r, &shmem[warp()*16]);
    // 	    __syncthreads();

    // 	    // rank 0 adds across warps
    // 	    T q = T(0);
    // 	    if (rank() == 0) {
    // 		int n = threads::num_warps();
    // 		q = add(n*16, shmem, 16);
    // 	    }
    // 	    __syncthreads();
    // 	    return q;
    // 	}

    // 	template<size_t I>
    // 	struct int_ { static const int value = I; };

    // 	template<size_t I, class T>
    // 	__device__ void reduce(const T (&A)[I], T *shmem, int_<I>) {
    // 	}

    // 	template<size_t I, class T, size_t N>
    // 	__device__ void reduce(const T (&A)[N], T *shmem, int_<I>) {
    // 	    T t = reduce(A[I], shmem + I);
    // 	    if (rank() == 0) shmem[I] = t;
    // 	    reduce(A, shmem, int_<I+1>());
    // 	}


    // 	template<size_t N, typename T>
    // 	__device__ void reduce(const T (&A)[N], T *shmem) {
    // 	    reduce(A, shmem, int_<0>());
    // 	}

    // }


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
}

#endif /* _CUDAXX_DEVICE_DEVICE_HPP_ */

