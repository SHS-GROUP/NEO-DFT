#ifndef CUDA_DEVICE_ALGORITHM_HPP
#define CUDA_DEVICE_ALGORITHM_HPP

namespace cuda {
namespace device {

    namespace detail {
	 
	template<int N> struct int_ {};
   
	template<typename T, size_t N, int I>
	__device__
	void fill(T (&A)[N], const T &value, int_<I>) {
	    A[I] = value;
	    fill(A, value, int_<I+1>());
	}

	template<typename T, size_t N>
	__device__
	void fill(T (&A)[N], const T &value, int_<N>) { }

    }


    template<typename T, size_t N>
    __device__
    void fill(T (&A)[N], const T &value) {
	using namespace detail;
	fill(A, value, int_<0>());
    }

}
}

#endif // CUDA_DEVICE_ALGORITHM_HPP
