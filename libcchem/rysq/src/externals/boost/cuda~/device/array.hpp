#ifndef BOOST_CUDA_DEVICE_ARRAY_HPP
#define BOOST_CUDA_DEVICE_ARRAY_HPP

#include <cuda.h>
#include <cuda_runtime.h>

// #include <boost/utility/enable_if.hpp>
// #include <boost/type_traits/is_void.hpp>

#include "boost/cuda/runtime.hpp"

namespace boost {
namespace cuda {
namespace device {

	template<typename T, size_t N, class Enable = void>
	struct array_base {
	    typedef cuda::device_ptr<T> device_ptr;
	    typedef cuda::device_ptr<const T> const_device_ptr;
	    array_base() : data_(NULL), size_(0) {}
	    // array_base(size_t size) 
	    // 	: data_(new_<T>(N*size).get()), size_(size) {}
	    array_base(device_ptr data, size_t size) : data_(data.get()), size_(size) {}
	    // array_base( ::cuda::host::array_ref<T,N> &ref);
	    __host__ __device__ size_t size() const { return size_; }
	    __host__ __device__ bool empty() const { return size_ == 0; }
	    __host__ device_ptr data() { return wrap(data_); }
	    __host__ device_ptr data() const { return wrap(data_); }
	protected:
	    T* data_;
	    size_t size_; 
	};	    


	template<typename T, size_t N = 1>
	struct const_array : array_base<const T,N> {
	    typedef array_base<const T,N> base;
	    typedef typename base::device_ptr device_ptr;
	    const_array() {}
	    //explicit const_array(size_t size) : base(size) {}
	    const_array(device_ptr data, size_t size = 0) : base(data, size) {}
	    const_array(array<const T,N> &ref)
		: base(ref.begin(), ref.size()) {}
	    __device__ const T* operator[](int i) const { return base::data_ + i*N; }
	};

	template<>
	struct const_array<int4,1 > : array_base<const int4,1> {
	    typedef array_base<const int4,1> base;
	    typedef base::device_ptr device_ptr;
	    const_array() {}
	    //explicit const_array(size_t size) : base(size) {}
	    const_array(device_ptr data, size_t size = 0) : base(data, size) {}
	    const_array(array<const int,4> &ref)
		: base(reinterpret<const int4>(ref.begin()), ref.size()) {}
	    __device__ const int4& operator[](int i) const { return base::data_[i]; }
	};

	template<typename T>
	struct const_array<T,1> : array_base<const T,1,T&> {
	    typedef array_base<const T,1,T&> base;
	    typedef typename base::device_ptr device_ptr;
	    const_array() {}
	    //explicit const_array(size_t size) : base(size) {}
	    const_array(device_ptr data, size_t size = 0) : base(data, size) {}
	    const_array(array<const T,1> &ref)
		: base(ref.begin(), ref.size()) {}
	    __device__ const T& operator[](int i) const { return base::data_[i]; }
	};

	// template<size_t N >
	// struct device_array< void, N > {
	//     //device_array(host::device_ptr<T> &data) : data_(data) {}
	//     device_array(size_t size = 0) : data_(N*size) {}
	//     __host__ __device__ size_t size() const { return size_; }
	//     __host__ __device__ bool empty() const { return size_ == 0; }
	//     __host__ __device__ operator  void*() { return data_; }
	//     __host__ __device__ operator const void*() const { return data_; }
	// private:
	//     size_t size_;
	//     host::device_ptr<void> data_;
	// };

    
}
}
}

#endif /* BOOST_CUDA_DEVICE_ARRAY_HPP */
