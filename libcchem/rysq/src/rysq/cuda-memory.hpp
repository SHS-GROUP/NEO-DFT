#ifndef _RYSQ_CUDA_MEMORY_HPP_
#define _RYSQ_CUDA_MEMORY_HPP_

#include <stdlib.h>
#include <boost/cuda/runtime.hpp>

namespace rysq {
    namespace cuda {

	using boost::cuda::stream;
	using boost::cuda::synchronous;
	// using boost::cuda::new_;
	// using boost::cuda::delete_;

	// template<typename T, class = void>
	// struct device_ptr {
	//     T *data_;
	//     device_ptr() : data_(NULL) {}
	//     explicit device_ptr(T *data) : data_(data) {}
	//     operator bool()const { return data_; }
	//     T* data() { return data_; }
	//     const T* data() const { return data_; }
	// };

	// template<typename T>
	// struct device_ptr<const T>
	//     : device_ptr <const T, device_ptr<const T> > {
	//     typedef device_ptr<const T, device_ptr> base;
	//     device_ptr() {}
	//     explicit device_ptr(T *data) : base(data) {}
	//     device_ptr(device_ptr<T> ptr) : base(ptr.data()) {} 
	// };

	// template<>
	// struct device_ptr<const void>
	//     : device_ptr <const void, device_ptr<const void> > {
	//     typedef device_ptr<const void,device_ptr> base;
	//     device_ptr() :  base(NULL) {}
	//     explicit device_ptr(void *data) : base(data) {}
	//     template<typename U>
	//     device_ptr(device_ptr<U> ptr) : base(ptr.data()) {}
	// };

	// template<>
	// struct device_ptr<void> : device_ptr <void, device_ptr<void> > {
	//     typedef device_ptr<void,device_ptr> base;
	//     device_ptr() :  base(NULL) {}
	//     explicit device_ptr(void *data) : base(data) {}
	//     template<typename U>
	//     device_ptr(device_ptr<U> ptr) : base(ptr.data()) {}
	// };

	// template<typename T, class U>
	// device_ptr<T,U> operator+(device_ptr<T,U> ptr, size_t size) {
	//     return device_ptr<T,U>(ptr.data() + size);
	// }
	

	// device_ptr<void> malloc(size_t size);
	// void free(device_ptr<void> ptr);

	// void memcpy(device_ptr<const void> from, void *to, size_t size,
	// 	    const stream &stream = synchronous);
	// void memcpy(const void *from, device_ptr<void> to, size_t size,
	// 	    const stream &stream = synchronous);
	// void memset(device_ptr<void> from, char byte, size_t size);

	// template<typename T>
	// device_ptr<T> malloc( size_t size) {
	//     void *data = cuda::malloc(size*sizeof(T)).data();
	//     return device_ptr<T>(reinterpret_cast<T*>(data));
	// }		

	// template<typename T>
	// void free(device_ptr<T> ptr) {
	//     cuda::free(device_ptr<void>(ptr.data()));
	// }		

	// template<typename T>
	// void memcpy(device_ptr<const T> from, T *to, size_t size,
	// 	    const stream &stream = synchronous) {
	//     memcpy(device_ptr<const void>(from), to, size*sizeof(T), stream);
	// }

	// template<typename T>
	// void memcpy(const T *from, device_ptr<T> to, size_t size,
	// 	    const stream &stream = synchronous) {
	//     memcpy(from, device_ptr<void>(to), size*sizeof(T), stream);
	// }

	using boost::cuda:: device_ptr;

    }
}

#endif /* _RYSQ_CUDA_MEMORY_HPP_ */
