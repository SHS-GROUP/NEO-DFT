#ifndef BOOST_CUDA_ALLOCATOR_HPP
#define BOOST_CUDA_ALLOCATOR_HPP

#include <memory>
#include "boost/cuda/device_ptr.hpp"
#include "boost/cuda/mapped_ptr.hpp"
#include "boost/cuda/runtime.hpp"

//#include <boost/thread/thread.hpp>

namespace boost {
namespace cuda {


    template<typename T>
    struct device_allocator {
	typedef device_ptr<T> pointer;
	pointer allocate(size_t size) {
	    pointer p = (size) ? device_malloc<T>(size) : pointer();
	    // if (p) std::cout << boost::this_thread::get_id()
	    // 		     << ": malloc " << p.get() << std::endl;
	    return p;
	}
	void deallocate(pointer ptr) {
	    // if (ptr) std::cout << boost::this_thread::get_id()
	    // 		       << ": free " << ptr.get() << std::endl;
	    if (ptr) device_free(ptr);
	}
    };

    template<typename T>
    struct mapped_allocator {
	typedef T value_type;
	typedef mapped_ptr<T> pointer;
	typedef mapped_ptr<const T> const_pointer;
	typedef T& reference;
	typedef const T& const_reference;
	template<typename U>
	struct rebind {
	    typedef mapped_allocator<U> other;
	};
	detail::mapped_ptr<T>* allocate(size_t size) {
	    void *ptr = (size) ? detail::runtime::malloc_mapped(size*sizeof(T)) : 0;
	    // std::cout << ptr <<  " allocate " <<  size << std::endl;
	    return reinterpret_cast<detail::mapped_ptr<T>*>(ptr);
	}
	void deallocate(pointer ptr, size_t n) {
	    // std::cout << ptr.get() <<  " deallocate " <<  n << std::endl;
	    if (ptr) mapped_free(ptr);
	}
	void construct(T *p, const_reference value) {
	    // std::cout << "construct " <<  p << std::endl;
	  assert(p);
	  new (p) T(value);
	  //std::allocator<T>().construct(p, value);
	}
	void destroy(T *p) {
	    std::allocator<T>().destroy(p);
	}
	size_t max_size() const { return (1u << 30)/sizeof(T); }
    };
    
}
}

#endif // BOOST_CUDA_ALLOCATOR_HPP
