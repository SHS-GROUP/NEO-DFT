#ifndef BOOST_CUDA_MAPPED_PTR_HPP
#define BOOST_CUDA_MAPPED_PTR_HPP

#include "boost/cuda/device_ptr.hpp"
#include "boost/cuda/detail/runtime.hpp"
#include <iostream>
namespace boost {
namespace cuda {
	

    namespace detail {
	template<typename T>
	struct mapped_ptr;
    }

    template<typename T>
    struct mapped_ptr;

    template<typename T, class U>
    struct mapped_ptr_base;

    template<typename T>
    struct mapped_ptr_base<T, void> : detail::ptr_base<T, mapped_ptr>
    {
	operator T*() const { return this->get(); } 
	operator device_ptr<T>() {
	    T *ptr = reinterpret_cast<T*>(detail::runtime::mapped_to_device(*this));
	    return device_ptr<T>::reinterpret(ptr);
	}
	static mapped_ptr<T> reinterpret(T *data) {
	    mapped_ptr<T> ptr;
	    ptr.data_ = data;
	    return ptr;
	}
    protected:
	mapped_ptr_base(T *data) { this->data_ = data; }
	mapped_ptr_base(detail::mapped_ptr<T> *data = NULL) {
	    this->data_ = reinterpret_cast<T*>(data);
	}
    };

    template<typename T, typename U>
    struct mapped_ptr_base : mapped_ptr_base<T, void> {
	typedef mapped_ptr_base<T, void> base_type;
	typedef typename base_type::difference_type difference_type;
	typedef T& reference;
	typedef const U& const_reference;
	reference operator[](difference_type difference) {
	    // std::cout << this->data_ << " operator[]" << std::endl;
	    return *(this->data_ + difference);
	}
	const_reference operator[](difference_type difference) const {
	    return *(this->data_ + difference);
	}
    protected:
	mapped_ptr_base(T *data) : base_type(data) {}
	mapped_ptr_base(detail::mapped_ptr<T> *data = NULL)
	    : base_type(data) {}
    };

    /**
     * @brief device pointer wrapper
     * @tparam T datatype
     */
    template<typename T>
    struct mapped_ptr : mapped_ptr_base<T, T> {
	mapped_ptr(detail::mapped_ptr<T> *data = 0)
	    : mapped_ptr_base<T, T>(data) {}
    };


    /**
     * @brief device pointer wrapper
     * @tparam T datatype
     */
    template<typename T>
    struct mapped_ptr<const T> : mapped_ptr_base<const T, T> {
	mapped_ptr(detail::mapped_ptr<const T> *data = 0)
	    : mapped_ptr_base<const T, T>(data) {}
	mapped_ptr(const mapped_ptr<T> &ptr)
	    : mapped_ptr_base<const T, T>(ptr.get()) {}
    };


    /**
     * @brief device pointer wrapper
     * @tparam void datatype
     */
    template<>
    struct mapped_ptr<const void> : mapped_ptr_base<const void, void> {
	mapped_ptr(detail::mapped_ptr<const void> *data = 0)
	    : mapped_ptr_base<const void, void>(data) {}
	template<typename T>
	mapped_ptr(const mapped_ptr<T> &ptr)
	    : mapped_ptr_base<const void, void>(ptr.get()) {}
    };

    /**
     * @brief device pointer wrapper
     * @tparam void datatype
     */
    template<>
    struct mapped_ptr<void> : mapped_ptr_base<void, void> {
	mapped_ptr(detail::mapped_ptr<void> *data = 0) {}
	template<typename T>
	mapped_ptr(const mapped_ptr<T> &ptr)
	    : mapped_ptr_base<void, void>(ptr.get()) {}
    private:
	template<typename T>
	mapped_ptr(const mapped_ptr<const T> &ptr);
    };


}
}


#endif // BOOST_CUDA_MAPPED_PTR_HPP
