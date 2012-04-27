#ifndef BOOST_CUDA_DEVICE_PTR_HPP
#define BOOST_CUDA_DEVICE_PTR_HPP

#include <cstdlib>
#include <iterator>

namespace boost {
namespace cuda {


    namespace detail {

    /**
     * @brief device pointer wrapper
     * @tparam T datatype
     */
	template<typename T, template<typename> class P>
	struct ptr_base {
	    typedef std::random_access_iterator_tag iterator_category;
	    typedef T value_type;
	    typedef ptrdiff_t difference_type;
	    typedef P<T> pointer;
	    struct traits {
		typedef P<T> pointer;
		typedef P<const T> const_pointer;
	    };
	    T* get() const { return data_; }
	    P<T> operator+(size_t offset) const {
		return P<T>::reinterpret(data_ + offset);
	    }
	    P<T> operator-(int offset) const {
		return P<T>::reinterpret(data_ - offset);
	    }
	    ptrdiff_t operator-(const P<T> &rhs) const {
		return this->data_ - rhs.data_;
	    }
	    bool operator!=(const P<T> &rhs) const {
		return this->data_ != rhs.data_;
	    }
	    bool operator==(const P<T> &rhs) const {
		return this->data_ == rhs.data_;
	    }
	    P<T>& operator++() {
		++data_;
		return *static_cast<P<T>*>(this);
	    }
	    P<T>& operator--() {
		--data_;
		return *static_cast<P<T>*>(this);
	    }
	    operator bool() const { return data_; }
	protected:
	    ptr_base(T *data = NULL) : data_(data) {}
	    T* data_;		/**< device pointer */
	};

    } // namespace detail

    template<typename T>
    struct device_ptr;

    template<typename T>
    struct device_ptr_base : detail::ptr_base<T, device_ptr> {
	operator T*() { return this->get(); } 
	static device_ptr<T> reinterpret(T *data) {
	    device_ptr<T> ptr;
	    ptr.data_ = data;
	    return ptr;
	}
    protected:
	device_ptr_base(T *data = NULL) {
	    this->data_ = data;
	}
    };

    /**
     * @brief device pointer wrapper
     * @tparam T datatype
     */
    template<typename T>
    struct device_ptr : device_ptr_base<T> { };


    /**
     * @brief device pointer wrapper
     * @tparam void datatype
     */
    template<>
    struct device_ptr<const void> : device_ptr_base<const void> {
	device_ptr() {}
	template<typename T>
	device_ptr(const device_ptr<T> &ptr)
	    : device_ptr_base<const void>(ptr.get()) {}
    };

    /**
     * @brief device pointer wrapper
     * @tparam void datatype
     */
    template<>
    struct device_ptr<void> : device_ptr_base<void> {
	device_ptr() {}
	template<typename T>
	device_ptr(const device_ptr<T> &ptr)
	    : device_ptr_base<void>(ptr.get()) {}
    private:
	template<typename T>
	device_ptr(const device_ptr<const T> &ptr);
    };

    /**
     * @brief device pointer wrapper
     * @tparam T datatype
     */
    template<typename T>
    struct device_ptr<const T> : device_ptr_base<const T> {
	device_ptr() {}
	device_ptr(const device_ptr<T> &ptr)
	    : device_ptr_base<const T>(ptr.get()) {}
    };

    template<typename T, typename U>
    device_ptr<T> reinterpret(device_ptr<U> ptr) {
	return device_ptr<T>::reinterpret(reinterpret_cast<T*>(ptr.get()));
}
    template<typename T>
    T* pointer(device_ptr<T> ptr) {
	return ptr.get();
    }

}
}


#endif // BOOST_CUDA_DEVICE_PTR_HPP
