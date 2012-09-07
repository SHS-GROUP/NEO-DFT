#ifndef BOOST_CUDA_RUNTIME_HPP
#define BOOST_CUDA_RUNTIME_HPP

#include <string>
#include <vector>

#include "boost/cuda/forward.hpp"
#include "boost/cuda/device_ptr.hpp"
#include "boost/cuda/mapped_ptr.hpp"
#include "boost/cuda/stream.hpp"
#include "boost/cuda/detail/runtime.hpp"

#include <cuda_runtime.h>

namespace boost {
namespace cuda {

    inline void initialize(int device = 0) {
	//std::cout << "initialize " << device << std::endl;
	detail::runtime::initialize(device);
	// std::cout << "initialize " << int(f) << std::endl;
	//detail::runtime::set(f);
    }

    // inline bool is_active() {
    // 	return detail::runtime::is_active();
    // }

    inline void reset() {
	//std::cout << "reset" << std::endl;
	detail::runtime::reset();
    }

    inline std::vector<int> devices(float capability = 0) {
	const std::vector<int> &d = detail::runtime::devices(capability);
	std::vector<int> devices;
	for (size_t i = 0; i < d.size(); ++i) {
	    devices.push_back(d.at(i));
	}
	return devices;
    }

    struct thread {
	static const int disabled = -1;
	explicit thread(int device, flags::flag f = flags::flag(0)) {
	    enabled_ = (device != disabled);
	    if (enabled_) {
		initialize(device);
		boost::cuda::flags::set(f);
	    }
	}
	//~thread() { if (enabled_) thread::exit(); }
	bool enabled() const { return enabled_; }
	static void synchronize() {
	    detail::runtime::thread_synchronize();
	}
	// static void exit() {
	//     reset();
	// }
    private:
	bool enabled_;
    };


    inline void flags::set(flag value) {
	detail::runtime::set_flags(value);
    }

    inline void cache::set(config value) {
	detail::runtime::set_cache_config(value);
    }


    namespace detail {

	/**
	 * void type traits
	 * non-void type traits are partially specialized.
	 * @tparam T datatype
	 * @tparam R reference type
	 */
	template<typename T>
	struct type_traits {
	    typedef T& reference;
	    static const size_t size = sizeof(T);
	};

	/**
	 * non-void type traits.
	 * @tparam T datatype
	 */
	template<>
	struct type_traits<const void> {
	    static const size_t size = 1; /**< type size */
	};

	/**
	 * non-void type traits.
	 * @tparam T datatype
	 */
	template<>
	struct type_traits<void> {
	    static const size_t size = 1; /**< type size */
	};

	template<typename T>
	size_t sizeof_() { return type_traits<T>::size; }

	template<typename T>
	size_t sizeof_(size_t size) { return sizeof_<T>()*size; }

    } // namespace detail


    template<typename T>
    device_ptr<T> device_malloc(size_t size) {
	size_t bytes = detail::sizeof_<T>(size);
	T* ptr = reinterpret_cast<T*>(detail::runtime::malloc_device(bytes));
	return device_ptr<T>::reinterpret(ptr);
    }

    template<typename T>
    void device_free(device_ptr<T> ptr) {
	detail::runtime::free_device(ptr.get());
    }

    template<typename T>
    mapped_ptr<T> mapped_malloc(size_t size) {
	size_t bytes = detail::sizeof_<T>(size);
	T *ptr = reinterpret_cast<T*>(detail::runtime::malloc_mapped(bytes));
	return mapped_ptr<T>::reinterpret(ptr);
    }

    template<typename T>
    void mapped_free(mapped_ptr<T> ptr) {
	detail::runtime::free_host(ptr.get());
    }

    template<typename T>
    void copy(const T *from, device_ptr<T> to, size_t size,
	      const stream &stream = synchronous) {
	detail::runtime::copy(from, to, detail::sizeof_<T>(size), stream);
    }

    template<typename T>
    void copy(const T *begin, const T *end, device_ptr<T> to,
	      const stream &stream = synchronous) {
	size_t size = end - begin;
	detail::runtime::copy(begin, to, detail::sizeof_<T>(size), stream);
    }


    template<typename T>
    void copy(device_ptr<const T> from, T *to, size_t size,
	      const stream &stream = synchronous) {
	detail::runtime::copy(from, to, detail::sizeof_<T>(size), stream);
    }

    template<typename T>
    void copy(const device_ptr<T> &from, T *to, size_t size,
	      const stream &stream = synchronous) {
	detail::runtime::copy(from, to, detail::sizeof_<T>(size), stream);
    }

    template<typename T>
    void copy(device_ptr<const T> from, const std::string &symbol, size_t size,
	      stream stream = synchronous) {
	detail::runtime::copy(from, symbol, detail::sizeof_<T>(size), stream);
    }

    template<typename T>
    void fill(device_ptr<T> ptr, size_t size, char value) {
	detail::runtime::fill(ptr, size*sizeof(T), value);
    }

    template<typename T>
    struct device_reference {
	static device_reference bind(device_ptr< T> data) {
	    return device_reference(data);
	}
	operator T() const {
	    return get(data_);
	} 
	device_reference& operator=(const T&value) {
	    set(data_, value);
	    return *this;
	} 
    private:
	device_ptr<T> data_;
	device_reference(device_ptr< T> data) : data_(data) {}
    };


}
}

#endif /* BOOST_CUDA_RUNTIME_HPP */
