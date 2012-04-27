#ifndef BOOST_CUDA_VECTOR_HPP
#define BOOST_CUDA_VECTOR_HPP

/**
 * @file 
 * @brief  Basic CUDA C++ template types and operations
 */

#include <iostream>
#include <vector>
#include <string>

#include <boost/array.hpp>
#include "boost/cuda/runtime.hpp"
#include "boost/cuda/allocator.hpp"

namespace boost {
namespace cuda {

    namespace detail {

	template<typename T, size_t N, class P>
	struct array {
	    typedef typename P::pointer pointer;
	    typedef typename P::const_pointer const_pointer;
	    size_t size() const { return size_; }
	    pointer begin() { return data_; }
	    const_pointer begin() const { return data_; }
	    pointer operator+(size_t offset) {
		return pointer::reinterpret(data_ + offset*N);
	    }
	    const_pointer operator+(size_t offset) const {
		return pointer::reinterpret(data_ + offset*N);
	    }
	protected:
	    array() : data_(), size_(0) {}
	    array(pointer data, size_t size)
		: data_(data), size_(size) {}
	protected:
	    pointer data_;
	    size_t size_;
	};	

	template<class T, class P>
	void copy(const std::vector<T> &from, array<T,1,P> &to) {
	    if (!from.empty())
		copy(&from.front(), to.data(), from.size());
	}
	
	template<class T, size_t N, class P>
	void copy(const std::vector<boost::array<T,N> > &from, array<T,N,P> &to) {
	    if (!from.empty())
		copy(&from.front(), to.data(), N*from.size());
	}

	template<typename T, size_t N, class P>
	void copy(const array<T,N,P> &from, const std::string &symbol,
		  stream stream = synchronous) {
	    copy(from.begin(), symbol, from.size()*N, stream);
	}

    } // namespace detail

    template<typename T, class A = device_allocator<T> >
    struct vector : detail::array<T, 1, typename device_ptr<T>::traits> {
	typedef detail::array<T, 1, typename device_ptr<T>::traits> base;
	typedef A allocator;
	vector() { initialize(0); }
	explicit vector(size_t size) { initialize(size); }
	template<class V>
	explicit vector(const V &v) {
	    initialize(v.size());
	    assign(v);
	}
	vector(const vector &rhs) {
	    initialize(rhs.size());
	    operator=(rhs);
	}
	~vector() {
	    allocator_.deallocate(this->data_);
	}
	vector& operator=(const vector &rhs) {
	    if (this == &rhs) return *this;
	    if (rhs.size())
		throw std::runtime_error("cuda::vector::operator=");
	    this->size_ = 0;
	    return *this;
	}
	template<class A2>
	void assign(const std::vector<T,A2> &v, const stream &stream = synchronous) {
	    resize(v.size());
	    if (!v.empty())
		copy(&v.front(), this->begin(), v.size(), stream);
	}
	void reserve(size_t size){
	    if (this->capacity_ < size) {
		allocator_.deallocate(this->data_);
		this->data_ = allocator_.allocate(size);
		this->capacity_ = size;
	    }
	}
	void resize(size_t size) {
	    reserve(size);
	    this->size_ = size;
	}
	void resize(size_t size, char value) {
	    size_t previous = this->size_;
	    resize(size);
	    if (size > previous) {
		fill(this->data_ + previous, (size - previous), value);
	    }
	}
	void clear() { this->size_ = 0; }
	void reset() {
	    allocator_.deallocate(this->data_);
	    initialize();
	}
    private:
	using base::operator+;
	void initialize(size_t size = 0) {
	    this->data_ = allocator_.allocate(size);
	    this->size_ = size;
	    this->capacity_ = size;
	}
	allocator allocator_;
	size_t capacity_;
    };


	// template<typename T>
	// struct vector {
	//     vector() : size_(0), capacity_(0), data_() {}
	//     explicit vector(const std::vector<T> &data)
	// 	: size_(0), capacity_(0), data_() {
	// 	assign(data);
	//     }
	//     explicit vector(size_t size)
	// 	: size_(0), capacity_(0), data_() {
	// 	resize(size);
	//     }
	//     vector(const vector & rhs) { operator=(rhs); }
	//     vector& operator=(const vector &rhs) {
	// 	if (this == &rhs) return *this;
	// 	if (rhs.data_)
	//     	    throw std::runtime_error("cuda::vector::operator=");
	//     	size_ = 0;
	//     	capacity_ = 0;
	//     	return *this;
	//     }		
	//     ~vector() { if (data_) delete_<T>(data_); }
	//     size_t size() const { return size_; }
	//     void reserve(size_t size){
	// 	//std::cout << "reserve" << std::endl;
	// 	if (capacity_ < size) {
	// 	    if (data_) delete_<T>(data_);
	// 	    data_ = device_ptr<T>();
	// 	    //std::cout << "free" << std::endl;
	// 	}
	// 	if (!data_) {
	// 	    capacity_ = size;
	// 	    data_ = new_<T>(size);
	// 	    //std::cout << " malloc " << size << std::endl;
	// 	}
	//     }
	//     size_t capacity() const { return capacity_; }
	//     void resize(size_t size) {
	// 	reserve(size);
	// 	size_ = size;
	//     }
	//     void resize(size_t size, char value) {
	// 	size_t previous = size_;
	// 	resize(size);
	// 	if (size > previous) {
	// 	    fill(data_ + previous, value, sizeof(T)*(size - previous));
	// 	}
	//     }
	//     void assign(const std::vector<T> &v, const stream &stream = synchronous) {
	// 	resize(v.size());
	// 	if (!v.empty())
	// 	    copy(&v.front(), data_, v.size(), stream);
	//     }
	//     void clear() { size_ = 0 ; }	    
	//     T* data() { return data_.data(); }
	//     const T* data() const { return data_.data(); }
	//     operator device_ptr<T>() { return data_; }
	//     operator device_ptr<const T>() const { return data_; }
	//     operator device_ptr<void>() { return data_; }
	//     operator device_ptr< const void>()  const { return data_; }
	// private:
	//     size_t size_;
	//     size_t capacity_;
	//     device_ptr<T> data_;
	// };


    template<typename T,size_t N = 1>
    struct array
	: detail::array<T, N, typename device_ptr<T>::traits>
    {
    	typedef detail::array<T, N, typename device_ptr<T>::traits> base;
    	array() {}
    	array& operator+(size_t offset) const;
	void reset() { *this = array(); }
    };

    template<typename T, size_t N>
    struct array<const T, N>
	: detail::array<const T, N, typename device_ptr<const T>::traits>
    {
    	typedef detail::array<const T, N, typename device_ptr<const T>::traits> base;
    	array() {}
    	array(device_ptr<const T> data, size_t size) : base(data, size) {}
    	// array(const base &a) : base(a.data(), a.size()) {}
    	array& operator+(size_t offset) const;
	void reset() { *this = array(); }
    };


}
}

#endif /* BOOST_CUDA_VECTOR_HPP */

