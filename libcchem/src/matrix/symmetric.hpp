#ifndef _MATRIX_SYMMETRIC_HPP_
#define _MATRIX_SYMMETRIC_HPP_

#include <algorithm>

namespace matrix {

    static int triangular_index(int i, int j) { return i + (j*j + j)/2; }

    template<typename T>
    struct symmetric_adapter {
	typedef size_t size_type;
	typedef int index_type;
	typedef T& reference;
	typedef const T& const_reference;
	symmetric_adapter(T *data, size_type N)
	    : data_(data), size_(N) {}
	size_type size1() const { return size_; }
	size_type size2() const { return size_; }
	reference operator()(index_type i, index_type j) {
	    return data_[triangular_index(std::min(i,j), std::max(i,j))];
	}
	const_reference operator()(index_type i, index_type j) const {
	    return data_[triangular_index(std::min(i,j), std::max(i,j))];
	}
    private:
	T *data_;
	size_type size_;
    };

    template<typename T>
    struct const_symmetric_adapter {
	typedef size_t size_type;
	typedef int index_type;
	typedef T& reference;
	typedef const T& const_reference;
	const_symmetric_adapter(const T *data, size_type N)
	    : data_(data), size_(N) {}
	size_type size1() const { return size_; }
	size_type size2() const { return size_; }
	const_reference operator()(index_type i, index_type j) const {
	    return data_[triangular_index(std::min(i,j), std::max(i,j))];
	}
    private:
	const T *data_;
	size_type size_;
    };

}

#endif /* _MATRIX_SYMMETRIC_HPP_ */
