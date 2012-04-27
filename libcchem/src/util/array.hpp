#ifndef _UTIL_ARRAY_HPP_
#define _UTIL_ARRAY_HPP_

#include <algorithm>
#include <functional>
#include <iostream>
#include <boost/array.hpp>

namespace util {

    template<typename T, size_t N>
    struct array_base : boost::array<T,N > {
	typedef boost::array<T,N> base;
	typedef T array[N];
	typedef const array const_array;
	array_base() {}
	array_base(T t) { this->assign(t); }
	array_base(const T *p) { *this = p; }
	operator array&() { return this->elems; }
	operator const_array&() const { return this->elems; }
	void operator=(const T &t) { this->assign(t); }
	void operator=(const T *p) { std::copy(p, p+N, this->begin()); }
	template<typename U>
	array_base(const boost::array<U,N> &a) {
	    std::copy(a.begin(), a.end(), this->begin());
	}
    };

    template< typename T, size_t N>
    struct Array : array_base<T,N> {
	typedef array_base<T,N> base;
	Array() {}
	Array(T t) : base(t) {}
	Array(const T *p) : base(p) {}
	template<typename U> Array(const boost::array<U,N> &a) : base(a) {}
    };

    template<typename T>
    struct Array<T,2> : array_base<T,2> {
	static const size_t N = 2;
	typedef array_base<T,2> base;
	Array() {}
	Array(T t) : base(t) {}
	Array(const T *p) : base(p) {}
	Array(T t0, T t1) { T a[] = { t0, t1 }; *this = a; }
	template<typename U> Array(const boost::array<U,N> &a) : base(a) {}
    };

    //     template< typename T>
    //     Array<T,2> operator,(T t0,  T t1) {
    // 	return Array<T,2>(t0, t1);
    //     }


    template<typename T, size_t N, class F>
    Array<T,N> binary_operator(const Array<T,N> &A, const T &t, const F &f) {
	Array<T,N> B;
	for (uint i = 0; i < N; ++i) { B[i] = f(A[i], t); }
	return B;
    }

    template<typename T, typename U, size_t N, class F>
    Array<T,N> binary_operator(const Array<T,N> &A, const Array<U,N> &B, const F &f) {
	Array<T,N> C;
	for (uint i = 0; i < N; ++i) { C[i] = f(A[i], B[i]); }
	return C;
    }

    template<typename T, size_t N>
    Array<T,N> operator-(const Array<T,N> &A, const T &t) {
	return binary_operator(A, t, std::minus<T>());
    }

    template<typename T, size_t N>
    Array<T,N> operator+(const Array<T,N> &A, const T &t) {
	return binary_operator(A, t, std::plus<T>());
    }

    template<typename T, typename U, size_t N>
    Array<T,N> operator+(const Array<T,N> &A, const Array<U,N> &B) {
	return binary_operator(A, B, std::plus<T>());
    }

    template<typename T, typename U, size_t N>
    Array<T,N> operator-(const Array<T,N> &A, const Array<U,N> &B) {
	return binary_operator(A, B, std::minus<T>());
    }



    template<size_t N, typename T>
    Array<T,N> array_of(const T *p) {
	return Array<T,N>(p);
    }

    template<typename T>
    Array<T,2> array_of(const T &t0, const T &t1) {
	return Array<T,2>(t0, t1);
    }

    template<typename T, size_t N>    
    std::ostream& operator<<(std::ostream &os, const Array<T,N> &A) {
	os << "{ ";
	std::copy(A.begin(), A.end() - 1, std::ostream_iterator<T>(os, ", "));
	os << A.back() << " }";
	return os;
    }

}

#endif /* _UTIL_ARRAY_HPP_ */
