#ifndef _CXX_ARRAY_HPP_
#define _CXX_ARRAY_HPP_


#include <iostream>
#include <boost/array.hpp>

namespace boost {

    // template<typename  T,size_t N>
    // inline std::ostream& operator<<(std::ostream &os, const  T (&a)[N]) {
    // 	os <<  "{ ";
    // 	for (size_t i = 0; i < N-1; ++i) { os << a[i] << ", "; }
    // 	return os << a[N-1] << " }";
    // }

    template<typename  T,size_t N>
    inline std::ostream& operator<<(std::ostream &os, const boost::array<T,N> &a) {
	os <<  "{ ";
	for (size_t i = 0; i < N-1; ++i) { os << a[i] << ", "; }
	return os << a[N-1] << " }";
    }

    template<typename  T,size_t N>
    boost::array<T,N> make_array(const T (&a)[N]) {
	boost::array<T,N> b;
	for (size_t i = 0; i < N; ++i) { b[i] = a[i]; }
	return b ;
    }

    template<typename  T>
    boost::array<T,3> make_array(const T &a0, const T &a1, const T &a2) {
	boost::array<T,3> b = {{ a0, a1, a2 }};
	return b ;
    }

    template<typename  T>
    boost::array<T,4> make_array(const T &a0, const T &a1, const T &a2, const T &a3) {
	boost::array<T,4> a = {{ a0, a1, a2, a3 }};
	return a;
    }

}

#endif /* _CXX_ARRAY_HPP_ */
