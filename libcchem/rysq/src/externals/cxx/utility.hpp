#ifndef _CXX_UTILITY_HPP_
#define _CXX_UTILITY_HPP_

#include <algorithm>

#include "cxx/utility/bitmask.hpp"
#include "cxx/utility/permute.hpp"

namespace  cxx {
    namespace utility {



	template<typename  T>
	T ceiling2(const T &v) {
	    T r = T(1);
	    while (r < v) r *= 2;
	    return r;
	}

	template<typename T>
	T qceiling(const T &m, const T &n) { return m/n + (m%n > 0); }

	template<typename T>
	T max(const T &a, const T &b) {
	    return std::max(a,b);
	}

	template<typename T>
	T max(const T &a, const T &b, const T &c) {
	    return max(max(a,b),c);
	}


    }
}

#endif /* _CXX_UTILITY_HPP_ */
