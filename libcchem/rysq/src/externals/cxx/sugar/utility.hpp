#ifndef _SUGAR_MINMAX_HPP_
#define _SUGAR_MINMAX_HPP_

#include <algorithm>

namespace sugar {

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

#endif /* _SUGAR_MINMAX_HPP_ */
