#ifndef INTEGRAL_FORWARD_HPP
#define INTEGRAL_FORWARD_HPP

#include "basis/basis.hpp"
#include <boost/array.hpp>

namespace integral {

    typedef int index;
    typedef Basis::Center Center;

    // template<typename T>
    // struct quartet : boost::array<T,4> {
    // 	typedef boost::array<T,4> array_type;
    // 	static array_type
    // 	array(const T &t0, const T &t1, const T &t2, const T &t3) {
    // 	    array_type array = {{ t0, t1, t2, t3 }};
    // 	    return array;
    // 	}
    // 	quartet(const T &t0, const T &t1, const T &t2, const T &t3)
    // 	    : array_type(array(t0, t1, t2, t3)) {}
    // };

    class Integral;

}

#endif // INTEGRAL_FORWARD_HPP
