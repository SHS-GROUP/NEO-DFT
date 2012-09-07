#ifndef _RYSQ_QUARTET_HPP_
#define _RYSQ_QUARTET_HPP_

// #include "cxx-sugar/boost_array.hpp"
#include <rysq/core.hpp>

#include "transpose.hpp"
#include "util.h"

namespace rysq {

    template<typename T>
    inline void pack(T const &T0, T const &T1, T const &T2, T const &T3, T (&Tn)[4]) {
	Tn[0] = T0;
	Tn[1] = T1;
	Tn[2] = T2;
	Tn[3] = T3;
    }

    template<class C, typename T>
    inline void unpack(const Quartet<C> &Q, T &T0, T &T1, T &T2, T &T3) {
	T0 = Q.template get<T>(0);
	T1 = Q.template get<T>(1);
	T2 = Q.template get<T>(2);
	T3 = Q.template get<T>(3);
    }

}

#endif /* _RYSQ_QUARTET_HPP_ */
