#ifndef _DISTRIBUTED_TRAITS_HPP_
#define _DISTRIBUTED_TRAITS_HPP_

#include "distributed/forward.hpp"

namespace distributed {

    template<class A>
    struct array_traits;

    template<typename T, size_t N_>
    struct array_traits<array<T, N_> > {
	static const size_t N = N_;
    };

    template<typename T>
    struct array_traits<matrix<T> > {
	static const size_t N = 2;
    };

}

#endif // _DISTRIBUTED_TRAITS_HPP_

