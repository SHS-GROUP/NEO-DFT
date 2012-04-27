#ifndef CC_TRIPLES_CLEAR_HPP
#define CC_TRIPLES_CLEAR_HPP

#include <algorithm>
#include <boost/typeof/typeof.hpp>

#ifdef HAVE_CUDA
#include <cuda.h>
#endif // HAVE_CUDA

#include "array/adapter.hpp"

namespace cc {
namespace triples {
namespace detail {

    template<class A>
    struct clear_impl {
	static void apply(A &a, int t, int nt) {
	    int n = a.num_elements();
	    BOOST_AUTO(begin, a.origin() + t*(n/nt));
	    BOOST_AUTO(end, begin + n/nt);
	    if (t == nt-1) end += n%nt;	
	    std::fill(begin, end, 0);
	}
    };


#ifdef HAVE_CUDA

    template<typename T, size_t N>
    struct clear_impl< array::adapter<T, N, array::device_tag> > {
	static void apply(array::adapter<T, N, array::device_tag> &a,
			  int t, int nt) {
	    if (t != 0) return;
	    cudaMemset(a.origin(), 0, (a.num_elements())*sizeof(T));
	}
    };

#endif // HAVE_CUDA


    template<class A>
    void clear(A &a, int t, int nt) {
	clear_impl<A>::apply(a, t, nt);
    }



} // namespace detail
} // namespace triples
} // namespace cc

#endif // CC_TRIPLES_CLEAR_HPP
