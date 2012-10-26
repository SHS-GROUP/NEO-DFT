#ifndef MULTI_ARRAY_HPP
#define MULTI_ARRAY_HPP

#include "multi_array/detail.hpp"

#include <cstdlib>
#include <cassert>
#include <algorithm>

#include <boost/multi_array.hpp>
#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>

 
template<typename T, size_t N>
struct multi_array_view :  detail::multi_array_ref<T,N> {
    typedef detail::multi_array_ref<T,N> base_type;
    typedef size_t size_type;
    multi_array_view(const size_type (&dims)[N], T *data)
	: base_type(dims, data, (size_type*)0) {}
    template<typename U>
    multi_array_view(const size_type (&dims)[N], T *data,
		     const U *strides = NULL)
	: base_type(dims, data, strides) {}
    multi_array_view(boost::detail::multi_array::multi_array_view<T,N> view)
	: base_type(view.shape(), view.origin(), view.strides()) {}
    
    
};


#endif // MULTI_ARRAY_HPP
