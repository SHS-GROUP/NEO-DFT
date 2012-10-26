#ifndef TENSOR_STORAGE_HPP
#define TENSOR_STORAGE_HPP

#include "tensor/index.hpp"
#include "tensor/lambda.hpp"
// #include "tensor/functional.hpp"

#include <boost/multi_array.hpp>

#include <boost/mpl/identity.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/mpl/count_if.hpp>
#include <boost/mpl/eval_if.hpp>
#include <boost/mpl/size.hpp>
#include <boost/mpl/assert.hpp>

#include <boost/fusion/include/transformation.hpp>
#include <boost/fusion/include/intrinsic.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/boost_array.hpp>
#include <boost/fusion/include/io.hpp>
#include <boost/type_traits.hpp>


namespace tensor {
namespace storage {

    template<class A, class enable = void>
    struct array_traits;

    template<typename T, size_t N = 0>
    struct array {
	typedef boost::multi_array<T, N> type;
    };

    template<class A>
    struct array<A> : array_traits<A> {};

    template<class A>
    struct array_ref {
	typedef typename array<A>::array_ref type;
	static type generate(A &data, int i) {
	    return array<A>::template generate<type>(data, i);
	}
    };

    template<class A>
    struct const_array_ref {
	typedef typename array<A>::const_array_ref type;
	static type generate(A &data, int i) {
	    return array<A>::template generate<type>(data, i);
	}
    };


    template<class A, size_t N>
    struct array_view {
	typedef typename array<A>::template view<N>::type type;
	template<class R>
	static type generate(A &data, const R &ranges) {
	    return array<A>::template generate_view<type>(data, ranges);
	}
   };


    template<class A, size_t N>
    struct const_array_view {
	typedef typename array<A>::template const_view<N>::type type;
	template<class R>
	static type generate(const A &data, const R &ranges) {
	    return array<A>::template generate_view<type>(data, ranges);
	}
    };

    template<class A>
    typename array<A>::pointer
    data(A &a) {
	return array<A>::data(a);
    }

    template<class A>
    typename array<A>::iterator
    begin(A &a) {
	return array<A>::begin(a);
    }

    template<class A>
    typename array<A>::iterator
    end(A &a) {
	return array<A>::end(a);
    }


    template<class A>
    typename array<A>::size_type
    size(const A &a, size_t i) {
	return array<A>::size(a, i);
    }

    template<class A>
    typename array<A>::stride_type
    stride(const A &a, size_t i) {
	return array<A>::stride(a, i);
    }


}
}

#endif // TENSOR_STORAGE_HPP
