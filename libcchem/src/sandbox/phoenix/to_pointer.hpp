#ifndef TO_POINTER_HPP
#define TO_POINTER_HPP

#include <boost/spirit/home/phoenix/function/function.hpp>

namespace boost {namespace phoenix {

	struct to_pointer_impl {
	    template<class A>
	    struct result { typedef typename A::value_type *type; };
	    template<class A>
	    typename A::value_type*operator()(A &a) const {
		return  &a.data()[0];
	    }
	};

	function<to_pointer_impl> to_pointer;

    } }

#endif // TO_POINTER_HPP
