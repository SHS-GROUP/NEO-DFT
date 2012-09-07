#ifndef BOOST_FUSION_INTERSECTION_HPP
#define BOOST_FUSION_INTERSECTION_HPP

#include <boost/fusion/include/filter_if.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/mpl/contains.hpp>

namespace boost {
namespace fusion {

    namespace result_of {

	template<class A, class B>
	struct intersection {
	    typedef mpl::contains<const A, mpl::_> in_A;
	    typedef typename filter_if<const B, in_A>::type type;
	};

    } // namespace result_of

    template<class A, class B>
    typename result_of::intersection<A, B>::type
    intersection(const A &a, const B &b) {
	typedef typename result_of::intersection<A, B>::in_A in_A;
	return filter_if<in_A>(b);
    }

}
}
	
#endif // BOOST_FUSION_INTERSECTION_HPP
