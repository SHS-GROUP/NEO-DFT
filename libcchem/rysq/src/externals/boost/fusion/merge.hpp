#ifndef BOOST_FUSION_MERGE_HPP
#define BOOST_FUSION_MERGE_HPP

#include <boost/fusion/include/remove_if.hpp>
#include <boost/fusion/include/join.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/typeof/typeof.hpp>

namespace boost {
namespace fusion {

    namespace result_of {

	template<class A, class B>
	struct merge {
	    typedef mpl::contains<const A, mpl::_> in_A;
	    typedef typename join<
		const A, const typename remove_if<const B, in_A>::type
		>::type type;
	};

    } // namespace result_of

    template<class A, class B>
    typename result_of::merge<A, B>::type
    merge(const A &a, const B &b) {
	typedef typename result_of::merge<A, B>::in_A in_A;
	return join(a, remove_if<in_A>(b));
    }

}
}
	
#endif // BOOST_FUSION_MERGE_HPP
