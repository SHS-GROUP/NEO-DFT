#ifndef BOOST_FUSION_KEYS_HPP
#define BOOST_FUSION_KEYS_HPP

#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/include/set.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/mpl/transform.hpp>

namespace boost {
namespace fusion {

    namespace result_of {

	template<class M>
	struct keys {
	    typedef typename result_of::as_set<
		typename mpl::transform<M, first<mpl::_> >::type
		>::type type;
	};

    } // namespace result_of

    template<class M>
    typename result_of::keys<M>::type
    keys(const M &map) {
	return as_set(typename result_of::keys<M>::type());
    }

}
}
	
#endif // BOOST_FUSION_KEYS_HPP
