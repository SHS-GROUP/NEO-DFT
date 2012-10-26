#ifndef BOOST_FUSION_SWAP_HPP
#define BOOST_FUSION_SWAP_HPP

#include <boost/fusion/replace_at.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/typeof/typeof.hpp>

namespace boost {
namespace fusion {

    namespace result_of {

	template<class S, class I, class J>
	struct swap {
	    typedef typename result_of::begin<S>::type S_;
	    typedef typename result_of::distance<S_, I>::type M;
	    typedef typename result_of::distance<S_, J>::type N;

	    typedef typename replace_at<
		S, M, typename value_of<J>::type>::type _;
	    typedef typename replace_at<
		typename as_vector<_>::type,
		N, typename value_of<I>::type
		>::type type;
	};

    } // namespace result_of


    template<class S, class I, class J>
    typename result_of::swap<S, I, J>::type
    swap(const S &s, const I &i, const J &j) {
	typedef typename result_of::begin<S>::type S_;
	typedef typename result_of::distance<S_, I>::type M;
	typedef typename result_of::distance<S_, J>::type N;

	BOOST_AUTO(const &_, replace_at<M>(s, deref(j)));
	return (replace_at<N>(as_vector(_), deref(i)));
    }

}
}
	
#endif // BOOST_FUSION_SWAP_HPP


