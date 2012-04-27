#ifndef BOOST_FUSION_REPLACE_AT_HPP
#define BOOST_FUSION_REPLACE_AT_HPP

#include <boost/fusion/include/iterator.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/mpl/contains.hpp>
#include <boost/typeof/typeof.hpp>

namespace boost {
namespace fusion {

    namespace result_of {

	template<class S, class M, class T>
	struct replace_at {
	    typedef typename advance<typename begin<S>::type, M>::type I;

	    typedef typename as_vector<
		typename erase<S, I>::type
	    		>::type const S_;
	    typedef typename insert<
		S_, 
		typename advance<typename begin<S_>::type, M>::type,
		T
		>::type type;
	};

    } // namespace result_of


    template<class M, class S, class T>
    typename result_of::replace_at<const S, M, T>::type
    replace_at(const S &s, const T &t) {
	BOOST_AUTO(const &s_, as_vector(erase(s, advance<M>(begin(s)))));
	return (insert(s_, advance<M>(begin(s_)), t));
    }

}
}
	
#endif // BOOST_FUSION_REPLACE_AT_HPP


