#ifndef TENSOR_EXCHANGE_HPP
#define TENSOR_EXCHANGE_HPP

#include "tensor/index.hpp"
#include "tensor/view.hpp"
#include "tensor/swap.hpp"

#include <boost/utility/swap.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/sort.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/assert.hpp>

namespace tensor {
namespace detail {

    struct exchange_impl {

	template<class S, class T, class U>
	static
	typename boost::fusion::result_of::remove<
	    typename boost::fusion::result_of::insert<
		S, typename boost::fusion::result_of::find<S,T>::type,
		U>::type,
	    T>::type
	replace_index(const S &s, const T &t, const U &u) {
	    namespace fusion = boost::fusion;
	    return (fusion::remove
		    (fusion::insert
		     (fusion::size<T>(s), u), t));
	}

	template<int N, class S>
	static
	typename boost::fusion::result_of::deref<
	    typename boost::fusion::result_of::advance_c<
		typename boost::fusion::result_of::begin<const S>::type, N
		>::type
	    >::type
	index_at(const S &s) {
	    namespace fusion = boost::fusion;
	    return (fusion::deref
		    (fusion::advance_c<N>
		     (fusion::begin(s))));
	}

	template<int I, int J, int K,
		 class T, class S, class P>
	static
	typename boost::enable_if_c<(K != I && K != J)>::type
	apply(T &a, const S &indices, const P &p) {
	    namespace fusion = boost::fusion;
	    for (int k = 0; k < int(a.size(K)); ++k) {
		BOOST_AUTO(const &K_, index_at<K>(indices));
		apply<I,J,K-1>(a, replace_index(indices, K_, k), p);
	    }
	}

	template<int I, int J, int K,
		 class T, class S, class P>
	static
	typename boost::enable_if_c<(K == J)>::type
	apply(T &a, const S &indices, const P &p) {
	    namespace mpl = boost::mpl;
	    for (int j = 0; j < int(a.size(J)); ++j) {
		apply<I,J,K-1>(a, indices, j);
	    }
	}
	
	template<int I, int J, int K, class T, class S>
	static typename boost::enable_if_c<K == I>::type
	apply(T &a, const S &indices, int j) {
	    namespace fusion = boost::fusion;
	    typedef typename fusion::result_of::as_vector<S>::type V;
	    BOOST_AUTO(a_, a(::tensor::indices<V>(fusion::as_vector(indices))));
	    for (int i = 0; i < j; ++i) {
		boost::swap(a_[i][j], a_[j][i]);
	    }
	}
    };

}
}

namespace tensor {

    template<int I, int J, class A, class R>
    void exchange(tensor_view<A,R> a) {
	namespace mpl = boost::mpl;

	typedef tensor_view<A,R> T;
	static const int K = T::rank-1;
	static const int I_ = (I < J) ? I : J;
	static const int J_ = (J < I) ? I : J;

	BOOST_MPL_ASSERT_MSG(((I_ >= 0) && (I_ < J_) && (J_ <= K)),
			     TENSOR_EXCHANGE_INVALID_INDICES,
			     (T, mpl::int_<I>, mpl::int_<J>));
	BOOST_ASSERT(a.size(I) == a.size(J));
	detail::exchange_impl::apply<I_,J_,K>(a, a.indices(), mpl::void_());
    }



}

#endif // TENSOR_EXCHANGE_HPP
