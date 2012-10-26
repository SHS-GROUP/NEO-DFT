#ifndef BOOST_PP_IS_ITERATING

#ifndef TENSOR_PERMUTE_HPP
#define TENSOR_PERMUTE_HPP

#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/preprocessor/cat.hpp>

#include "tensor/config.hpp"
#include "tensor/index.hpp"
#include "tensor/view.hpp"
#include "tensor/exchange.hpp"

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/replace.hpp>
#include <boost/mpl/for_each.hpp>
#include <boost/mpl/sort.hpp>
#include <boost/mpl/unique.hpp>
#include <boost/mpl/comparison.hpp>


namespace tensor {
namespace detail {

    struct permute_impl {

	template<class TU, class T, class U>
	struct swap {
	    struct _ {};
	    typedef typename boost::mpl::replace<TU, T, _>::type _U;
	    typedef typename boost::mpl::replace<_U, U, T>::type _T;
	    typedef typename boost::mpl::replace<_T, _, U>::type type;
	};

	template<class S, int K = boost::mpl::size<S>::value-1>
	struct pairs {
	    typedef boost::mpl::integral_c<int,K> J;
	    typedef boost::mpl::integral_c<
		int,
		boost::mpl::distance<
		    typename boost::mpl::begin<S>::type,
		    typename boost::mpl::find<S, J>::type
		    >::type::value
		> I;
	    typedef typename boost::mpl::push_front<
		typename pairs<typename swap<S, I, J>::type, K-1>::type,
		boost::mpl::pair<I,J>
		>::type type;
	};
	
	template<class S>
	struct pairs<S,0> {
	    typedef boost::mpl::vector<> type;
	};

	template<class T>
	struct functor {
	    T &a;
	    template<class I>
	    void operator()(boost::mpl::pair<I,I>) const { }		
	    template<class I, class J>
	    void operator()(boost::mpl::pair<I,J>) const {
		exchange<I::value, J::value>(a);
	    }		
	};

	template<class S, class T>
	static void apply(T &a) {
	    typedef typename pairs<S>::type P;
	    functor<T> f = { a };
	    boost::mpl::for_each<typename pairs<S>::type>(f);
	}

    };

} // namespace tensor
} // namespace detail

// generate specializations
#define BOOST_PP_ITERATION_LIMITS (2, TENSOR_MAX_RANK)
#define BOOST_PP_FILENAME_1       "tensor/permute.hpp" // this file
#include BOOST_PP_ITERATE()

#endif // TENSOR_PERMUTE_HPP


#else // BOOST_PP_IS_ITERATING

#define N BOOST_PP_ITERATION()

namespace tensor {

    template<BOOST_PP_ENUM_PARAMS(N, int I), class A, class R>
    void permute(tensor_view<A,R> a) {
	namespace mpl = boost::mpl;
	typedef tensor_view<A,R> T;
	typedef typename mpl::vector_c<int, BOOST_PP_ENUM_PARAMS(N, I)>::type S;
	typedef typename mpl::unique<typename mpl::sort<S>::type,
	    boost::is_same<mpl::_1,mpl::_2> >::type U;

	BOOST_MPL_ASSERT_MSG((N == T::rank),
			     BOOST_PP_CAT(TENSOR_PERMUTE_EXPECTED_RANK_, N),
			     (T, S));
	
	BOOST_MPL_ASSERT_MSG(((mpl::size<U>::value == N) && 
			      (mpl::front<U>::type::value == 0) &&
			      (mpl::back<U>::type::value == N-1)),
			     TENSOR_PERMUTE_INVALID_INDICES_,
			     (T, S));

	detail::permute_impl::apply<S>(a);
    }

}


#undef N

#endif // BOOST_PP_IS_ITERATING
