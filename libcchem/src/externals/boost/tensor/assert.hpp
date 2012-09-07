#ifndef TENSOR_ASSERT_HPP
#define TENSOR_ASSERT_HPP

#include "tensor/index.hpp"

#include <boost/type_traits.hpp>
#include <boost/mpl/equal.hpp>
#include <boost/mpl/sort.hpp>
#include <boost/mpl/unique.hpp>
#include <boost/mpl/assert.hpp>

#include <boost/fusion/include/mpl.hpp>

#include <boost/preprocessor/cat.hpp>


namespace tensor {
namespace detail {
namespace assert {

    namespace mpl = boost::mpl;

    template<class A, class B>
    struct same_rank {
	typedef void type;
	BOOST_MPL_ASSERT_MSG((mpl::bool_<(A::rank == B::rank)>::value),
			     TENSOR_NOT_SAME_RANK, (A, B));
    };

    template<class A, class B>
    struct same_keys {
	typedef void type;
	typedef typename A::keys_type A_;
	typedef typename B::keys_type B_;
	BOOST_MPL_ASSERT_MSG((mpl::equal<
			      typename mpl::sort<typename A::keys_type>::type,
			      typename mpl::sort<typename B::keys_type>::type
			      >::value),
			     TENSOR_INDEX_NOT_SAME_KEYS, (A_, B_));
    };

    template<class S>
    struct unique_keys {
	typedef void type;
	typedef typename mpl::sort<S>::type S_;
	typedef typename mpl::unique<
	    S_, boost::is_same<mpl::_1, mpl::_2> >::type U;
	BOOST_MPL_ASSERT_MSG((mpl::equal<S_,U>::value),
			     TENSOR_INDEX_NOT_UNIQUE_KEYS, (S));
    };

}
}
}

#define TENSOR_ASSERT_SAME_RANK(A,B)					\
    typedef typename ::tensor::detail::assert::same_rank<A,B>::type	\
    BOOST_PP_CAT(TENSOR_ASSERT_SAME_RANK_, __LINE__) 

#define TENSOR_ASSERT_SAME_KEYS(A,B)					\
    typedef typename ::tensor::detail::assert::same_keys<A,B>::type	\
    BOOST_PP_CAT(TENSOR_ASSERT_SAME_KEYS_, __LINE__) 


#define TENSOR_ASSERT_UNIQUE_KEYS(A)					\
    typedef typename ::tensor::detail::assert::unique_keys<A>::type	\
    BOOST_PP_CAT(TENSOR_ASSERT_UNIQUE_KEYS_, __LINE__) 


#endif // TENSOR_ASSERT_HPP
