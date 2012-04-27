#ifndef BOOST_MPL_SWAP_HPP
#define BOOST_MPL_SWAP_HPP

#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_void.hpp>

#include <boost/mpl/advance.hpp>
#include <boost/mpl/distance.hpp>
#include <boost/mpl/erase.hpp>
#include <boost/mpl/insert.hpp>


namespace boost {
namespace mpl {
namespace detail {

    struct swap {

	template<template<class, class, class> class F, class S,
		 class I,
		 class A = typename mpl::next<
		     typename mpl::advance<
			 typename mpl::begin<S>::type, I>::type
		     >::type
	>
	struct apply_at {
	    typedef typename F<
		S, 
		typename mpl::advance<
		    typename mpl::begin<S>::type, I>::type,
		A
		>::type type;
	};

	template<class S, class I, class T>
	struct replace_at {
	    typedef typename apply_at<
		mpl::insert,
		typename apply_at<mpl::erase, S, I>::type,
		I,
		T
		>::type type;
	};

    }; // struct swap

} // namespace detail


    template<class S, class I, class J>
    struct swap {
	typedef typename mpl::distance<typename mpl::begin<S>::type, I>::type M;
	typedef typename mpl::distance<typename mpl::begin<S>::type, J>::type N;
	typedef typename detail::swap::replace_at<
	    typename detail::swap::replace_at<
		S, M, typename mpl::deref<J>::type>::type,
	    N, typename mpl::deref<I>::type
	    >::type type;
    };

} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_SWAP_HPP
