#ifndef BOOST_MPL_NEXT_PERMUTATION_HPP
#define BOOST_MPL_NEXT_PERMUTATION_HPP

#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/is_void.hpp>

#include "boost/mpl/swap.hpp"
#include "boost/mpl/identity.hpp"
#include <boost/mpl/vector.hpp>
#include <boost/mpl/reverse.hpp>
#include <boost/mpl/insert_range.hpp>
#include <boost/mpl/comparison.hpp>

namespace boost {
namespace mpl {
namespace detail {

    template<class S>
    struct next_permutation_index {

	template<class K, class I = void>
	struct index {

	    template<class I_>
	    struct get_iterator {
		typedef typename I_::iterator type;
	    };

	    typedef K previous;
	    typedef typename mpl::eval_if<
		boost::is_void<I>,
		mpl::identity<K>, get_iterator<I>
		>::type iterator;
	};

	template<class _1, class _2>
	struct test :
	    mpl::greater<
	    typename mpl::deref<_2>::type,
	    typename mpl::deref<typename _1::previous>::type
	    >::type {};

	typedef typename mpl::begin<S>::type begin;

	typedef typename mpl::iter_fold<
	    S, index<begin>,
	    mpl::if_<
	    	test<mpl::_1, mpl::_2>,
	    	index<mpl::_2>,
	    	index<mpl::_2, mpl::_1>
	    	>
	    >::type::iterator next_iterator;


	typedef mpl::not_<boost::is_same<next_iterator, begin> > status;
	typedef typename mpl::if_<
	    status,
	    typename mpl::prior<next_iterator>::type,
	    next_iterator
	    >::type iterator;
    };

} // namespace detail

    template<class S>
    struct next_permutation {
	typedef detail::next_permutation_index<S> I;
	BOOST_STATIC_CONSTANT(bool, value = I::status::value);

	typedef typename mpl::iterator_range<
	    typename I::next_iterator,
	    typename mpl::end<S>::type
	    >::type T;
	typedef detail::next_permutation_index<T> J;

	typedef typename mpl::swap<
	    S, typename I::iterator, typename J::iterator>::type U;

	typedef typename mpl::distance<
	    typename mpl::begin<S>::type,
	    typename I::next_iterator
	    >::type N;

	typedef typename mpl::iterator_range<
	    typename mpl::begin<U>::type,
	    typename mpl::advance<typename mpl::begin<U>::type, N>::type
	    >::type P;

	typedef typename mpl::iterator_range<
	    typename mpl::advance<typename mpl::begin<U>::type, N>::type,
	    typename mpl::end<U>::type
	    >::type R;

	typedef typename mpl::insert_range<
	    mpl::vector<>,
	    typename mpl::begin<mpl::vector<> >::type,
	    typename mpl::joint_view<
		P,
		typename mpl::reverse_copy<
		    R, mpl::back_inserter<mpl::vector<> >
		    >::type
		>
	    >::type type;
    };


} // namespace mpl
} // namespace boost

#endif // BOOST_MPL_NEXT_PERMUTATION_HPP
