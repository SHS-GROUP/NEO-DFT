#ifndef BOOST_PP_IS_ITERATING

#ifndef TENSOR_PERMUTATION_HPP
#define TENSOR_PERMUTATION_HPP

#include <boost/preprocessor/repetition.hpp>
#include <boost/preprocessor/iteration/iterate.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/preprocessor/cat.hpp>

#include "tensor/config.hpp"
#include "tensor/index.hpp"
#include "tensor/expression.hpp"

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/next_permutation.hpp>
#include <boost/mpl/for_each.hpp>

#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/io.hpp>

namespace tensor {
namespace detail {

    template<class S, bool = true>
    struct permutation_vector {
	typedef boost::mpl::next_permutation<S> N;
	typedef typename boost::mpl::push_back<
	    // boost::mpl::vector<
	    // 	typename boost::mpl::next_permutation<
	    // 	    typename N::type>::type
	    // 	    >,
	    typename permutation_vector<typename N::type, N::value>::type,
	    typename N::type
	    >::type type;
    };

    template<class S>
    struct permutation_vector<S, false> {
	//typedef typename boost::mpl::print<S>::type _;
	typedef boost::mpl::vector<> type;
    };


    template<class A, class S, size_t N = boost::mpl::size<S>::value>
    struct permutation_expression {

	typedef typename permutation_vector<
	    typename boost::mpl::copy<
		boost::mpl::range_c<int,0,N>,
		boost::mpl::back_inserter<boost::mpl::vector<> >
		>::type
	    >::type V;

	struct expand {
	    template<class F>
	    struct result;

	    template<class E, class P>
	    struct result<expand(const E&, P)> {
		typedef binary_expression<
		    BOOST_TYPEOF(lambda::_1 + lambda::_2), E, A> type; 
	    };

	    template<class E, class P>
	    typename result<expand(const E&, P)>::type
	    operator()(const E &e, const P &p) const;
	};

	typedef typename boost::fusion::result_of::accumulate<
	    typename boost::fusion::result_of::pop_front<V>::type,
	    A, expand
	    >::type type;

	template<class T>
	void operator ()(const T &t)const {
	    //typename boost::mpl::print<T>::type _;
	    std::cout << boost::mpl::at_c<T,0>::type::value
		      << boost::mpl::at_c<T,1>::type::value
		      << std::endl;
	}
	static type generate(const A &a) {
	    boost::mpl::for_each<V>(permutation_expression());
	}
    };


} // namespace detail
} // namespace tensor


namespace tensor {

    template<class S>
    struct permutation_generator {
	S indices;
	explicit
	permutation_generator(const S &indices) : indices(indices) {}
	template<class A>
	typename detail::permutation_expression<A,S>::type
	operator()(const A &a) const {
	    return detail::permutation_expression<A,S>::generate(a);
	}
    };

}


// generate specializations
#define BOOST_PP_ITERATION_LIMITS (1, TENSOR_MAX_RANK)
#define BOOST_PP_FILENAME_1       "tensor/permutation.hpp" // this file
#include BOOST_PP_ITERATE()

#endif // TENSOR_PERMUTATION_HPP

#else // BOOST_PP_IS_ITERATING

#define N BOOST_PP_ITERATION()

namespace tensor {

    template<BOOST_PP_ENUM_PARAMS(N, class I)>
    permutation_generator<boost::fusion::vector<BOOST_PP_ENUM_PARAMS(N, I)> >
    permutation(BOOST_PP_ENUM_BINARY_PARAMS(N, const I, &i)) {
	typedef boost::fusion::vector<BOOST_PP_ENUM_PARAMS(N, I)> S;
	return permutation_generator<S>(S(BOOST_PP_ENUM_PARAMS(N, i)));
    }

}

#undef N

#endif // BOOST_PP_IS_ITERATING
