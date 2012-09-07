#ifndef TENSOR_FUNCTIONAL_HPP
#define TENSOR_FUNCTIONAL_HPP

#include "tensor/forward.hpp"
#include "tensor/expression.hpp"
#include "tensor/index.hpp"
#include "tensor/assert.hpp"
#include "tensor/lambda.hpp"

#include <boost/static_assert.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/utility/enable_if.hpp>

#include <boost/mpl/int.hpp>

#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/pair.hpp>
#include <boost/fusion/include/intrinsic.hpp>
#include <boost/fusion/include/algorithm.hpp>

#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/boost_array.hpp>

#include <stdexcept>

namespace tensor {
namespace detail {

    struct functional {


	template<size_t N, class A, class F, class I, class S>
	static typename boost::enable_if_c<(N == 0)>::type
	apply(fusion::vector<A&> v, const F &f, const I &indices, const S &map) {
	    lambda::apply(f, fusion::at_c<0>(v)[map]);
	}

	template<size_t N, class A, class B, class F, class I, class S>
	static typename boost::enable_if_c<(N == 0)>::type
	apply(fusion::vector<A&,B&> v, const F &f,
	      const I &indices, const S &map) {
	    lambda::apply(f, fusion::at_c<0>(v)[map], fusion::at_c<1>(v)[map]);
	}

	template<size_t N, class Tie, class F, class I, class S>
	static
	typename boost::enable_if_c<(N > 0)>::type
	apply(Tie tie, const F &f, const I &indices, const S &map) {
	    namespace fusion = boost::fusion;
	    BOOST_AUTO(const &kr, fusion::at_c<N-1>(indices));
	    typedef BOOST_TYPEOF(fusion::at_c<0>(kr)) K;
	    BOOST_AUTO(const &r, fusion::at_c<1>(kr));
	    BOOST_AUTO(it, r.begin());
	    while (it < r.end()) {
		BOOST_AUTO(const &i, *it);
		typedef fusion::pair<K, const BOOST_TYPEOF(i)&> pair;
		apply<N-1>(tie, f, indices, fusion::push_front(map, pair(i)));
		++it;
	    }
	}

	template<class A, class B>
	struct index_ranges_gen {
	    const A &a;
	    const B &b;
	    mutable index_range *ranges;
	    index_ranges_gen(const A &a, const B &b, index_range *ranges)
		: a(a), b(b), ranges(ranges) {}
	    template<class K>
	    void operator()(const K&) const {
		namespace fusion = boost::fusion;
		static const bool ka = fusion::result_of::has_key<A,K>::type::value;
		static const bool kb = fusion::result_of::has_key<B,K>::type::value;
		check<K>(boost::mpl::bool_<(ka && kb)>());
		(apply<K>(a, *ranges, boost::mpl::bool_<ka>()) ||
		 apply<K>(b, *ranges, boost::mpl::bool_<kb>()));
		++ranges;
	    }

	    template<class K>
	    void check(boost::mpl::true_) const {
		// BOOST_AUTO(const &ra, boost::fusion::at_key<K>(a));
		// BOOST_AUTO(const &rb, boost::fusion::at_key<K>(b));
		// if (ra.zero_based() != rb.zero_based()) {
		//     throw std::range_error("invalid range");
		// }
	    }
	    template<class K>
	    void check(boost::mpl::false_) const {}

	    template<class K, class A_>
	    static bool apply(A_ &a, index_range& r, boost::mpl::true_) {
		r = boost::fusion::at_key<K>(a);//.zero_based();
		return true;
	    }
	    template<class K, class A_>
	    static bool apply(A_ &a, index_range& r, boost::mpl::false_) {
		return false;
	    }
	};

	template<class A, class B, class K>
	static boost::array<index_range, boost::mpl::size<K>::value>
	index_ranges(const A &a, const B &b, const K &k) {
 	    static const size_t N = boost::mpl::size<K>::value;
 	    boost::array<index_range, N> ranges;
 	    index_ranges_gen<A,B> g(a, b, ranges.begin());
	    boost::fusion::for_each(k, g);
	    return ranges;
	}    

	template<class A, class B, class F, class S, class K>
	static void apply(A &a, B &b, const F &f, const S &map, const K &k) {
	    static const size_t N = boost::mpl::size<K>::value;
	    BOOST_AUTO(const &ranges, index_ranges(a.indices(), b.indices(), k));
	    BOOST_AUTO(const &indices, boost::fusion::zip(k, ranges));
	    apply<N>(fusion::vector_tie(a, b), f,
		     boost::fusion::as_vector(indices), map);
	}

	template<class A, class B, class F, class S>
	static void apply(A &a, B &b, const F &f, const S &map) {
	    namespace mpl = boost::mpl;
	    namespace fusion = boost::fusion;
	    namespace result_of = boost::fusion::result_of;

	    typedef typename A::indices_type I;
	    typedef result_of::has_key<S, mpl::_> key_in_S;
	    typedef result_of::has_key<I, mpl::_> key_in_A;

	    BOOST_AUTO(const &keys,
	    	       (fusion::remove_if<key_in_S>
	    		 (fusion::join
	    		  (a.keys(),
			   fusion::remove_if<key_in_A>(b.keys())))));
	    apply(a, b, f, map, keys);
	}

    };

} // namespace detail
} // namespace tensor


namespace tensor {

    template<class F, class A, class S>
    typename boost::result_of<const F(typename A::value_type, S)>::type
    reduce(const F &f, const expression<A> &a, const S &s) {
	typename boost::result_of<const F
	    (typename A::value_type, S)>::type state(s);
	using namespace lambda;
	detail::functional::apply(a(), (ref(state) = bind(f, (ref(state), _1))));
    }
	

    template<class A, class B, class F, class S>
    void contract(const F &f, const expression<A> &a, const expression<B> &b,
		  const S &map) {
	detail::functional::apply(a(), b(), f, boost::fusion::as_map(map));
    }

    template<class F, class A, class B>
    void transform(const F &f, expression<A> &a, const expression<B> &b) {
	TENSOR_ASSERT_SAME_RANK(A,B);
	TENSOR_ASSERT_SAME_KEYS(A,B);
	detail::functional::apply(a(), b(), f,
				  boost::fusion::map<>(), a().keys());
    }

    template<class A>
    void assign(expression<A> &a, const typename A::value_type &s) {
	transform((lambda::_1 = lambda::val(s)), a);
    }

    template<class A>
    void times_assign(expression<A> &a, const typename A::value_type &s) {
	transform((lambda::_1 *= lambda::val(s)), a);
    }



#define TENSOR_OPERATION(NAME, FUNC)						\
    template<class A, class B>							\
    void NAME(expression<A> &a, const expression<B> &b) {			\
	transform(FUNC, a, b);						\
    }

    TENSOR_OPERATION(plus_assign, (lambda::_1 += lambda::_2))
    TENSOR_OPERATION(assign, (lambda::_1 = lambda::_2));
}

#endif // TENSOR_FUNCTIONAL_HPP
