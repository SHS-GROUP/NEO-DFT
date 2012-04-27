#ifndef TENSOR_OPERATOR_HPP
#define TENSOR_OPERATOR_HPP

#include "tensor/expression.hpp"
#include "tensor/lambda.hpp"

#include <boost/utility/enable_if.hpp>
#include <boost/mpl/logical.hpp>

#include <boost/fusion/swap.hpp>
#include <boost/fusion/iterator.hpp>
#include <boost/fusion/mpl.hpp>

#include <boost/type_traits/is_base_of.hpp>
#include <boost/type_traits/remove_const.hpp>

namespace tensor {
namespace traits {

    template<class E>
    struct remove_cr {
	typedef typename boost::remove_const<
	    typename boost::remove_reference<E>::type
	    >::type T;
    };

    template<class E>
    struct is_expression {
	typedef typename boost::is_base_of<expression<E>, E>::type type;
	BOOST_STATIC_CONSTANT(bool, value = type::value);
    };

    template<class E>
    struct is_expression<expression<E> > : boost::mpl::true_ {};

    template<class E>
    struct is_expression<const E> : is_expression<E> {};


    template<class E>
    struct expression_ {
	BOOST_MPL_ASSERT((is_expression<E>));
	typedef E type;
    };

    template<class E>
    struct expression_<expression<E> > {
	typedef E type;
    };

    template<class E>
    struct expression_<const E> : expression_<E> {};


}
}

namespace tensor {

    template<class A, class B, class = void>
    struct result_of_plus;

    template<class A, class B>
    struct result_of_plus<
	A, B, typename boost::enable_if<
		  boost::mpl::and_<
		      traits::is_expression<typename boost::remove_reference<A>::type>,
		      traits::is_expression<typename boost::remove_reference<B>::type>
		      >
		  >::type
	>
    {
	typedef typename boost::remove_reference<A>::type A_;
	typedef typename boost::remove_reference<B>::type B_;
	typedef binary_expression<
	    BOOST_TYPEOF(lambda::_1 + lambda::_2),
	    typename traits::expression_<A_>::type,
	    typename traits::expression_<B_>::type
	    > type;
    };

    template<class A, class B, class = void>
    struct result_of_minus;

    template<class A, class B>
    struct result_of_minus<
	A, B, typename boost::enable_if<
		  boost::mpl::and_<
		      traits::is_expression<typename boost::remove_reference<A>::type>,
		      traits::is_expression<typename boost::remove_reference<B>::type>
		      >
		  >::type
	>
    {
	typedef typename boost::remove_reference<A>::type A_;
	typedef typename boost::remove_reference<B>::type B_;
	typedef binary_expression<
	    BOOST_TYPEOF(lambda::_1 - lambda::_2),
	    typename traits::expression_<A_>::type,
	    typename traits::expression_<B_>::type
	    > type;
    };


}


namespace tensor {
namespace Operator {

    template<class F>
    struct operator_expression;

    template<class F, class G>
    struct operator_sum {
	typedef operator_sum self_type;
	operator_sum(const F &f, const G &g) : f(f), g(g) {}

	template<class O>
	struct result;

	template<class E>
	struct result<self_type(E)> {
	    typedef typename result_of_plus<
		typename F::template result<F(E)>::type,
		typename G::template result<G(E)>::type
		>::type type;
	};

	template<class E>
	typename result<self_type(const expression<E>&)>::type
	operator()(const expression<E> &e) const {
	    return (f(e) + g(e));
	}
    private:
	F f; G g;
    };

    template<class F, class G>
    struct operator_difference {
	typedef operator_difference self_type;
	operator_difference(const F &f, const G &g) : f(f), g(g) {}

	template<class O>
	struct result;

	template<class E>
	struct result<self_type(E)> {
	    typedef typename result_of_minus<
		typename F::template result<F(E)>::type,
		typename G::template result<G(E)>::type
		>::type type;
	};

	template<class E>
	typename result<self_type(const expression<E>&)>::type
	operator()(const expression<E> &e) const {
	    return (f(e) - g(e));
	}
    private:
	F f; G g;
    };


    template<class F, class G>
    struct operator_product {
	typedef operator_product self_type;
	operator_product(const F &f, const G &g) : f(f), g(g) {}

	template<class O>
	struct result;

	template<class E>
	struct result<self_type(E)> {
	    typedef typename F::template result<
		F(typename G::template result<G(E)>::type)
		>::type type;
	};

	template<class E>
	typename result<self_type(const expression<E>&)>::type
	operator()(const expression<E> &e) const {
	    return f(g(e));
	}
    private:
	F f; G g;
    };


    template<class F>
    struct operator_expression {
	typedef operator_expression self_type;
	explicit operator_expression(const F &f = F()) : f(f) {}

	template<class O>
	struct result;// : F::template result<O> {};

	template<class E>
	struct result<self_type(const expression<E>&)> {
	    typedef typename F::template result<
		F(const expression<E>&)>::type type;
	};

	template<class G>
	struct result<self_type(const operator_expression<G>&)> {
	    typedef operator_expression<operator_product<F,G> > type;
	};

     	template<class E>
	typename result<self_type(const expression<E>&)>::type
	operator()(const expression<E> &e) const {
	    return f(e);
	}

     	template<class G>
	typename result<self_type(const operator_expression<G>&)>::type
	operator()(const operator_expression<G> &g) const {
	    typedef operator_product<F,G> P;
	    return operator_expression<P>(P(f,g));
	}

	operator const F&() const { return f; }
    private:
	F f;
    };

    template<class F, class G>
    operator_expression<operator_difference<F,G> >
    operator-(const operator_expression<F> &f,
	      const operator_expression<G> &g) {
	typedef operator_difference<F,G> S;
	return operator_expression<S>(S(f,g));
    }

    template<class F, class G>
    operator_expression<operator_sum<F,G> >
    operator+(const operator_expression<F> &f,
	      const operator_expression<G> &g) {
	typedef operator_sum<F,G> S;
	return operator_expression<S>(S(f,g));
    }

    struct identity_operator {
	template<class F>
	struct result;

	template<class F, class E>
	struct result<F(E)> {
	    typedef E type;
	};

	template<class E>
	typename result<identity_operator(const expression<E>&)>::type
	operator()(const expression<E> &e) const { return e; }
    };
    typedef operator_expression<identity_operator> Identity;
    const Identity I = Identity();
    

    template<class I, class J>
    struct permutation_operator {

	template<class S>
	struct permute {

	    typedef typename boost::fusion::result_of::as_vector<S>::type V;
	    BOOST_MPL_ASSERT((boost::mpl::contains<V,I>));
	    BOOST_MPL_ASSERT((boost::mpl::contains<V,J>));
	    
	    typedef typename boost::fusion::result_of::find<V,I>::type M;
	    typedef typename boost::fusion::result_of::find<V,J>::type N;

	    BOOST_MPL_ASSERT((boost::mpl::not_<boost::is_same<M,N> >));

	    typedef boost::mpl::bool_<
		(boost::fusion::result_of::distance<M, N>::value < 0)
		    > X;
	    // typedef typename boost::mpl::print<X>::type _;

	    typedef typename boost::fusion::result_of::as_vector<
		typename boost::mpl::eval_if<
		    X,
		    boost::fusion::result_of::swap<V, M, N>,
		    boost::mpl::identity<V>
		    >::type
		>::type type;

	    static type evaluate(const S &s) {
		return evaluate(boost::fusion::as_vector(s), X());
	    }
	private:
	    template<class V>
	    static type evaluate(const V &s, boost::mpl::false_) {
		return s;
	    }
	    template<class V>
	    static type evaluate(const V &s, boost::mpl::true_) {
		namespace fusion = boost::fusion;
		return (fusion::as_vector
			(fusion::swap
			 (s, fusion::find<I>(s), fusion::find<J>(s)))); 
	    }
	};

	template<class F>
	struct result;

	template<class F, class E>
	struct result<F(E)> {
	    typedef typename E::template result<
		const E(indices<typename permute<typename E::indices_type>::type>)
		>::type type;
	};

	template<class F, class E>
	struct result<F(const expression<E>&)> : result<F(const E)> {};

	template<class E>
	typename result<permutation_operator(const expression<E>&)>::type
	operator()(const expression<E> &e) const {
	    typedef permute<typename E::indices_type> P;
	    return e()(indices<typename P::type>(P::evaluate(e().indices())));
	}
    };


    struct permutation_generator {

	template<class I, class J>
	operator_expression<permutation_operator<I,J> >
	operator()(const I& i, const J& j) const {
	    return operator_expression<permutation_operator<I,J> >();
	}

    };
    const permutation_generator P = permutation_generator();


}
}

#endif // TENSOR_OPERATOR_HPP
