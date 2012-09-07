#ifndef TENSOR_EXPRESSION_HPP
#define TENSOR_EXPRESSION_HPP

#include "tensor/config.hpp"
#include "tensor/generator.hpp"
#include "tensor/lambda.hpp"
#include "tensor/index.hpp"
#include "tensor/functional.hpp"
#include "tensor/assert.hpp"

#include <boost/mpl/contains.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/transformation.hpp>
#include <boost/fusion/include/mpl.hpp>


namespace tensor {

    template<class E>
    struct expression {
	operator E&() { return operator()(); }
	operator const E&() const { return operator()(); }
	E& operator()() { return *static_cast<E*>(this); }
	const E& operator()() const { return *static_cast<const E*>(this); }
    };

    template<class F, class T, size_t N = T::rank>
    struct unary_expression;

    template<class F, class A, class B, size_t N = A::rank>
    struct binary_expression;

}


namespace tensor {
namespace detail {


    template<class F, class A, size_t N = A::rank>
    struct unary_expression_generator {
	typedef unary_expression_generator self_type;
	static const size_t rank = N;

	template<class O>
	struct result : detail::fail_instantiate<O> {};

	template<class O, class I>
	struct result<O(indices<I>)> {
	    typedef typename A::template result<const A(indices<I>)> A_;
	    static const size_t rank = A_::rank;

	    template<class E>
	    struct expression {
		typedef unary_expression<F,E> type;
	    };

	    typedef typename boost::mpl::eval_if_c<
		rank,
		expression<typename A_::type>,
		boost::result_of<F(typename A_::type)>
		>::type type;
	};

	template<class G, class I>
	typename result<const self_type(indices<I>)>::type
	operator()(const G &g, const indices<I> &i) const {
	    typedef unary_expression<F,A> E;
	    return (static_cast<const E&>(g))(g, i);
	}
    };


    template<class F, class A, class B, size_t N = A::rank>
    struct binary_expression_generator {
	typedef binary_expression_generator self_type;
	static const size_t rank = N;

	template<class O>
	struct result : detail::fail_instantiate<O> {};

	template<class O, class I>
	struct result<O(indices<I>)> {
	    typedef typename A::template result<const A(indices<I>)> A_;
	    typedef typename B::template result<const B(indices<I>)> B_;

	    TENSOR_ASSERT_SAME_RANK(A_, B_);
	    static const size_t rank = A_::rank;

	    template<class A_, class B_>
	    struct expression {
		typedef binary_expression<F,A_,B_> type;
	    };

	    typedef typename boost::mpl::eval_if_c<
		rank,
		expression<typename A_::type, typename B_::type>,
		boost::result_of<F(typename A_::type, typename B_::type)>
		>::type type;
	};

	template<class G, class I>
	typename result<const self_type(indices<I>)>::type
	operator()(const G &g, const indices<I> &i) const {
	    typedef binary_expression<F,A,B> E;
	    return (static_cast<const E&>(g))(g, i);
	}
    };


} // namespace detail
} // namespace tensor


namespace tensor {

    template<class F, class A, size_t N>
    struct unary_expression :
	detail::const_generator<detail::unary_expression_generator<F,A> >,
	expression<unary_expression<F,A,N> >
    {
	typedef unary_expression self_type;
	typedef detail::const_generator<
	    detail::unary_expression_generator<F,A> > generator;

	static const size_t rank = N;
	typedef typename A::indices_type indices_type;
	typedef typename A::keys_type keys_type;

	typedef typename boost::result_of<
	    F(typename A::const_result_type)>::type const_result_type;
	
	unary_expression(const F &f, const expression<A> &a)
	    : f(f), a(a) {}

	indices_type indices() const { return indices_type(); }

	using detail::const_generator<
	    detail::unary_expression_generator<F,A> >::operator();

	template<class O>
	struct result : generator::template result<O> {};

	template<class I>
	typename boost::disable_if_c<
	    tensor::indices<I>::rank,
	    typename result<const self_type(tensor::indices<I>)>::type
	    >::type
	operator()(const generator&, const tensor::indices<I> &i) const {
	    return f(a(i));
	}

	template<class I>
	typename boost::enable_if_c<
	    tensor::indices<I>::rank,
	    typename result<const self_type(tensor::indices<I>)>::type
	    >::type
	operator()(const generator&, const tensor::indices<I> &i) const {
	    return unary_expression(f, a(i));
	}

    private:
	F f;
	A a;
    };


    template<class F, class A, class B, size_t N>
    struct binary_expression :
	detail::const_generator<detail::binary_expression_generator<F,A,B> >,
	expression<binary_expression<F,A,B,N> >
    {
	TENSOR_ASSERT_SAME_RANK(A,B);
	TENSOR_ASSERT_SAME_KEYS(A,B);

	typedef binary_expression self_type;
	typedef detail::const_generator<
	    detail::binary_expression_generator<F,A,B> > generator;

	static const size_t rank = N;
	typedef typename A::indices_type indices_type;
	typedef typename A::keys_type keys_type;

	// typedef typename boost::result_of<
	//     F(typename A::const_result_type)>::type const_result_type;
	
	binary_expression(const F &f,
			  const expression<A> &a,
			  const expression<B> &b)
	    : f(f), a(a), b(b) {}

	indices_type indices() const { return indices_type(); }

	using generator::operator();

	template<class O>
	struct result : generator::template result<O> {};

	template<class I>
	typename result<const self_type(tensor::indices<I>)>::type
	operator()(const generator&, const tensor::indices<I> &i) const {
	    typedef typename result<const self_type(tensor::indices<I>)>::type R;
	    return evaluate<R>(i);
	}

    private:
	F f;
	A a;
	B b;

    private:
	template<class R, class I>
	typename boost::disable_if_c<tensor::indices<I>::rank, R>::type
	evaluate(const tensor::indices<I> &i) const {
	    return f(a(i), b(i));
	}

	template<class R, class I>
	typename boost::enable_if_c<tensor::indices<I>::rank, R>::type
	evaluate(const tensor::indices<I> &i) const {
	    return R(f, a(i), b(i));
	}
    };


    template<class F, class A>
    unary_expression<F,A>
    apply(const F &f, const expression<A> &a) {
	return unary_expression<F,A>(f, a);
    }

    template<class F, class A, class B>
    binary_expression<F,A,B>
    apply(const F &f, const expression<A> &a, const expression<B> &b) {
	return binary_expression<F,A,B>(f, a, b);
    }


} // namespace tensor

#endif // TENSOR_EXPRESSION_HPP

