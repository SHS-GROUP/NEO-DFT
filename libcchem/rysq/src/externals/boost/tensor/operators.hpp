#ifndef TENSOR_OPERATORS_HPP
#define TENSOR_OPERATORS_HPP

#include "tensor/expression.hpp"
#include "tensor/lambda.hpp"
#include "tensor/product.hpp"

namespace tensor {


    template<class A>
    unary_expression<BOOST_TYPEOF_TPL(-lambda::_1), A>
    operator-(const expression<A> &a) {
	return apply((-lambda::_1), a);
    }


    namespace detail {

	template<typename T, class A,
		 class = typename boost::enable_if<
		     boost::is_convertible<T, typename A::value_type> >::type
		 >
	struct enable_if_convertible;

	template<typename T, class A, class = enable_if_convertible<T,A> >
	struct unary_multiply_expression {
	    typedef unary_expression<
		BOOST_TYPEOF(lambda::val(T())*lambda::_1), A> type;
	};

	template<typename T, class A, class = enable_if_convertible<T,A> >
	struct unary_divide_expression {
	    typedef unary_expression<
		BOOST_TYPEOF(lambda::val(T())/lambda::_1), A> type;
	};

    }


#define TENSOR_UNARY_OPERATOR(op, expr, A_, B_)		\
    template<typename T, class A>			\
    typename expr<T,A>::type				\
    operator op (A_, B_) {				\
	return apply((lambda::val(t)*lambda::_1), a);	\
    }

    TENSOR_UNARY_OPERATOR(*, detail::unary_multiply_expression,
			  const T &t, const expression<A> &a)

    TENSOR_UNARY_OPERATOR(*, detail::unary_multiply_expression,
			  const expression<A> &a, const T &t)

    TENSOR_UNARY_OPERATOR(/, detail::unary_divide_expression,
			  const T &t, const expression<A> &a)

    TENSOR_UNARY_OPERATOR(/, detail::unary_divide_expression,
			  const expression<A> &a, const T &t)


#undef TENSOR_UNARY_OPERATOR

}

namespace tensor {

    template<class A, class B>
    binary_expression<BOOST_TYPEOF(lambda::_1 + lambda::_2), A, B>
    operator+(const expression<A> &a, const expression<B> &b) {
	return apply((lambda::_1 + lambda::_2), a, b);
    }

    template<class A, class B>
    binary_expression<BOOST_TYPEOF(lambda::_1 - lambda::_2), A, B>
    operator-(const expression<A> &a, const expression<B> &b) {
	return apply((lambda::_1 - lambda::_2), a, b);
    }

}

#endif // TENSOR_OPERATORS_HPP
