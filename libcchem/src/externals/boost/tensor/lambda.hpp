#ifndef TENSOR_LAMBDA_HPP
#define TENSOR_LAMBDA_HPP

#include <boost/spirit/include/phoenix_core.hpp>
#include <boost/spirit/home/phoenix/core/value.hpp>
#include <boost/spirit/include/phoenix_operator.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/scope/lambda.hpp>

namespace tensor {
    namespace lambda {

	using namespace boost::phoenix;
	using namespace boost::phoenix::arg_names;


	template<class F, class A>
	void apply(const F &f, A &a) {
	    f(a);
	}

	template<class F, class A, class B>
	void apply(const F &f, A &a, const B &b) {
	    f(a, b);
	}

	template<class F, class A, class B>
	void apply(const F &f, const A &a, const B &b) {
	    f(a, b);
	}

    }
}

#endif // TENSOR_LAMBDA_HPP
