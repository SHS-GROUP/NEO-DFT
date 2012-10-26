#ifndef TENSOR_CONTRACT_HPP
#define TENSOR_CONTRACT_HPP

#include "tensor/expression.hpp"
#include "tensor/detail/contract.hpp"
#include <boost/utility/result_of.hpp>

namespace tensor {


    template<typename Alpha, class A, class B, typename Beta, class C>
    void contract(Alpha alpha, A a, B b, Beta beta, C c) {
	detail::contract_impl(alpha, a, b, beta, c);
    }

    template<class A, class B, typename U = int>
    struct contract_expression :
	expression<contract_expression<A,B,U> >
    {
	contract_expression(A a, B b, U alpha = U(1))
	    : a(a), b(b), alpha_(alpha) {}

	template<typename T>
	struct scaled {
	    typedef contract_expression<A, B, BOOST_TYPEOF_TPL(T()*U())> type;
	};

	template<typename T>
	typename scaled<T>::type
	operator*(T alpha) const {
	    return typename scaled<T>::type(a, b, alpha*alpha_);
	}
	template<typename T, class C>
	void evaluate(T alpha, T beta, expression<C> &c) const {
	    ::tensor::contract(alpha*alpha_, a, b, beta, c());
	} 
    private:
	A a; B b;
	U alpha_;
    };

    template<typename T, class A, class B, typename U>
    typename contract_expression<A,B,U>::template scaled<T>::type
    operator*(T alpha, contract_expression<A,B,U> e) {
	return e*alpha;
    }

#define TENSOR_CONTRACT_ASSIGN(NAME, ALPHA, BETA)				\
    template<class C, class A, class B, typename U>				\
    void									\
    NAME(expression<C> &c, const expression<contract_expression<A,B,U> > &e) {	\
	e().evaluate(ALPHA, BETA, c);						\
    }

    TENSOR_CONTRACT_ASSIGN(assign, 1, 0)
    TENSOR_CONTRACT_ASSIGN(plus_assign, 1, 1)
    TENSOR_CONTRACT_ASSIGN(minus_assign, -1, 1)

    template<class A, class B>
    contract_expression<A,B> contract(A a, B b) {
	return contract_expression<A,B>(a, b);
    }

}

#endif // TENSOR_CONTRACT_HPP
