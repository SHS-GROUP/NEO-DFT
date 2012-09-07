#ifndef PHOENIX_UBLAS_HPP
#define PHOENIX_UBLAS_HPP

#include <boost/spirit/home/phoenix/core.hpp>
#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/spirit/home/phoenix/operator/self.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/mpl/bool.hpp>

namespace boost {
namespace phoenix {
namespace ublas {

    struct data_eval {
	template<class A>
	struct result {
	    typedef typename A::array_type& type;
	};
	template<class A>
	struct result<const A> {
	    typedef const typename A::array_type& type;
	};
	template<class A>
	typename result<A>::type
	operator()(A &a) const {
	    return a.data();
	}
    };
    const function<data_eval> data = function<data_eval>();

    struct trans_eval {
	template<class E>
	struct result {
	    typedef typename numeric::ublas::scalar_identity<
		typename E::value_type> scalar_identity;
	    typedef typename numeric::ublas::matrix_unary2_traits<
		const E, scalar_identity>::result_type type;
	};
	template<class E>
	typename result<E>::type
	operator()(const numeric::ublas::matrix_expression<E> &e) const {
	    return numeric::ublas::trans(e);
	}
    };

    const function<trans_eval> trans = function<trans_eval>();


    struct prod_eval {
	template<class E1, class E2>
	struct result {
	    typedef typename numeric::ublas::matrix_matrix_binary_traits<
		typename E1::value_type, E1, typename E2::value_type, E2
		>::result_type type;
	};
	template<class E1, class E2>
	typename result<E1,E2>::type
	operator()(const numeric::ublas::matrix_expression<E1> &e1,
		   const numeric::ublas::matrix_expression<E2> &e2) const {
	    return numeric::ublas::prod(e1,e2);
	}
    };

    const function<prod_eval> prod = function<prod_eval>();

}

template<class C, class E>
struct result_of_assign<numeric::ublas::noalias_proxy<C>,E> {
    typedef typename numeric::ublas::noalias_proxy<C>::closure_type& type;
};

namespace ublas {

    struct noalias_eval {
    	template<class C>
    	struct result {
    	    typedef numeric::ublas::noalias_proxy<C> type;
    	};
    	template<class C>
    	typename result<C>::type
    	operator()(C &lvalue) const {
    	    return numeric::ublas::noalias(lvalue);
    	}
    };
    const function<noalias_eval> noalias = function<noalias_eval>();

}
}
}

#endif // PHOENIX_UBLAS_HPP
