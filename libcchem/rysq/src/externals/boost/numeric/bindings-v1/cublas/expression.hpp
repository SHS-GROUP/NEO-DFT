#ifndef BOOST_BINDINGS_CUBLAS_EXPRESSION_HPP
#define BOOST_BINDINGS_CUBLAS_EXPRESSION_HPP

#include <boost/numeric/bindings/cublas/forward.hpp>
#include <boost/numeric/bindings/cublas/cublas3.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

    template<class E>
    struct blas_transpose : matrix_expression<blas_transpose<E> > {
	typedef typename E::const_closure_type expression_closure_type;
	explicit blas_transpose(E &e) : e_(e) {}
	const E& operator()() const { return e_; } 
	size_t size1() const { return e_.size1(); }
	size_t size2() const { return e_.size2(); }
    private:
	E &e_;
    };

    template<class E>
    struct matrix_expression<blas_transpose<E> > {
	const blas_transpose<E>& operator()() const {
	    return *static_cast<const blas_transpose<E>*>(this);
	}
    };

    template<class E>
    blas_transpose<E> trans(E &e) {
	return blas_transpose<E>(e);
    }
    

    template<class E>
    struct blas_matrix_expression
	: matrix_expression<blas_matrix_expression<E> >
    {
	typedef E expression_type;
	const expression_type& operator()() const {
	    return *static_cast<const expression_type*>(this);
	}
	size_t size1() const { return operator()().size1(); }
	size_t size2() const { return operator()().size2(); }
	template<typename T, class C>
	void operator()(T beta, cublas::matrix_expression<C> &c) const {
	    operator()()(beta, c);
	}
    };

    template<class E>
    struct matrix_expression<blas_matrix_expression<E> > {
	typedef blas_matrix_expression<E> expression_type;
	const expression_type& operator()() const {
	    return *static_cast<const expression_type*>(this);
	}
	size_t size1() const { return operator()().size1(); }
	size_t size2() const { return operator()().size2(); }
    };


    template<class A, class B>
    void matrix_assign(cublas::matrix_expression<A> &a,
		       const cublas::matrix_expression<
		       blas_matrix_expression<B> > &b) {
	b()(typename A::value_type(0), a);
    }

    template<class A, class B>
    void matrix_plus_assign(cublas::matrix_expression<A> &a,
			    const cublas::matrix_expression<
			    blas_matrix_expression<B> > &b) {
	b()(typename A::value_type(1), a);
    }


    template<class A, class B>
    struct gemm_expression : blas_matrix_expression<gemm_expression<A,B> > {
	gemm_expression(const A &a, const B &b) : a_(a), b_(b) {}
	size_t size1() const { return a_.size1(); }
	size_t size2() const { return b_.size2(); }
	template<typename T, class C>
	void operator()(T beta, cublas::matrix_expression<C> &c) const {
	    gemm(T(1), a_, b_, beta, c);
	}
    private:
	const A &a_;
	const B &b_;
    };

    template<class A, class B>
    gemm_expression<A, B>
    gemm(const A &a, const B &b) {
      // std::cout << a.size1() << std::endl;
	return gemm_expression<A, B>(a, b);
    }

}
}
}
}

#endif // BOOST_BINDINGS_CUBLAS_EXPRESSION_HPP
