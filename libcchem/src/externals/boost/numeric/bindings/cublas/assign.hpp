#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_ASSIGN_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_ASSIGN_HPP

#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/begin.hpp>

#include <boost/numeric/bindings/cublas/vector.hpp>
#include <boost/numeric/bindings/cublas/matrix.hpp>
#include <boost/numeric/bindings/cublas/level1.hpp>
#include <cuda_runtime.h>
#include <cublas.h>

#include <boost/numeric/ublas/vector_expression.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>

#include <boost/static_assert.hpp>


namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {


namespace detail {

    struct device_to_host {};
    struct host_to_device {};
    struct device_to_device {};

    template<class A, class B>
    bool check_vector_size(const A &a, const B &b) {
	bool value = (bindings::size(a) == bindings::size(b));
	return value;
    }

    template<class A, class B>
    bool check_matrix_size(const A &a, const B &b) {
	bool value = ((bindings::size_row(a) == bindings::size_row(b)) && 
		      (bindings::size_column(a) == bindings::size_column(b)));
	return value;
    }

    template<class F, class X, class Y>
    void vector_assign(device_to_host, X &x, const Y &y) {
	check_vector_size(x, y);
	if (!bindings::size(x)) return;
	cublasGetVector(bindings::size(x),
			bindings::begin_value(y), bindings::stride(y),
			bindings::begin_value(x), bindings::stride(x));
	BOOST_NUMERIC_BINDINGS_CUBLAS_CHECK_STATUS( );
    }

    template<class F, class X, class Y>
    void vector_assign(host_to_device, X &x, const Y &y) {
	check_vector_size(x, y);
	if (!bindings::size(x)) return;
	cublasSetVector(bindings::size(x),
			bindings::begin_value(y), bindings::stride(y),
			bindings::begin_value(x), bindings::stride(x));
	BOOST_NUMERIC_BINDINGS_CUBLAS_CHECK_STATUS( );
    }

    template<class A, class B>
    void matrix_assign(device_to_host, A &a, const B &b) {
	BOOST_STATIC_ASSERT((is_same<typename A::orientation_category,
			     typename B::orientation_category>::value));
	BOOST_STATIC_ASSERT((is_same<typename A::value_type,
			     typename B::value_type>::value));
	typedef typename A::value_type value_type;
	assert(check_matrix_size(a,b));
	const size_t m = bindings::size_row(b);
	const size_t n = bindings::size_column(b);
	if (!(m && n)) return;

	cublasGetMatrix(m, n, sizeof(value_type),
			bindings::begin_value(b), bindings::stride_major(b),
			bindings::begin_value(a), bindings::stride_major(a));
	BOOST_NUMERIC_BINDINGS_CUBLAS_CHECK_STATUS( );
    }

    template<class A, class B>
    void matrix_assign(host_to_device, A &a, const B &b) {
	BOOST_STATIC_ASSERT((is_same<typename A::orientation_category,
			     typename B::orientation_category>::value));
	BOOST_STATIC_ASSERT((is_same<typename A::value_type,
			     typename B::value_type>::value));
	typedef typename A::value_type value_type;
	assert(check_matrix_size(a,b));
	const size_t m = bindings::size_row(a);
	const size_t n = bindings::size_column(a);
	if (!(m && n)) return;

	cublasSetMatrix(m, n, sizeof(value_type),
			bindings::begin_value(b), bindings::stride_major(b),
			bindings::begin_value(a), bindings::stride_major(a));
	BOOST_NUMERIC_BINDINGS_CUBLAS_CHECK_STATUS( );
    }

}


    template<class X, class Y>
    void vector_assign(ublas::vector_expression<X> &x,
		       const cublas::vector_expression<Y> &y) {
	detail::vector_assign(detail::host_to_device(), x(), y);
    }

    template<class X, class Y>
    void vector_assign(cublas::vector_expression<X> &x,
		       const ublas::vector_expression<Y> &y) {
	detail::vector_assign(detail::device_to_host(), x, y());
    }

    template<class X, class Y>
    void vector_assign(cublas::vector_expression<X> &x,
		       const cublas::vector_expression<Y> &y) {
	cublas::copy(y, x);
    }


    template<class A, class B>
    void matrix_assign(ublas::matrix_expression<A> &a,
		       const cublas::matrix_expression<B> &b) {
	detail::matrix_assign(detail::device_to_host(), a(), b());
    }

    template<class A, class B>
    void matrix_assign(cublas::matrix_expression<A> &a,
		       const ublas::matrix_expression<B> &b) {
	detail::matrix_assign(detail::host_to_device(), a(), b());
    }

    template<class A, class B>
    void matrix_assign(cublas::matrix_expression<A> &a,
		       const cublas::matrix_expression<B> &b) {
	detail::check_matrix_size(a,b);
	for (int j = 0; j < bindings::size_column(a); ++j) {
	    cublas::column(a(), j) = cublas::column(b(), j);
	}
    }

    template<typename T, class L>
    void matrix_clear(const cublas::matrix_adaptor<T,L> &a);// {
    // 	cudaMemset(a.data().begin(), 0, a.data().size()*sizeof(T));
    // 	BOOST_NUMERIC_BINDINGS_CUBLAS_CHECK_STATUS( );
    // }

    template<typename T, class A>
    void matrix_clear(cublas::matrix<T,A> &a) {
	cudaMemset(a.data().begin(), 0, a.data().size()*sizeof(T));
	BOOST_NUMERIC_BINDINGS_CUBLAS_CHECK_STATUS( );
    }

    // template<class A, class B>
    // void  matrix_plus_assign(cublas::matrix_expression<A> &a,
    // 			     const cublas::matrix_expression<B> &b) {
    // }

}
}
}
}

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_ASSIGN_HPP
