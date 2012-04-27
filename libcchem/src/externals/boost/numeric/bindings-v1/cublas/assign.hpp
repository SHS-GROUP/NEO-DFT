#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_ASSIGN_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_ASSIGN_HPP

#include <boost/numeric/bindings/cublas/traits.hpp>
#include <boost/numeric/bindings/cublas/forward.hpp>
#include <boost/numeric/bindings/cublas/cublas.hpp>
#include <boost/numeric/bindings/cublas/cublas1.hpp>
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
	bool value = (traits::vector_size(a) == traits::vector_size(b));
	return value;
    }

    template<class A, class B>
    bool check_matrix_size(const A &a, const B &b) {
	bool value = ((traits::matrix_size1(a) == traits::matrix_size1(b)) && 
		      (traits::matrix_size2(a) == traits::matrix_size2(b)));
	return value;
    }

    template<class F, class X, class Y>
    void vector_assign(device_to_host, X &x, const Y &y) {
	check_vector_size(x, y);
	if (!traits::vector_size(x)) return;
	cublasGetVector(traits::vector_size(x),
			traits::vector_storage(y), traits::vector_stride(y),
			traits::vector_storage(x), traits::vector_stride(x));
    }

    template<class F, class X, class Y>
    void vector_assign(host_to_device, X &x, const Y &y) {
	check_vector_size(x, y);
	if (!traits::vector_size(x)) return;
	cublasSetVector(traits::vector_size(x),
			traits::vector_storage(y), traits::vector_stride(y),
			traits::vector_storage(x), traits::vector_stride(x));
    }

    template<class A, class B>
    void matrix_assign(device_to_host,
		       A &a, const cublas::matrix_expression<B> &b) {
	BOOST_STATIC_ASSERT((is_same<typename A::orientation_category,
			     typename B::orientation_category>::value));
	BOOST_STATIC_ASSERT((is_same<typename A::value_type,
			     typename B::value_type>::value));
	typedef typename A::value_type value_type;
	assert(check_matrix_size(a,b));
	const size_t m = traits::matrix_size1(a);
	const size_t n = traits::matrix_size2(a);
	if (!(m && n)) return;


	size_t ma = traits::leading_dimension(a);
	size_t mb = traits::leading_dimension(b);
	if (m == ma && m == mb) {
	    cudaMemcpy(traits::matrix_storage(a), traits::matrix_storage(b),
		       m*n*sizeof(value_type), cudaMemcpyDeviceToHost);
	}
	else {
	    cublasGetMatrix(m, n, sizeof(value_type),
			    traits::matrix_storage(b), mb,
			    traits::matrix_storage(a), ma);
	}
    }

    template<class A, class B>
    void matrix_assign(host_to_device,
		       cublas::matrix_expression<A> &a, const B &b) {
	BOOST_STATIC_ASSERT((is_same<typename A::orientation_category,
			     typename B::orientation_category>::value));
	BOOST_STATIC_ASSERT((is_same<typename A::value_type,
			     typename B::value_type>::value));
	typedef typename A::value_type value_type;
	assert(check_matrix_size(a,b));
	const size_t m = traits::matrix_size1(a);
	const size_t n = traits::matrix_size2(a);
	if (!(m && n)) return;

	size_t ma = traits::leading_dimension(a);
	size_t mb = traits::leading_dimension(b);
	if (m == ma && m == mb) {
	    cudaMemcpy(traits::matrix_storage(a), traits::matrix_storage(b),
		       m*n*sizeof(value_type), cudaMemcpyHostToDevice);
	}
	else {
	    cublasSetMatrix(m, n, sizeof(value_type),
			    traits::matrix_storage(b), mb,
			    traits::matrix_storage(a), ma);
	}
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
	detail::matrix_assign(detail::device_to_host(), a(), b);
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
	for (int j = 0; j < traits::matrix_size2(a); ++j) {
	    cublas::column(a(), j) = cublas::column(b(), j);
	}
    }

    template<typename T, class A>
    void matrix_clear(cublas::matrix<T,A> &a) {
	cudaMemset(a.data().begin(), 0, a.data().size()*sizeof(T));
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
