#ifndef STORAGE_HPP
#define STORAGE_HPP

#include "externals/boost/numeric/ublas/storage_adaptors.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/remove_const.hpp>

namespace bindings {
namespace detail {

    namespace ublas = boost::numeric::ublas;


    template<typename T>
    struct array_adaptor {
	typedef typename boost::remove_const<T>::type value_type;
	typedef typename
	boost::mpl::if_<boost::is_const<T>,
			ublas::const_array_ref_adaptor<value_type>,
			ublas::array_ref_adaptor<value_type> >::type type;
    };

    template<typename T>
    struct vector {
	typedef typename array_adaptor<T>::type array_type;
	typedef ublas::vector<T, array_type> type;
    };

    template<typename T>
    struct fortran_matrix {
	typedef typename array_adaptor<T>::type array_type;
	typedef ublas::matrix<T, ublas::column_major, array_type> type;
    };

    template<typename T>
    struct symmetric_fortran_matrix {
	typedef typename array_adaptor<T>::type array_type;
	typedef ublas::symmetric_matrix<T, ublas::upper,
					ublas::column_major, array_type> type;
    };

}
}

namespace bindings {

    template<typename T>
    struct vector : detail::vector<T>::type {
	typedef detail::vector<T> adapter;
	typedef typename adapter::type ublas_type;
	typedef typename ublas_type::array_type array_type;
	vector(size_t m, T* data)
	    : ublas_type(m, array_type(0, NULL))
	{
	    array_type array(m, data);
	    this->data().swap(array);
	}
    };


    template<typename T>
    struct fortran_matrix
	:  detail::fortran_matrix<T>::type {
	typedef detail::fortran_matrix<T> adapter;
	typedef typename adapter::type ublas_type;
	typedef typename ublas_type::array_type array_type;
	fortran_matrix(size_t m, size_t n, T* data)
	    : ublas_type(m, n, array_type(0, NULL))
	{
	    array_type array(m*n, data);
	    this->data().swap(array);
	}
	boost::numeric::ublas::matrix_range<ublas_type>
	range(size_t m, size_t n) {
	    namespace ublas = boost::numeric::ublas;
	    ublas::range r1(0,m), r2(0,n);
	    return ublas::matrix_range<ublas_type>(*this, r1, r2);
	}
    };

    template<typename T>
    struct symmetric_fortran_matrix
	:  detail::symmetric_fortran_matrix<T>::type {
	typedef detail::symmetric_fortran_matrix<T> adapter;
	typedef typename adapter::type ublas_type;
	typedef typename ublas_type::array_type array_type;
	symmetric_fortran_matrix(size_t m, T* data)
	    : ublas_type(m, array_type(0, NULL))
	{
	    array_type array((m*m + m)/2, data); 
	    this->data().swap(array);
	}
    };

    namespace detail {
	inline size_t triangular(size_t n) {
	    return (n*n + n)/2;
	}
    }

    template<class S, class L, typename T>
    boost::numeric::ublas::symmetric_matrix<
	T, S, L, boost::numeric::ublas::array_ref_adaptor<T> >
    make_matrix(size_t m, T *data, S = S(), L = L()) {
	namespace ublas = boost::numeric::ublas;
	typedef ublas::array_ref_adaptor<T> A;
	return ublas::symmetric_matrix<T, S, L, A>(m, A(detail::triangular(m), data));
    }


}

#endif // STORAGE_HPP

