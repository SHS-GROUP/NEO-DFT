#ifndef BOOST_NUMERIC_UBLAS_ADAPTOR_HPP
#define BOOST_NUMERIC_UBLAS_ADAPTOR_HPP

#include "boost/numeric/ublas/storage_adaptors.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/mpl/if.hpp>
#include <boost/type_traits/is_const.hpp>
#include <boost/type_traits/remove_const.hpp>

namespace boost {
namespace numeric {
namespace ublas {
namespace detail {

    template<typename T>
    struct array_adaptor {
	typedef typename boost::remove_const<T>::type value_type;
	typedef typename
	boost::mpl::if_<boost::is_const<T>,
			ublas::const_array_ref_adaptor<value_type>,
			ublas::array_ref_adaptor<value_type> >::type type;
    };

}
}
}
}

namespace boost {
namespace numeric {
namespace ublas {

    template<typename T>
    struct vector_adaptor {
	typedef typename detail::array_adaptor<T>::type array_type;
	typedef ublas::vector<T, array_type> type;
    };

    template<typename T, class F>
    struct matrix_adaptor {
	typedef typename detail::array_adaptor<T>::type array_type;
	typedef ublas::matrix<T, F, array_type> type;
    };

    template<typename T, class L, class F>
    struct symmetric_matrix_adaptor {
	typedef typename detail::array_adaptor<T>::type array_type;
	typedef ublas::symmetric_matrix<T, L, F, array_type> type;
    };

    template<typename T>
    typename detail::vector_adaptor<T>::type
    make_vector(size_t m, T* data) {
	typedef typename detail::vector_adaptor<T>::type adaptor;
	typedef typename adaptor::array_type array_type;
	adaptor A(m, array_type(0, NULL));
	array_type a(m, data);
	A.data().swap(a);
	return A;
    };

    template<class L, typename T>
    typename detail::matrix_adaptor<T,L>::type
    make_matrix(size_t m, size_t n, T* data, const L& = L()) {
	typedef typename detail::matrix_adaptor<T,L>::type adaptor;
	typedef typename adaptor::array_type array_type;
	adaptor A(m, n, array_type(0, NULL));
	array_type a(m*n, data);
	A.data().swap(a);
	return A;
    };

    template<class L, class F, typename T>
    typename detail::symmetric_matrix_adaptor<T,L,F>::type
    make_matrix(size_t m, T* data, const L& = L(), const F& = F()) {
	typedef typename detail::symmetric_matrix_adaptor<T,L,F>::type adaptor;
	typedef typename adaptor::array_type array_type;
	adaptor A(m, array_type(0, NULL));
	array_type a((m*m+m)/2, data);
	A.data().swap(a);
	return A;
    };

}
}
}

#endif // BOOST_NUMERIC_UBLAS_ADAPTOR_HPP

