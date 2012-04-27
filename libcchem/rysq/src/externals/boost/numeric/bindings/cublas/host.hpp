#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_HOST_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_HOST_HPP

#include <boost/numeric/bindings/cublas/forward.hpp>
#include <boost/numeric/ublas/matrix.hpp>


namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

    template<class C>
    struct host_matrix
	: ublas::matrix<typename C::value_type, ublas::column_major>
    {
	explicit host_matrix(const cublas::matrix_expression<C> &m) {
	    this->resize(m().size1(), m().size2(), false);
	    matrix_assign(*this, m);
	}
    };

    template<typename C>
    host_matrix<C> host(const cublas::matrix_expression<C> &m) {
	return host_matrix<C>(m);
    }


    template<class U>
    struct host_reference {
	typedef U matrix_type;
	explicit host_reference(matrix_type &m) : data_(m) {}
	template<class C>
	matrix_type& operator=(const cublas::matrix_expression<C> &m) {
	    //data_.clear();
	    matrix_assign(data_, m);
	    return data_;
	}
    private:
	matrix_type &data_;
    };


    template<class E>
    host_reference<E> host(ublas::matrix_expression<E> &m) {
	return host_reference<E>(m());
    }




}
}
}
}

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_HOST_HPP
