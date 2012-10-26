#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_VECTOR_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_VECTOR_HPP

#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/vector_proxy.hpp>

#include <boost/numeric/bindings/cublas/forward.hpp>
#include <boost/numeric/bindings/cublas/assign.hpp>
#include <boost/numeric/bindings/cublas/storage.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

    template<class E>
    struct vector_expression : protected ublas::vector_expression<E> {
	// friend class cublas::vector_expression<E>;
	typedef cublas::vector_expression<E> base;
	typedef typename base::expression_type expression_type;
	expression_type& operator()() {
	    return base::operator()(); }
	const expression_type& operator()() const { return base::operator()(); }
    };

}
}
}
}

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_VECTOR_HPP
