#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_FORWARD_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_FORWARD_HPP

#include <cstdlib>
#include <boost/numeric/ublas/storage.hpp>



namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

    using ublas::range;

    // forward declaration

    template<typename T>
    struct unbounded_array;

    template<typename T, class A = unbounded_array<T> >
    struct matrix;

    template<typename T, class L = ublas::column_major>
    struct matrix_adaptor;

    template<class M, class U>
    struct matrix_vector;

    template<class M>
    struct matrix_column;

    template<class M>
    struct matrix_row;

    template<class A>
    matrix_column<A> column(A &a, size_t index) {
	return matrix_column<A>(a, index);
    }

    template<class A>
    matrix_row<A> row(A &a, size_t index) {
	return matrix_row<A>(a, index);
    }

    template<class E>
    struct matrix_expression;

    template<class E>
    struct matrix_range;

    template<class E>
    struct vector_expression;

    template<class E>
    struct cublas_matrix_expression;

    template<class E>
    struct cublas_transpose;

}
}
}
}

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_FORWARD_HPP
