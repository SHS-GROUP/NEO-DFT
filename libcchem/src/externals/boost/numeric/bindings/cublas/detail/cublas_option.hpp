//
// Copyright (c) 2009 Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS_OPTION_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS_OPTION_HPP

#include <boost/numeric/bindings/tag.hpp>
#include <boost/numeric/bindings/cublas/cublas.h>

namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {
namespace detail {

template< cublasOperation_t op >
struct cublas_operation_t {
    static const cublasOperation_t value = op;
};

template< typename Tag >
struct cublas_option {};

template<>
struct cublas_option< tag::transpose >: cublas_operation_t< CUBLAS_OP_T > {};

template<>
struct cublas_option< tag::no_transpose >: cublas_operation_t< CUBLAS_OP_N > {};

template<>
struct cublas_option< tag::conjugate >: cublas_operation_t< CUBLAS_OP_C > {};

// template<>
// struct cublas_option< tag::upper >: mpl::char_< 'U' > {};

// template<>
// struct cublas_option< tag::lower >: mpl::char_< 'L' > {};

// template<>
// struct cublas_option< tag::unit >: mpl::char_< 'U' > {};

// template<>
// struct cublas_option< tag::non_unit >: mpl::char_< 'N' > {};

// template<>
// struct cublas_option< tag::left >: mpl::char_< 'L' > {};

// template<>
// struct cublas_option< tag::right >: mpl::char_< 'R' > {};

} // namespace detail
} // namespace cublas
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
