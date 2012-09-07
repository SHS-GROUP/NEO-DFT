//
// Copyright (c) 2002--2010
// Toon Knapen, Karl Meerbergen, Kresimir Fresl,
// Thomas Klimpel and Rutger ter Borg
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//
// THIS FILE IS AUTOMATICALLY GENERATED
// PLEASE DO NOT EDIT!
//

#ifndef BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_GETRI_HPP
#define BOOST_NUMERIC_BINDINGS_LAPACK_COMPUTATIONAL_GETRI_HPP

#include <boost/assert.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/data_order.hpp>
#include <boost/numeric/bindings/detail/array.hpp>
#include <boost/numeric/bindings/is_complex.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/is_real.hpp>
#include <boost/numeric/bindings/lapack/workspace.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/traits/detail/utils.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_const.hpp>
#include <boost/utility/enable_if.hpp>

//
// The LAPACK-backend for getri is selected by defining a pre-processor
// variable, which can be one of
// * for ATLAS's CLAPACK, define BOOST_NUMERIC_BINDINGS_LAPACK_CLAPACK
// * netlib-compatible LAPACK is the default
//
#if defined BOOST_NUMERIC_BINDINGS_LAPACK_CLAPACK
#include <boost/numeric/bindings/lapack/detail/clapack.h>
#include <boost/numeric/bindings/lapack/detail/clapack_option.hpp>
#else
#include <boost/numeric/bindings/lapack/detail/lapack.h>
#include <boost/numeric/bindings/lapack/detail/lapack_option.hpp>
#endif

namespace boost {
namespace numeric {
namespace bindings {
namespace lapack {

//
// The detail namespace contains value-type-overloaded functions that
// dispatch to the appropriate back-end LAPACK-routine.
//
namespace detail {

#if defined BOOST_NUMERIC_BINDINGS_LAPACK_CLAPACK
//
// Overloaded function for dispatching to
// * ATLAS's CLAPACK backend, and
// * float value-type.
//
template< typename Order >
inline std::ptrdiff_t getri( Order, const int n, float* a, const int lda,
        const int* ipiv, float* work, const int lwork ) {
    return clapack_sgetri( clapack_option< Order >::value, n, a, lda, ipiv );
}

//
// Overloaded function for dispatching to
// * ATLAS's CLAPACK backend, and
// * double value-type.
//
template< typename Order >
inline std::ptrdiff_t getri( Order, const int n, double* a, const int lda,
        const int* ipiv, double* work, const int lwork ) {
    return clapack_dgetri( clapack_option< Order >::value, n, a, lda, ipiv );
}

//
// Overloaded function for dispatching to
// * ATLAS's CLAPACK backend, and
// * complex<float> value-type.
//
template< typename Order >
inline std::ptrdiff_t getri( Order, const int n, std::complex<float>* a,
        const int lda, const int* ipiv, std::complex<float>* work,
        const int lwork ) {
    return clapack_cgetri( clapack_option< Order >::value, n, a, lda, ipiv );
}

//
// Overloaded function for dispatching to
// * ATLAS's CLAPACK backend, and
// * complex<double> value-type.
//
template< typename Order >
inline std::ptrdiff_t getri( Order, const int n, std::complex<double>* a,
        const int lda, const int* ipiv, std::complex<double>* work,
        const int lwork ) {
    return clapack_zgetri( clapack_option< Order >::value, n, a, lda, ipiv );
}

#else
//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * float value-type.
//
template< typename Order >
inline std::ptrdiff_t getri( Order, const fortran_int_t n, float* a,
        const fortran_int_t lda, const fortran_int_t* ipiv, float* work,
        const fortran_int_t lwork ) {
    BOOST_STATIC_ASSERT( (is_same<Order, tag::column_major>::value) );
    fortran_int_t info(0);
    LAPACK_SGETRI( &n, a, &lda, ipiv, work, &lwork, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * double value-type.
//
template< typename Order >
inline std::ptrdiff_t getri( Order, const fortran_int_t n, double* a,
        const fortran_int_t lda, const fortran_int_t* ipiv, double* work,
        const fortran_int_t lwork ) {
    BOOST_STATIC_ASSERT( (is_same<Order, tag::column_major>::value) );
    fortran_int_t info(0);
    LAPACK_DGETRI( &n, a, &lda, ipiv, work, &lwork, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<float> value-type.
//
template< typename Order >
inline std::ptrdiff_t getri( Order, const fortran_int_t n,
        std::complex<float>* a, const fortran_int_t lda,
        const fortran_int_t* ipiv, std::complex<float>* work,
        const fortran_int_t lwork ) {
    BOOST_STATIC_ASSERT( (is_same<Order, tag::column_major>::value) );
    fortran_int_t info(0);
    LAPACK_CGETRI( &n, a, &lda, ipiv, work, &lwork, &info );
    return info;
}

//
// Overloaded function for dispatching to
// * netlib-compatible LAPACK backend (the default), and
// * complex<double> value-type.
//
template< typename Order >
inline std::ptrdiff_t getri( Order, const fortran_int_t n,
        std::complex<double>* a, const fortran_int_t lda,
        const fortran_int_t* ipiv, std::complex<double>* work,
        const fortran_int_t lwork ) {
    BOOST_STATIC_ASSERT( (is_same<Order, tag::column_major>::value) );
    fortran_int_t info(0);
    LAPACK_ZGETRI( &n, a, &lda, ipiv, work, &lwork, &info );
    return info;
}

#endif
} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to getri.
//
template< typename Value, typename Enable = void >
struct getri_impl {};

//
// This implementation is enabled if Value is a real type.
//
template< typename Value >
struct getri_impl< Value, typename boost::enable_if< is_real< Value > >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename VectorIPIV, typename WORK >
    static std::ptrdiff_t invoke( MatrixA& a, const VectorIPIV& ipiv,
            detail::workspace1< WORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixA >::type order;
        BOOST_STATIC_ASSERT( (bindings::is_mutable< MatrixA >::value) );
        BOOST_ASSERT( bindings::size(ipiv) >= bindings::size_column(a) );
        BOOST_ASSERT( bindings::size(work.select(real_type())) >=
                min_size_work( bindings::size_column(a) ));
        BOOST_ASSERT( bindings::size_column(a) >= 0 );
        BOOST_ASSERT( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        BOOST_ASSERT( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(a)) );
        return detail::getri( order(), bindings::size_column(a),
                bindings::begin_value(a), bindings::stride_major(a),
                bindings::begin_value(ipiv),
                bindings::begin_value(work.select(real_type())),
                bindings::size(work.select(real_type())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename VectorIPIV >
    static std::ptrdiff_t invoke( MatrixA& a, const VectorIPIV& ipiv,
            minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixA >::type order;
        bindings::detail::array< real_type > tmp_work( min_size_work(
                bindings::size_column(a) ) );
        return invoke( a, ipiv, workspace( tmp_work ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename VectorIPIV >
    static std::ptrdiff_t invoke( MatrixA& a, const VectorIPIV& ipiv,
            optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixA >::type order;
#if defined BOOST_NUMERIC_BINDINGS_LAPACK_CLAPACK
        bindings::detail::array< real_type > tmp_work( 0 );
#else
        real_type opt_size_work;
        detail::getri( order(), bindings::size_column(a),
                bindings::begin_value(a), bindings::stride_major(a),
                bindings::begin_value(ipiv), &opt_size_work, -1 );
        bindings::detail::array< real_type > tmp_work(
                traits::detail::to_int( opt_size_work ) );
#endif
        return invoke( a, ipiv, workspace( tmp_work ) );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
#if defined BOOST_NUMERIC_BINDINGS_LAPACK_CLAPACK
        return 0;
#else
        return std::max< std::ptrdiff_t >(1,n);
#endif
    }
};

//
// This implementation is enabled if Value is a complex type.
//
template< typename Value >
struct getri_impl< Value, typename boost::enable_if< is_complex< Value > >::type > {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;

    //
    // Static member function for user-defined workspaces, that
    // * Deduces the required arguments for dispatching to LAPACK, and
    // * Asserts that most arguments make sense.
    //
    template< typename MatrixA, typename VectorIPIV, typename WORK >
    static std::ptrdiff_t invoke( MatrixA& a, const VectorIPIV& ipiv,
            detail::workspace1< WORK > work ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixA >::type order;
        BOOST_STATIC_ASSERT( (bindings::is_mutable< MatrixA >::value) );
        BOOST_ASSERT( bindings::size(ipiv) >= bindings::size_column(a) );
        BOOST_ASSERT( bindings::size(work.select(value_type())) >=
                min_size_work( bindings::size_column(a) ));
        BOOST_ASSERT( bindings::size_column(a) >= 0 );
        BOOST_ASSERT( bindings::size_minor(a) == 1 ||
                bindings::stride_minor(a) == 1 );
        BOOST_ASSERT( bindings::stride_major(a) >= std::max< std::ptrdiff_t >(1,
                bindings::size_column(a)) );
        return detail::getri( order(), bindings::size_column(a),
                bindings::begin_value(a), bindings::stride_major(a),
                bindings::begin_value(ipiv),
                bindings::begin_value(work.select(value_type())),
                bindings::size(work.select(value_type())) );
    }

    //
    // Static member function that
    // * Figures out the minimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member function
    // * Enables the unblocked algorithm (BLAS level 2)
    //
    template< typename MatrixA, typename VectorIPIV >
    static std::ptrdiff_t invoke( MatrixA& a, const VectorIPIV& ipiv,
            minimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixA >::type order;
        bindings::detail::array< value_type > tmp_work( min_size_work(
                bindings::size_column(a) ) );
        return invoke( a, ipiv, workspace( tmp_work ) );
    }

    //
    // Static member function that
    // * Figures out the optimal workspace requirements, and passes
    //   the results to the user-defined workspace overload of the 
    //   invoke static member
    // * Enables the blocked algorithm (BLAS level 3)
    //
    template< typename MatrixA, typename VectorIPIV >
    static std::ptrdiff_t invoke( MatrixA& a, const VectorIPIV& ipiv,
            optimal_workspace ) {
        namespace bindings = ::boost::numeric::bindings;
        typedef typename result_of::data_order< MatrixA >::type order;
#if defined BOOST_NUMERIC_BINDINGS_LAPACK_CLAPACK
        bindings::detail::array< value_type > tmp_work( 0 );
#else
        value_type opt_size_work;
        detail::getri( order(), bindings::size_column(a),
                bindings::begin_value(a), bindings::stride_major(a),
                bindings::begin_value(ipiv), &opt_size_work, -1 );
        bindings::detail::array< value_type > tmp_work(
                traits::detail::to_int( opt_size_work ) );
#endif
        return invoke( a, ipiv, workspace( tmp_work ) );
    }

    //
    // Static member function that returns the minimum size of
    // workspace-array work.
    //
    static std::ptrdiff_t min_size_work( const std::ptrdiff_t n ) {
#if defined BOOST_NUMERIC_BINDINGS_LAPACK_CLAPACK
        return 0;
#else
        return std::max< std::ptrdiff_t >(1,n);
#endif
    }
};


//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. In
// addition, if applicable, they are overloaded for user-defined workspaces.
// Calls to these functions are passed to the getri_impl classes. In the 
// documentation, most overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for getri. Its overload differs for
// * User-defined workspace
//
template< typename MatrixA, typename VectorIPIV, typename Workspace >
inline typename boost::enable_if< detail::is_workspace< Workspace >,
        std::ptrdiff_t >::type
getri( MatrixA& a, const VectorIPIV& ipiv, Workspace work ) {
    return getri_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, ipiv, work );
}

//
// Overloaded function for getri. Its overload differs for
// * Default workspace-type (optimal)
//
template< typename MatrixA, typename VectorIPIV >
inline typename boost::disable_if< detail::is_workspace< VectorIPIV >,
        std::ptrdiff_t >::type
getri( MatrixA& a, const VectorIPIV& ipiv ) {
    return getri_impl< typename bindings::value_type<
            MatrixA >::type >::invoke( a, ipiv, optimal_workspace() );
}

} // namespace lapack
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif
