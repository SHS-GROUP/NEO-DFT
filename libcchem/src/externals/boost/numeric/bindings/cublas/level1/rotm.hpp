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

#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_LEVEL1_ROTM_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_LEVEL1_ROTM_HPP

#include <boost/assert.hpp>
#include <boost/numeric/bindings/begin.hpp>
#include <boost/numeric/bindings/is_mutable.hpp>
#include <boost/numeric/bindings/remove_imaginary.hpp>
#include <boost/numeric/bindings/size.hpp>
#include <boost/numeric/bindings/stride.hpp>
#include <boost/numeric/bindings/value_type.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/type_traits/remove_const.hpp>

//
// The BLAS-backend is selected by defining a pre-processor variable,
//  which can be one of
// * for CBLAS, define BOOST_NUMERIC_BINDINGS_BLAS_CBLAS
// * for CUBLAS, define BOOST_NUMERIC_BINDINGS_BLAS_CUBLAS
// * netlib-compatible BLAS is the default
//
//// #if defined BOOST_NUMERIC_BINDINGS_BLAS_CBLAS
//// #include <boost/numeric/bindings/blas/detail/cblas.h>
//// #include <boost/numeric/bindings/blas/detail/cblas_option.hpp>
//// #elif defined BOOST_NUMERIC_BINDINGS_BLAS_CUBLAS
#include <boost/numeric/bindings/blas/detail/cublas.h>
#include <boost/numeric/bindings/blas/detail/default_order.hpp>
#include <boost/numeric/bindings/cublas/exception.hpp>
#include <boost/numeric/bindings/blas/detail/blas_option.hpp>
//// #else
//// #include <boost/numeric/bindings/blas/detail/blas.h>
//// #include <boost/numeric/bindings/blas/detail/blas_option.hpp>
//// #endif

namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

//
// The detail namespace contains value-type-overloaded functions that
// dispatch to the appropriate back-end BLAS-routine.
//
namespace detail {
using blas::detail::blas_option;
using blas::detail::default_order;

//// #if defined BOOST_NUMERIC_BINDINGS_BLAS_CBLAS
//// //
//// // Overloaded function for dispatching to
//// // * CBLAS backend, and
//// // * float value-type.
//// //
//// inline void rotm( const int n, float* x, const int incx, float* y,
////         const int incy, float* param ) {
////     cblas_srotm( n, x, incx, y, incy, param );
//// }
//// 
//// //
//// // Overloaded function for dispatching to
//// // * CBLAS backend, and
//// // * double value-type.
//// //
//// inline void rotm( const int n, double* x, const int incx, double* y,
////         const int incy, double* param ) {
////     cblas_drotm( n, x, incx, y, incy, param );
//// }
//// 
//// #elif defined BOOST_NUMERIC_BINDINGS_BLAS_CUBLAS
//
// Overloaded function for dispatching to
// * CUBLAS backend, and
// * float value-type.
//
inline void rotm( const int n, float* x, const int incx, float* y,
        const int incy, float* param ) {
    cublasSrotm( n, x, incx, y, incy, param );
}

//
// Overloaded function for dispatching to
// * CUBLAS backend, and
// * double value-type.
//
inline void rotm( const int n, double* x, const int incx, double* y,
        const int incy, double* param ) {
    cublasDrotm( n, x, incx, y, incy, param );
}

//// #else
//// //
//// // Overloaded function for dispatching to
//// // * netlib-compatible BLAS backend (the default), and
//// // * float value-type.
//// //
//// inline void rotm( const fortran_int_t n, float* x, const fortran_int_t incx,
////         float* y, const fortran_int_t incy, float* param ) {
////     BLAS_SROTM( &n, x, &incx, y, &incy, param );
//// }
//// 
//// //
//// // Overloaded function for dispatching to
//// // * netlib-compatible BLAS backend (the default), and
//// // * double value-type.
//// //
//// inline void rotm( const fortran_int_t n, double* x, const fortran_int_t incx,
////         double* y, const fortran_int_t incy, double* param ) {
////     BLAS_DROTM( &n, x, &incx, y, &incy, param );
//// }
//// 
//// #endif

} // namespace detail

//
// Value-type based template class. Use this class if you need a type
// for dispatching to rotm.
//
template< typename Value >
struct rotm_impl {

    typedef Value value_type;
    typedef typename remove_imaginary< Value >::type real_type;
    typedef void result_type;

    //
    // Static member function that
    // * Deduces the required arguments for dispatching to BLAS, and
    // * Asserts that most arguments make sense.
    //
    template< typename VectorX, typename VectorY, typename VectorPARAM >
    static result_type invoke( VectorX& x, VectorY& y, VectorPARAM& param ) {
        namespace bindings = ::boost::numeric::bindings;
        BOOST_STATIC_ASSERT( (is_same< typename remove_const<
                typename bindings::value_type< VectorX >::type >::type,
                typename remove_const< typename bindings::value_type<
                VectorY >::type >::type >::value) );
        BOOST_STATIC_ASSERT( (is_same< typename remove_const<
                typename bindings::value_type< VectorX >::type >::type,
                typename remove_const< typename bindings::value_type<
                VectorPARAM >::type >::type >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_mutable< VectorX >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_mutable< VectorY >::value) );
        BOOST_STATIC_ASSERT( (bindings::is_mutable< VectorPARAM >::value) );
        detail::rotm( bindings::size(x), bindings::begin_value(x),
                bindings::stride(x), bindings::begin_value(y),
                bindings::stride(y), bindings::begin_value(param) );
    }
};

//
// Functions for direct use. These functions are overloaded for temporaries,
// so that wrapped types can still be passed and used for write-access. Calls
// to these functions are passed to the rotm_impl classes. In the 
// documentation, the const-overloads are collapsed to avoid a large number of
// prototypes which are very similar.
//

//
// Overloaded function for rotm. Its overload differs for
//
template< typename VectorX, typename VectorY, typename VectorPARAM >
inline typename rotm_impl< typename bindings::value_type<
        VectorX >::type >::result_type
rotm( VectorX& x, VectorY& y, VectorPARAM& param ) {
    rotm_impl< typename bindings::value_type<
            VectorX >::type >::invoke( x, y, param );
    BOOST_NUMERIC_BINDINGS_CUBLAS_CHECK_STATUS();
}

} // namespace cublas
} // namespace bindings
} // namespace numeric
} // namespace boost

#endif

