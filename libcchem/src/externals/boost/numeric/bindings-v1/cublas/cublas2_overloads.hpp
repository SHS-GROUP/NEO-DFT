//
//  Copyright (C) Toon Knapen 2003
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS2_OVERLOADS_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS2_OVERLOADS_HPP

#include <boost/numeric/bindings/cublas/cublas.h>
#include <boost/numeric/bindings/traits/type_traits.hpp>

namespace boost { namespace numeric { namespace bindings { namespace cublas { namespace detail {

  using namespace boost::numeric::bindings::traits ;

  inline
  void gemv( char TRANS, const int& m, const int& n,
	     const float & alpha,
	     const float * a_ptr, const int& lda,
	     const float * x_ptr, const int& incx,
	     const float & beta,
	     float * y_ptr, const int& incy )
  {
      cublasSgemv( TRANS, m, n, alpha, a_ptr, lda, x_ptr, incx, beta, y_ptr, incy );
  }

  inline
  void gemv( char TRANS, const int& m, const int& n,
	     const double & alpha,
	     const double * a_ptr, const int& lda,
	     const double * x_ptr, const int& incx,
	     const double & beta,
	     double * y_ptr, const int& incy )
  {
      cublasDgemv( TRANS, m, n, alpha, a_ptr, lda, x_ptr, incx, beta, y_ptr, incy );
  }

  inline
  void gemv( char TRANS, const int& m, const int& n,
	     const complex & alpha,
	     const complex * a_ptr, const int& lda,
	     const complex * x_ptr, const int& incx,
	     const complex & beta,
	     complex * y_ptr, const int& incy )
  {
      cublasCgemv( TRANS, m, n, alpha, a_ptr, lda, x_ptr, incx, beta, y_ptr, incy );
  }

  inline
  void gemv( char TRANS, const int& m, const int& n,
	     const double_complex & alpha,
	     const double_complex * a_ptr, const int& lda,
	     const double_complex * x_ptr, const int& incx,
	     const double_complex & beta,
	     double_complex * y_ptr, const int& incy )
  {
      cublasZgemv( TRANS, m, n, alpha, a_ptr, lda, x_ptr, incx, beta, y_ptr, incy );
  }

  inline
  void ger( const int& m, const int& n,
	    const float & alpha,
	    const float * x_ptr, const int& incx,
	    const float * y_ptr, const int& incy,
	    float * a_ptr, const int& lda )
  {
      cublasSger( m, n, alpha, x_ptr, incx, y_ptr, incy, a_ptr, lda );
  }

  inline
  void ger( const int& m, const int& n,
	    const double & alpha,
	    const double * x_ptr, const int& incx,
	    const double * y_ptr, const int& incy,
	    double * a_ptr, const int& lda )
  {
      cublasDger( m, n, alpha, x_ptr, incx, y_ptr, incy, a_ptr, lda );
  }

/*  
  inline
  void geru( const int& m, const int& n, const complex & alpha, const complex * x_ptr, const int& incx, complex * y_ptr, const int& incy, complex * a_ptr, const int& lda ) { cublascgeru( m, &n, alpha, x_ptr, &incx, y_ptr, &incy, a_ptr, &lda ) ; }
  inline
  void geru( const int& m, const int& n, const double_complex & alpha, const double_complex * x_ptr, const int& incx, double_complex * y_ptr, const int& incy, double_complex * a_ptr, const int& lda ) { cublasZgeru( m, &n, alpha, x_ptr, &incx, y_ptr, &incy, a_ptr, &lda ) ; }
*/
}}}}}

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS2_OVERLOADS_HPP

