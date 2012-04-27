//
//  Copyright (C) Toon Knapen 2003
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS1_OVERLOADS_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS1_OVERLOADS_HPP

#include <boost/numeric/bindings/cublas/cublas.h>
#include <boost/numeric/bindings/traits/type.hpp>
#include <boost/numeric/bindings/traits/type_traits.hpp>

namespace boost { namespace numeric { namespace bindings { namespace cublas { namespace detail {

  using namespace boost::numeric::bindings::traits ;

  // x *= alpha 
  inline void scal(const int& n, const float&     alpha, float*     x, const int& incx)
  {
      cublasSscal( n, alpha, x, incx );
  } 
  inline void scal(const int& n, const double&    alpha, double*    x, const int& incx)
  {
      cublasDscal( n, alpha, x, incx );
  }
  inline void scal(const int& n, const complex& alpha, complex* x, const int& incx)
  {
      cublasCscal( n, alpha, ( x ), incx );
  }


  // y += alpha * x 
  inline void axpy(const int& n, const float    & alpha, const float    * x, const int& incx, float    * y, const int& incy)
  {
      cublasSaxpy( n, alpha, x, incx, y, incy );
  }
  inline void axpy(const int& n, const double   & alpha, const double   * x, const int& incx, double   * y, const int& incy)
  {
      cublasDaxpy( n, alpha, x, incx, y, incy );
  }
  inline void axpy(const int& n, const complex& alpha, const complex* x, const int& incx, complex* y, const int& incy)
  {
      cublasCaxpy( n, alpha, ( x ), incx, ( y ), incy );
  }


  // x^T . y 
  inline float  dot(const int& n, const float * x, const int& incx, const float * y, const int& incy)
  {
      return cublasSdot( n, x, incx, y, incy );
  }
  inline double dot(const int& n, const double* x, const int& incx, const double* y, const int& incy)
  {
      return cublasDdot( n, x, incx, y, incy );
  }

  // x^T . y
  inline complex dotu(const int& n, const complex* x, const int& incx, const complex* y, const int& incy)
  {
      return cublasCdotu(n, ( x ), incx, ( y ), incy );
  }


  // x^H . y
  inline complex dotc(const int& n, const complex* x, const int& incx, const complex* y, const int& incy)
  {
      return cublasCdotc(n, ( x ), incx, ( y ), incy );
  }


  // euclidean norm
  inline float  nrm2(const int& n, const float*   x, const int& incx)
  {
      return cublasSnrm2( n, x, incx );
  }
  inline double nrm2(const int& n, const double*  x, const int& incx)
  {
      return cublasDnrm2( n, x, incx );
  }
  inline float  nrm2(const int& n, const complex*   x, const int& incx)
  {
      return cublasScnrm2( n, (x), incx );
  }

  
  // 1-norm
  inline float  asum(const int& n, const float*   x, const int& incx)
  {
      return cublasSasum( n, x, incx );
  }
  inline double asum(const int& n, const double*  x, const int& incx)
  {
      return cublasDasum( n, x, incx );
  }
  inline float  asum(const int& n, const complex*   x, const int& incx)
  {
      return cublasScasum( n, (x), incx );
  }

  
  // copy
  inline void copy(const int& n, const float*     x, const int& incx, float*     y, const int& incy)
  {
      cublasScopy( n, x, incx, y, incy );
  }
  inline void copy(const int& n, const double*    x, const int& incx, double*    y, const int& incy)
  {
      cublasDcopy( n, x, incx, y, incy );
  }
  inline void copy(const int& n, const complex* x, const int& incx, complex* y, const int& incy)
  {
      cublasCcopy( n, (x), incx, (y), incy );
  }

#ifdef CUBLAS_HAVE_DOUBLE_COMPLEX

  inline void scal(const int& n, const double_complex &alpha, double_complex* x, const int& incx)
  {
      cublasZscal( n, alpha, ( x ), incx );
  }

  inline void axpy(const int& n, const double_complex& alpha, const double_complex* x, const int& incx, double_complex* y, const int& incy)
  {
      cublasZaxpy( n, alpha, ( x ), incx, ( y ), incy );
  }

  inline double_complex dotc(const int& n, const double_complex* x, const int& incx, const double_complex* y, const int& incy)
  {
      return cublasZdotc(n, ( x ), incx, ( y ), incy );
  }

  inline double_complex dotu(const int& n, const double_complex* x, const int& incx, const double_complex* y, const int& incy)
  {
      return cublasZdotu(n, ( x ), incx, ( y ), incy );
  }

  inline double nrm2(const int& n, const double_complex*  x, const int& incx)
  {
      return cublasDznrm2( n, (x), incx );
  }

  inline double asum(const int& n, const double_complex*  x, const int& incx)
  {
      return cublasDzasum( n, (x), incx );
  }

  inline void copy(const int& n, const double_complex* x, const int& incx, double_complex* y, const int& incy)
  {
      cublasZcopy( n, (x), incx, (y), incy );
  }

#endif // CUBLAS_HAVE_DOUBLE_COMPLEX

}}}}}

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS1_OVERLOADS_HPP

