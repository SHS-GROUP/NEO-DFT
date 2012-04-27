//
//  Copyright (C) Toon Knapen 2003
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS3_OVERLOADS_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS3_OVERLOADS_HPP

#include <boost/numeric/bindings/cublas/cublas.h>
#include <boost/numeric/bindings/traits/type_traits.hpp>

namespace boost { namespace numeric { namespace bindings { namespace cublas { namespace detail {

  using namespace boost::numeric::bindings::traits ;

  inline
  void gemm( char TRANSA, char TRANSB, const int& m, const int& n, const int& k,
	     const float & alpha,
	     const float * a_ptr, const int& lda,
	     const float * b_ptr, const int& ldb,
	     const float & beta,
	     float * c_ptr, const int& ldc )
  {
      cublasSgemm( TRANSA, TRANSB, m, n, k, alpha, a_ptr, lda, b_ptr, ldb, beta, c_ptr, ldc );
  }

  inline
  void gemm( char TRANSA, char TRANSB, const int& m, const int& n, const int& k,
	     const double & alpha,
	     const double * a_ptr, const int& lda,
	     const double * b_ptr, const int& ldb,
	     const double & beta,
	     double * c_ptr, const int& ldc )
  {
      cublasDgemm( TRANSA, TRANSB, m, n, k, alpha, a_ptr, lda, b_ptr, ldb, beta, c_ptr, ldc );
  }

  inline
  void gemm( char TRANSA, char TRANSB, const int& m, const int& n, const int& k,
	     const complex & alpha,
	     const complex * a_ptr, const int& lda,
	     const complex * b_ptr, const int& ldb,
	     const complex & beta,
	     complex * c_ptr, const int& ldc )
  {
      cublasCgemm( TRANSA, TRANSB, m, n, k,
		   alpha, a_ptr, lda, b_ptr, ldb,
		   beta, c_ptr, ldc );
  }

  //
  // syrk
  //
  inline
  void syrk( char uplo, char trans, const int& n, const int& k, const float& alpha,
             const float* a_ptr, const int lda, const float& beta, float* c_ptr,
             const int& ldc)
  {
     cublasSsyrk( uplo, trans, n, k, alpha, a_ptr, lda, beta, c_ptr, ldc);
  }

  inline
  void syrk( char uplo, char trans, const int& n, const int& k, const double& alpha,
             const double* a_ptr, const int lda, const double& beta, double* c_ptr,
             const int& ldc)
  {
     cublasDsyrk( uplo, trans, n, k, alpha, a_ptr, lda, beta, c_ptr, ldc);
  }

  inline
  void syrk( char uplo, char trans, const int& n, const int& k, const complex& alpha,
             const complex* a_ptr, const int lda, const complex& beta, complex* c_ptr,
             const int& ldc)
  {
     cublasCsyrk( uplo, trans, n, k, alpha, a_ptr, lda, beta, c_ptr, ldc);
  }

  //
  // HERK
  //
  inline
  void herk( char uplo, char trans, const int& n, const int& k, const float& alpha,
             const float* a_ptr, const int lda, const float& beta, float* c_ptr,
             const int& ldc)
  {
     cublasSsyrk( uplo, trans, n, k, alpha, a_ptr, lda, beta, c_ptr, ldc);
  }

  inline
  void herk( char uplo, char trans, const int& n, const int& k, const double& alpha,
             const double* a_ptr, const int lda, const double& beta, double* c_ptr,
             const int& ldc)
  {
     cublasDsyrk( uplo, trans, n, k, alpha, a_ptr, lda, beta, c_ptr, ldc);
  }


  inline
  void herk( char uplo, char trans, const int& n, const int& k,
	     const float& alpha, const complex* a_ptr, const int lda,
	     const float& beta, complex* c_ptr, const int& ldc)
  {
      struct value {
	  value(const float& data) : data_(data) {}
	  operator float() const { return data_; }
	  operator complex() const { return make_cuFloatComplex(data_, 0); }
      private: float data_;
      };
      cublasCherk( uplo, trans, n, k,
		   value(alpha), a_ptr, lda,
		   value(beta), c_ptr, ldc);
  }

  //
  // trsm
  //
  inline
  void trsm( char side, char uplo, char transa, char diag, int m, int n,
             float const& alpha, float const* a_ptr, int lda,
             float* b_ptr, int ldb )
  {
     cublasStrsm( side, uplo, transa, diag, m, n, alpha, a_ptr, lda, b_ptr, ldb ) ;
  }

  inline
  void trsm( char side, char uplo, char transa, char diag, int m, int n,
             double const& alpha, double const* a_ptr, int lda,
             double* b_ptr, int ldb )
  {
     cublasDtrsm( side, uplo, transa, diag, m, n, alpha, a_ptr, lda, b_ptr, ldb ) ;
  }

  inline
  void trsm( char side, char uplo, char transa, char diag, int m, int n,
             complex const& alpha, complex const* a_ptr, int lda,
             complex* b_ptr, int ldb )
  {
     cublasCtrsm( side, uplo, transa, diag, m, n, alpha, a_ptr, lda, b_ptr, ldb ) ;
  }


#ifdef CUBLAS_HAVE_DOUBLE_COMPLEX

  inline
  void gemm( char TRANSA, char TRANSB, const int& m, const int& n, const int& k,
	     const double_complex & alpha,
	     const double_complex * a_ptr, const int& lda,
	     const double_complex * b_ptr, const int& ldb,
	     const double_complex & beta,
	     double_complex * c_ptr, const int& ldc )
  {
      cublasZgemm( TRANSA, TRANSB, m, n, k,
		   alpha, a_ptr, lda, b_ptr, ldb,
		   beta, c_ptr, ldc );
  }

  inline
  void syrk( char uplo, char trans, const int& n, const int& k, const double_complex& alpha,
             const double_complex* a_ptr, const int lda, const double_complex& beta, double_complex* c_ptr,
             const int& ldc)
  {
     cublasZsyrk( uplo, trans, n, k, alpha, a_ptr, lda, beta, c_ptr, ldc);
  }

  inline
  void herk( char uplo, char trans, const int& n, const int& k, const double& alpha,
             const double_complex* a_ptr, const int lda, const double& beta, double_complex* c_ptr,
             const int& ldc)
  {
     cublasZherk( uplo, trans, n, k, alpha, a_ptr, lda, beta, c_ptr, ldc);
  }

  inline
  void trsm( char side, char uplo, char transa, char diag, int m, int n,
             double_complex const& alpha, double_complex const* a_ptr, int lda,
             double_complex* b_ptr, int ldb )
  {
     cublasZtrsm( side, uplo, transa, diag, m, n, alpha, a_ptr, lda, b_ptr, ldb ) ;
  }

#endif // CUBLAS_HAVE_DOUBLE_COMPLEX

}}}}}

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS3_OVERLOADS_HPP

