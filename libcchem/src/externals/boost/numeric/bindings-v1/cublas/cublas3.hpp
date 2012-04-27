//
//  Copyright Toon Knapen and Kresimir Fresl
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_BINDINGS_CUBLAS_CUBLAS3_HPP
#define BOOST_BINDINGS_CUBLAS_CUBLAS3_HPP

#include <boost/numeric/bindings/cublas/cublas3_overloads.hpp>
#include <boost/numeric/bindings/cublas/traits.hpp>
#include <boost/numeric/bindings/cublas/exception.hpp>

namespace boost { namespace numeric { namespace bindings { namespace cublas {

  // C <- alpha * op (A) * op (B) + beta * C 
  // op (X) == X || X^T || X^H
  template < typename value_type, typename matrix_type_a, typename matrix_type_b, typename matrix_type_c >
  // ! CAUTION this function assumes that all matrices involved are column-major matrices
  void gemm(const char TRANSA, const char TRANSB, 
	    const value_type& alpha,
	    const matrix_type_a &a,
	    const matrix_type_b &b,
	    const value_type &beta,
	    matrix_type_c &c
	    )
  {


    // std::cout << traits::matrix_size1( a ) << std::endl;
    // std::cout << "a " << traits::matrix_storage(a) << std::endl;
    // std::cout << "b " << traits::matrix_storage(b) << std::endl;
    // std::cout << "c " << traits::matrix_storage(c) << std::endl;
      

    const int m = TRANSA == traits::NO_TRANSPOSE ? traits::matrix_size1( a ) : traits::matrix_size2( a ) ;
    const int n = TRANSB == traits::NO_TRANSPOSE ? traits::matrix_size2( b ) : traits::matrix_size1( b );
    const int k = TRANSA == traits::NO_TRANSPOSE ? traits::matrix_size2( a ) : traits::matrix_size1( a ) ;

    assert( k ==  ( TRANSB == traits::NO_TRANSPOSE ? traits::matrix_size1( b ) : traits::matrix_size2( b ) ) ) ;

    assert( m == traits::matrix_size1( c ) ); 
    assert( n == traits::matrix_size2( c ) ); 

    if (!(m && n && k)) return;

    const int lda = traits::leading_dimension( a );
    const int ldb = traits::leading_dimension( b );
    const int ldc = traits::leading_dimension( c );

    detail::gemm(TRANSA, TRANSB, m, n, k,
    		 alpha,
    		 traits::matrix_storage(a), lda,
    		 traits::matrix_storage(b), ldb,
    		 beta,
    		 traits::matrix_storage(c), ldc);
    check_status();
  }


  // C <- alpha * A * B + beta * C 
  template<typename value_type, class A, class B, class C>
  void gemm(const value_type& alpha,
	    const cublas::matrix_expression<A> &a,
	    const cublas::matrix_expression<B> &b,
	    const value_type &beta,
	    cublas::matrix_expression<C> &c) {
      namespace transpose = traits::transpose;
      // std::cout << a().size1() << " " << a().size2() << std::endl;
      // std::cout << b().size1() << " " << b().size2() << std::endl;
      // // BOOST_AUTO(const &a_, transpose::no_transpose(a));
      // std::cout << char(transpose::option(a)) << " "
      // 		<< char(transpose::option(b)) << std::endl;
      gemm(transpose::option(a), transpose::option(b),
      	   alpha, transpose::no_transpose(a), transpose::no_transpose(b),
      	   beta, c);
  }


  // C <- A * B 
  // ! CAUTION this function assumes that all matrices involved are column-major matrices
  template<class A, class B, class C> 
  void gemm(const A &a, const B &b, C &c) {
      // std::cout << a.size1() << std::endl;
      gemm(1, a, b, 0, c);
  }


  // C <- alpha * A * A^T + beta * C
  // C <- alpha * A^T * A + beta * C
  template < typename value_type, typename matrix_type_a, typename matrix_type_c >
  void syrk( char uplo, char trans, const value_type& alpha, const matrix_type_a& a,
             const value_type& beta, matrix_type_c& c) {
     const int n = traits::matrix_size1( c );
     assert( n == traits::matrix_size2( c ) );
     const int k = trans == traits::NO_TRANSPOSE ? traits::matrix_size2( a ) : traits::matrix_size1( a ) ;
     assert( n == traits::NO_TRANSPOSE ? traits::matrix_size1( a ) : traits::matrix_size2( a ) );
     const int lda = traits::leading_dimension( a );
     const int ldc = traits::leading_dimension( c );

     const value_type *a_ptr = traits::matrix_storage( a ) ;
     value_type *c_ptr = traits::matrix_storage( c ) ;

     detail::syrk( uplo, trans, n, k, alpha, a_ptr, lda, beta, c_ptr, ldc );
  } // syrk()


  // C <- alpha * A * A^H + beta * C
  // C <- alpha * A^H * A + beta * C
  template < typename real_type, typename matrix_type_a, typename matrix_type_c >
  void herk( char uplo, char trans, const real_type& alpha, const matrix_type_a& a,
             const real_type& beta, matrix_type_c& c) {
     typedef typename matrix_type_c::value_type value_type ;

     const int n = traits::matrix_size1( c );
     assert( n == traits::matrix_size2( c ) );
     const int k = trans == traits::NO_TRANSPOSE ? traits::matrix_size2( a ) : traits::matrix_size1( a ) ;
     assert( n == traits::NO_TRANSPOSE ? traits::matrix_size1( a ) : traits::matrix_size2( a ) );
     const int lda = traits::leading_dimension( a );
     const int ldc = traits::leading_dimension( c );

     const value_type *a_ptr = traits::matrix_storage( a ) ;
     value_type *c_ptr = traits::matrix_storage( c ) ;

     detail::herk( uplo, trans, n, k, alpha, a_ptr, lda, beta, c_ptr, ldc );
  } // herk()

  // B <- alpha * op( A^-1 )
  // B <- alpha * B op( A^-1 )
  // op( A ) = A, A^T, A^H
  template < class T, class A, class B >
  void trsm( char side, char uplo, char transa, char diag, T const& alpha, A const& a, B& b ) {
     const int m = traits::matrix_size1( b ) ;
     const int n = traits::matrix_size2( b ) ;
     assert( ( side=='L' && m==traits::matrix_size2( a ) && m==traits::matrix_size1( a ) ) ||
             ( side=='R' && n==traits::matrix_size2( a ) && n==traits::matrix_size1( a ) ) ) ;
     assert( side=='R' || side=='L' ) ;
     assert( uplo=='U' || uplo=='L' ) ;
     assert( ( side=='L' && m==traits::matrix_size1( a ) ) || ( side=='R' && n==traits::matrix_size1( a ) ) ) ;
     detail::trsm( side, uplo, transa, diag, m, n, alpha,
                   traits::matrix_storage( a ), traits::leading_dimension( a ),
                   traits::matrix_storage( b ), traits::leading_dimension( b )
                 ) ;
  }

}}}}

#endif // BOOST_BINDINGS_CUBLAS_CUBLAS3_HPP
