//
//  Copyright Toon Knapen and Kresimir Fresl
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_BINDINGS_BLAS_BLAS3_HPP
#define BOOST_BINDINGS_BLAS_BLAS3_HPP

#include <boost/numeric/bindings/blas/blas3_overloads.hpp>
#include <boost/numeric/bindings/traits/traits.hpp>
#include <boost/numeric/bindings/traits/transpose.hpp>

namespace boost {
namespace numeric {
namespace bindings {
namespace blas {

    // C <- alpha * op (A) * op (B) + beta * C 
    // op (X) == X || X^T || X^H

    template < typename value_type, typename A,
	       typename B, typename C >
    // ! CAUTION this function assumes that all matrices involved are column-major matrices
    void gemm(const char TRANSA, const char TRANSB, 
	      const value_type& alpha,
	      const A &a,
	      const B &b,
	      const value_type &beta,
	      C &c
	      )
    {

	const int m = ((TRANSA == traits::NO_TRANSPOSE) ?
		       traits::matrix_size1( a ) : traits::matrix_size2( a ));
	const int n = ((TRANSB == traits::NO_TRANSPOSE) ?
		       traits::matrix_size2( b ) : traits::matrix_size1( b ));
	const int k = ((TRANSA == traits::NO_TRANSPOSE) ? 
		       traits::matrix_size2( a ) : traits::matrix_size1( a ));

	assert( k ==  ( TRANSB == traits::NO_TRANSPOSE ?
			traits::matrix_size1( b ) : traits::matrix_size2( b ) ) ) ;
	assert( m == traits::matrix_size1( c ) ); 
	assert( n == traits::matrix_size2( c ) ); 

	const int lda = traits::leading_dimension( a );
	const int ldb = traits::leading_dimension( b );
	const int ldc = traits::leading_dimension( c );

	const value_type *a_ptr = traits::matrix_storage( a ) ;
	const value_type *b_ptr = traits::matrix_storage( b ) ;
	value_type *c_ptr = traits::matrix_storage( c ) ;

	detail::gemm( TRANSA, TRANSB, m, n, k, alpha, a_ptr, lda, b_ptr, ldb, beta, c_ptr, ldc ) ;
    }


    // C <- alpha * A * B + beta * C 
    template < typename value_type, typename A, typename B, typename C >
    void gemm(const value_type& alpha,
	      const A &a,
	      const B &b,
	      const value_type &beta,
	      C &c
	      )
    {
	gemm( traits::NO_TRANSPOSE, traits::NO_TRANSPOSE, alpha, a, b, beta, c ) ;
    }


    // C <- A * B 
    // ! CAUTION this function assumes that all matrices involved are column-major matrices
    template<class A, class B, class C>
    void gemm(const A &a, const B &b, C &c) {
	typedef typename traits::matrix_traits<C>::value_type val_t; 
	gemm( traits::NO_TRANSPOSE, traits::NO_TRANSPOSE, (val_t) 1, a, b, (val_t) 0, c ) ;
    }


    // C <- alpha * A * A^T + beta * C
    // C <- alpha * A^T * A + beta * C
    template < typename value_type, typename A, typename C >
    void syrk( char uplo, char trans, const value_type& alpha, const A& a,
	       const value_type& beta, C& c) {
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
    template < typename real_type, typename A, typename C >
    void herk( char uplo, char trans, const real_type& alpha, const A& a,
	       const real_type& beta, C& c) {
	typedef typename C::value_type value_type ;

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

		
    template<typename T, class A, class B, class C>
    // ! CAUTION this function assumes that all matrices involved are column-major matrices
    void syr2k(const char uplo, const char trans, 
	       const T &alpha, const A &a, const B &b,
	       const T &beta, C &c) {

	const int n = traits::matrix_size1(c);
	const int k = ((trans == traits::NO_TRANSPOSE) ? 
		       traits::matrix_size2(a) : traits::matrix_size1(a));

	assert(n == (trans == traits::NO_TRANSPOSE ?
		     traits::matrix_size1(a) : traits::matrix_size2(a)));
	assert(n == (trans == traits::NO_TRANSPOSE ?
		     traits::matrix_size1(b) : traits::matrix_size2(b)));
	assert(k == (trans == traits::NO_TRANSPOSE ?
		     traits::matrix_size2(b) : traits::matrix_size1(b)));

	const int lda = traits::leading_dimension(a);
	const int ldb = traits::leading_dimension(b);
	const int ldc = traits::leading_dimension(c);

	typedef typename traits::matrix_traits<A>::value_type value_type;

	detail::syr2k(uplo, trans, n, k,
		      value_type(alpha),
		      traits::matrix_storage(a), lda,
		      traits::matrix_storage(b), ldb,
		      value_type(beta),
		      traits::matrix_storage(c), ldc);
    }		

    template<typename T, class A, class B, class C>
    void syr2k(const A &a, const B &b, T beta, C &c) {
	syr2k('U', traits::NO_TRANSPOSE, T(1), a, b, beta, c);   
    }


}
}
}
}

#endif // BOOST_BINDINGS_BLAS_BLAS3_HPP
