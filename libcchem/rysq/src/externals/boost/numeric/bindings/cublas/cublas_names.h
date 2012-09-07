//
//  Copyright (C) 2002, 2003 Si-Lab b.v.b.a and Toon Knapen 
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS_NAMES_H
#define BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS_NAMES_H

#include <boost/numeric/bindings/traits/fortran.h>

//
// level 1
//
#define CUBLAS_SSCAL FORTRAN_ID( sscal )
#define CUBLAS_DSCAL FORTRAN_ID( dscal )
#define CUBLAS_CSCAL FORTRAN_ID( cscal )
#define CUBLAS_ZSCAL FORTRAN_ID( zscal )

#define CUBLAS_SAXPY FORTRAN_ID( saxpy )
#define CUBLAS_DAXPY FORTRAN_ID( daxpy )
#define CUBLAS_CAXPY FORTRAN_ID( caxpy )
#define CUBLAS_ZAXPY FORTRAN_ID( zaxpy )

#define CUBLAS_SDOT  FORTRAN_ID( sdot )
#define CUBLAS_DDOT  FORTRAN_ID( ddot )

#define CUBLAS_CDOTU FORTRAN_ID( cdotu )
#define CUBLAS_ZDOTU FORTRAN_ID( zdotu )

#define CUBLAS_CDOTC FORTRAN_ID( cdotc )
#define CUBLAS_ZDOTC FORTRAN_ID( zdotc )

#define CUBLAS_SNRM2 FORTRAN_ID( snrm2 )
#define CUBLAS_DNRM2 FORTRAN_ID( dnrm2 )
#define CUBLAS_SCNRM2 FORTRAN_ID( scnrm2 )
#define CUBLAS_DZNRM2 FORTRAN_ID( dznrm2 )

#define CUBLAS_SASUM FORTRAN_ID( sasum )
#define CUBLAS_DASUM FORTRAN_ID( dasum )
#define CUBLAS_SCASUM FORTRAN_ID( scasum )
#define CUBLAS_DZASUM FORTRAN_ID( dzasum )

#define CUBLAS_SCOPY FORTRAN_ID( scopy )
#define CUBLAS_DCOPY FORTRAN_ID( dcopy )
#define CUBLAS_CCOPY FORTRAN_ID( ccopy )
#define CUBLAS_ZCOPY FORTRAN_ID( zcopy )

//
// level 2
//
#define CUBLAS_SGEMV FORTRAN_ID( sgemv )
#define CUBLAS_DGEMV FORTRAN_ID( dgemv )
#define CUBLAS_CGEMV FORTRAN_ID( cgemv )
#define CUBLAS_ZGEMV FORTRAN_ID( zgemv )

#define CUBLAS_SGER  FORTRAN_ID( sger )
#define CUBLAS_DGER  FORTRAN_ID( dger )

#define CUBLAS_CGERU FORTRAN_ID( cgeru )
#define CUBLAS_ZGERU FORTRAN_ID( zgeru )

#define CUBLAS_CGERC FORTRAN_ID( cgerc )
#define CUBLAS_ZGERC FORTRAN_ID( zgerc )

//
// level 3
//
#define CUBLAS_SGEMM FORTRAN_ID( sgemm )
#define CUBLAS_DGEMM FORTRAN_ID( dgemm )
#define CUBLAS_CGEMM FORTRAN_ID( cgemm )
#define CUBLAS_ZGEMM FORTRAN_ID( zgemm )

#define CUBLAS_SSYRK FORTRAN_ID( ssyrk )
#define CUBLAS_DSYRK FORTRAN_ID( dsyrk )
#define CUBLAS_CSYRK FORTRAN_ID( csyrk )
#define CUBLAS_ZSYRK FORTRAN_ID( zsyrk )
#define CUBLAS_CHERK FORTRAN_ID( cherk )
#define CUBLAS_ZHERK FORTRAN_ID( zherk )

#define CUBLAS_STRSM FORTRAN_ID( strsm )
#define CUBLAS_DTRSM FORTRAN_ID( dtrsm )
#define CUBLAS_CTRSM FORTRAN_ID( ctrsm )
#define CUBLAS_ZTRSM FORTRAN_ID( ztrsm )

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS_NAMES_H
