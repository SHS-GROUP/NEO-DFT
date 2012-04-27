/*
 * 
 * Copyright (c) Kresimir Fresl 2002 
 *
 * Distributed under the Boost Software License, Version 1.0.
 * (See accompanying file LICENSE_1_0.txt or copy at
 * http://www.boost.org/LICENSE_1_0.txt)
 *
 * Author acknowledges the support of the Faculty of Civil Engineering, 
 * University of Zagreb, Croatia.
 *
 */

#ifndef BOOST_NUMERIC_BINDINGS_BLAS_DETAIL_CBLAS_H
#define BOOST_NUMERIC_BINDINGS_BLAS_DETAIL_CBLAS_H


#if defined BOOST_NUMERIC_BINDINGS_BLAS_MKL
#include <boost/numeric/bindings/blas/detail/mkl_cblas.h> 

#else

extern "C" {
#include <cblas.h>
}

#endif

#endif 
