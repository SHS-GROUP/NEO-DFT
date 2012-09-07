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

#ifndef BOOST_NUMERIC_BINDINGS_BLAS_DETAIL_MKL_CBLAS_H
#define BOOST_NUMERIC_BINDINGS_BLAS_DETAIL_MKL_CBLAS_H

extern "C" {

#include <mkl_cblas.h>
#include <mkl_service.h>
#undef P4 // mkl_types.h defines P4 macro which breaks MPL

}

#endif 
