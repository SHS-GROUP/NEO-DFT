#ifndef BLAS_HPP
#define BLAS_HPP

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/storage_adaptors.hpp>

#define BOOST_NUMERIC_BINDINGS_BLAS_CBLAS
#ifdef HAVE_MKL
#define BOOST_NUMERIC_BINDINGS_BLAS_MKL
#endif

#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/matrix_expression.hpp>

#ifdef HAVE_CUBLAS
#include <boost/numeric/bindings/cublas/cublas.hpp>
#include <boost/numeric/bindings/cublas/vector.hpp>
#include <boost/numeric/bindings/cublas/matrix.hpp>
#endif

#ifdef HAVE_MKL
#ifndef mkl_free_buffers
#define mkl_free_buffers MKL_FreeBuffers
#endif
#endif

#ifdef _OPENMP
#include <omp.h>
#endif

namespace blas {

    using namespace boost::numeric::bindings::blas;

    inline void set_num_threads(size_t size) {
// #ifdef HAVE_MKL
// 	mkl_set_num_threads(size);
#if defined(_OPENMP)
	omp_set_num_threads(size);
#endif
    }

}

#ifdef HAVE_CUBLAS
namespace cublas {
    using namespace boost::numeric::bindings::cublas;
}
#endif

#endif // BLAS_HPP
