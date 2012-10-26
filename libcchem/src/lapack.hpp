#ifndef LAPACK_HPP
#define LAPACK_HPP

#ifdef HAVE_CONFIG_H
#include "config.h"
#endif

#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/ublas/matrix_expression.hpp>

#ifdef LIBCCHEM_WITH_INTEGER8
#define BIND_FORTRAN_INTEGER_8
#endif


#include <boost/numeric/bindings/lapack/driver.hpp>

namespace lapack = boost::numeric::bindings::lapack;

typedef fortran_int_t lapack_int;

#endif // LAPACK_HPP
