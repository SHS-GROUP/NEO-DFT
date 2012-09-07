#ifndef QC_UTIL_RECTANGULAR_H
#define QC_UTIL_RECTANGULAR_H


#define BOOST_UBLAS_SHALLOW_ARRAY_ADAPTOR
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/symmetric.hpp>
#include <boost/numeric/ublas/io.hpp>

typedef ublas::matrix<double, ublas::column_major> RectangularMatrix;

#endif
