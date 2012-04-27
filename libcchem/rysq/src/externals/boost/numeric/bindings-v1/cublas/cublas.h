//
//  Copyright (C) 2002, 2003 Si-Lab b.v.b.a and Toon Knapen 
//
// Distributed under the Boost Software License, Version 1.0.
// (See accompanying file LICENSE_1_0.txt or copy at
// http://www.boost.org/LICENSE_1_0.txt)
//

#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS_H
#define BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS_H

#include <cublas.h>
#include <boost/numeric/bindings/traits/type.h>
#include <complex>

namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

    typedef cuComplex complex;
    typedef cuDoubleComplex double_complex;

    // cuComplex complex(const std::complex<float> &f) {
    // 	return make_cuComplex(std::real(f), std::imag(f));
    // }

}
}
}
}

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_CUBLAS_H
