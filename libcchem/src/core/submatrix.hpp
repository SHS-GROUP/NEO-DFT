#ifndef SUBMATRIX_HPP
#define SUBMATRIX_HPP

#include "basis/basis.hpp"
#include "core/wavefunction.hpp"
#include <boost/numeric/ublas/matrix_proxy.hpp>

template<class M>
boost::numeric::ublas::matrix_range<M>
submatrix(M& A, const Basis::Shell& ao, const Wavefunction::Orbitals &mo) {
    namespace ublas = boost::numeric::ublas;
    ublas::range r1(ao.start(), ao.stop()), r2(mo.start(), mo.stop());
    return ublas::matrix_range<M>(A, r1, r2);
}

template<class M>
boost::numeric::ublas::matrix_range<M>
submatrix(M& A, const Basis &ao, const Wavefunction::Orbitals &mo) {
    namespace ublas = boost::numeric::ublas;
    ublas::range r1(0, ao.size()), r2(mo.start(), mo.stop());
    return ublas::matrix_range<M>(A, r1, r2);
}

template<class M>
boost::numeric::ublas::matrix_range<M>
submatrix(M& A, size_t m, const Basis::Shell& ao) {
    namespace ublas = boost::numeric::ublas;
    ublas::range r1(0, m), r2(ao.start(), ao.stop());
    return ublas::matrix_range<M>(A, r1, r2);
}


#endif // SUBMATRIX_HPP
