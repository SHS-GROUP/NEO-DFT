#include "bindings/bindings.h"
#include "bindings/bindings.hpp"
#include "bindings/storage.hpp"

#include "dft/dft.hpp"
#include "basis/normalize.hpp"

#include <boost/numeric/ublas/symmetric.hpp>
#include "boost/numeric/ublas/storage_adaptors.hpp"
#include <boost/typeof/typeof.hpp>

void DFT_fock(int wavefunction, double rcut, double ccut, int N,
	      const double *dr, int ldr, const double *w, int ldw,
	      double *F,
	      double *Xa, double *Xg, double *Ec, double *Rho,
	      const char *functional) {

    namespace ublas = boost::numeric::ublas;
    typedef ublas::symmetric_adaptor<dft::matrix_type, ublas::upper> symmetric;

    BOOST_AUTO(const &W, *bindings::any().find<Wavefunction*>(wavefunction));
    BOOST_AUTO(const &basis, W.basis());
    size_t size = basis.size();

    dft::matrix_type F_(size,size,0);

    //std::cout << rcut << std::endl;
    dft::Parameters p;
    p.rcut = rcut;
    p.ccut = ccut;
    dft::XC xc;
    xc = dft::fock(W, std::string(functional),
		   const_array_ref<double,3>(N, dr, ldr),
		   const_array_ref<double,1>(N, w, ldw),
		   F_, p);

    basis::normalize_matrix(W.basis(), F_);

    ublas::upper upper;
    ublas::column_major major;
    BOOST_AUTO(S, bindings::make_matrix(size, F, upper, major));

    typedef BOOST_TYPEOF(basis.index()) index_type;
    ublas::indirect_array<index_type> index(size);
    for (size_t i = 0; i < index.size(); ++i) {
	index[basis.index()[i]] = i;
    }
    S.plus_assign(ublas::project(symmetric(F_), index, index));

    *Xa += xc.Xa;
    *Xg += xc.Xg;
    *Ec += xc.Ec;
    *Rho += xc.rho;
}
