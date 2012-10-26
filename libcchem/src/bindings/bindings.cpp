#include "cchem.hpp"
#include "runtime.hpp"
#include "exception.hpp"

#include "bindings/bindings.h"
#include "bindings/fortran.h"
#include "bindings/bindings.hpp"

#include "basis/basis.hpp"
#include "basis/molecule.hpp"
#include "core/wavefunction.hpp"
#include "integrals/screening.hpp"
#include "mp2/mp2.hpp"

#include "utility/any.hpp"
#include "bindings/storage.hpp"

#include <iostream>
#include <iterator>
#include <algorithm>
#include <boost/numeric/ublas/adaptors.hpp>
#include <boost/numeric/ublas/io.hpp>

utility::any<int>& bindings::any() {
    static utility::any<int> any;
    return any;
}

namespace {

    template<typename T>
    void delete_(int key) {
	delete bindings::any().find<T>(key);
	bindings::any().erase(key);
    }

    template<typename T>
    void CXX_cout_impl(size_t size, const void *data) {
	const T *data_ = static_cast<const T*>(data);
	std::ostream_iterator<T> it(std::cout, " ");
	std::copy(data_, data_ + size, it);
	std::cout << std::endl;
    }

}


void CXX_cout(char type, size_t size, const void *data) {
    char t = tolower(type);
    if (t == 'd') CXX_cout_impl<double>(size, data);
    if (t == 'f') CXX_cout_impl<float>(size, data);
    if (t == 'i') CXX_cout_impl<integer_t>(size, data);
    if (t == 'c') CXX_cout_impl<char>(size, data);
}

void CChem_initialize() {
    cchem::initialize();
}

void CChem_runtime_set_int(const char *key, int value) {
    Runtime::rt().set(key, value);
}

void CChem_runtime_set_double(const char *key, double value) {
    //std::cout <<  key << ":" << value << std::endl;
    Runtime::rt().set(key, value);
}

void Delete_molecule(int key) {
    delete_<Molecule*>(key);
}

void Delete_basis(int key) {
    delete_<Basis*>(key);
}

void Delete_wavefunction(int key) {
    delete_<Wavefunction*>(key);
}

void Delete_integrals_screening(int key) {
    delete_<integrals::Screening*>(key);
}


int New_molecule(int size, int *Z, double *r3, double *Q) {
    Molecule* molecule = new Molecule(size, Z, r3, Q);
    return bindings::any().add(molecule);
}

int New_basis(int molecule) {
    const Molecule& m = *bindings::any().find<Molecule*>(molecule);
    return bindings::any().add(new Basis(m));
}

void Basis_sort(int basis) {
    bindings::any().find<Basis*>(basis)->sort();
}

void Basis_recenter(int basis, 
		    const double* x, int incx,
		    const double* y, int incy,
		    const double* z, int incz) {
      const double* r[] = { x,y,z };
      int inc[] = { incx, incy, incz };
      basis::recenter(*bindings::any().find<Basis*>(basis), r, inc);
}


int New_wavefunction(int basis) {
    const Basis *b = bindings::any().find<Basis*>(basis);
    assert(b);
    return bindings::any().add(new Wavefunction(*b));
}

void Wavefunction_set_virtuals(int wf, size_t start, size_t stop) {
    Wavefunction *w = bindings::any().find<Wavefunction*>(wf);
    w->set_virtuals(start, stop);
}

void Wavefunction_set_occupied(int wf, size_t start, size_t stop, size_t core) {
    Wavefunction *w = bindings::any().find<Wavefunction*>(wf);
    w->set_occupied(start, stop, core);
}

void Wavefunction_set_C(int wf, const char trans, size_t m, size_t n,
			const double* C, size_t ld) {
    namespace ublas = boost::numeric::ublas;
    BOOST_AUTO(c_, ublas::make_matrix<ublas::column_major>(ld, n, C));
    Wavefunction *w = bindings::any().find<Wavefunction*>(wf);
    ublas::range ra = ublas::range(0, w->basis().size());
    if (trans == 't' || trans == 'T') {
	ublas::range rm(0,m);
	BOOST_AUTO(c, ublas::project(c_, rm, ra));
	w->set_C(ublas::trans(c));
    }
    else {
	ublas::range rn(0,n);
	BOOST_AUTO(c, ublas::project(c_, ra, rn));
	w->set_C(c);
    }
}

void Wavefunction_normalize(int wf) {
    bindings::any().find<Wavefunction*>(wf)->normalize();
}

void Wavefunction_set_F(int wf, const double* F, size_t ldF) {
    using boost::numeric::ublas::make_matrix;
    using boost::numeric::ublas::column_major;
    Wavefunction *w = bindings::any().find<Wavefunction*>(wf);
    BOOST_AUTO(f, make_matrix<column_major>(ldF, w->virtuals().stop(), F));
    w->set_F(project(f, w->orbitals(), w->orbitals()));
}

void Wavefunction_set_e(int wf, const double* e) {
    using boost::numeric::ublas::make_vector;
    Wavefunction *w = bindings::any().find<Wavefunction*>(wf);
    w->set_e(make_vector(w->virtuals().stop(), e));
}


int New_integrals_screening(int basis, size_t N,
			   const double* K, double cutoff) {
    namespace ublas = boost::numeric::ublas;
    bindings::symmetric_fortran_matrix<const double> K_(N, K);
    typedef ublas::indirect_array<> index;
    index shells = bindings::any().find<Basis*>(basis)->shell_index<index>();

    // std::cout << N << std::endl;
    // for (size_t i = 0; i < shells.size(); ++i) {
    // 	std::cout << i << "->" << shells[i] << std::endl;
    // }

    integrals::Screening *screening =
	new integrals::Screening(ublas::project(K_, shells, shells), cutoff);
    return bindings::any().add(screening);
}

int CChem_mp2_energy(int wf, double *E) {
    Runtime &rt = Runtime::rt();
    try {
	*E = cchem::mp2::energy(*bindings::any().find<Wavefunction*>(wf), rt);
    }
    catch (const cchem::exception &e) {
	std::string what = e.what();
	if (!what.empty()) Parallel().cout() << what << std::endl;
	return 0;
    }
    return 1;
}

