#include "bindings/bindings.h"
#include "bindings/fortran.h"
#include "bindings/bindings.hpp"

#include "basis/basis.hpp"
#include "basis/molecule.hpp"
#include "core/wavefunction.hpp"
#include "core/integral.hpp"

#include "runtime.hpp"
#include "mp2/mp2.hpp"
#include "parallel/context.hpp"

#include "utility/any.hpp"
#include "bindings/storage.hpp"

#include <iostream>
#include <iterator>
#include <algorithm>

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
    if (t == 'i') CXX_cout_impl<Integer>(size, data);
    if (t == 'c') CXX_cout_impl<char>(size, data);
}

// void Runtime_set_int(const char *key, int value) {
//     //runtime::set_default(key, value);
// }

void Delete_molecule(int key) {
    delete_<Molecule*>(key);
}

void Delete_basis(int key) {
    delete_<Basis*>(key);
}

void Delete_wavefunction(int key) {
    delete_<Wavefunction*>(key);
}

void Delete_integral_screening(int key) {
    delete_<Integral::Screening*>(key);
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

int New_wavefunction(int basis,
		     size_t core, size_t active, size_t virtuals,
		     bool axm, const double* C, size_t ldC,
		     const double* E) {
    const Basis *basis_ = bindings::any().find<Basis*>(basis);
    assert(basis);

    namespace ublas = boost::numeric::ublas;
    typedef ublas::indirect_array<> index_type;
    index_type bf = basis_->index<index_type>();
    index_type mo(core + active + virtuals);
    for (size_t i = 0; i < mo.size(); ++i) mo(i) = i;
    
    bindings::fortran_matrix<const double> C_(ldC, (axm ? mo : bf).size(), C);
    bindings::vector<const double> E_(mo.size(), E);

    Wavefunction *wavefunction = axm ?
	new Wavefunction(*basis_, core, active, virtuals,
			 ublas::project(C_, bf, mo),
			 E_) :
	new Wavefunction(*basis_, core, active, virtuals,
			 ublas::trans(ublas::project(C_, mo, bf)),
			 E_);
    
    wavefunction->normalize();
    return bindings::any().add(wavefunction);
}

int New_integral_screening(int basis, size_t N,
			   const double* K, double cutoff) {
    namespace ublas = boost::numeric::ublas;
    bindings::symmetric_fortran_matrix<const double> K_(N, K);
    typedef ublas::indirect_array<> index;
    index shells = bindings::any().find<Basis*>(basis)->shell_index<index>();
    Integral::Screening *screening =
	new Integral::Screening(ublas::project(K_, shells, shells), cutoff);
    return bindings::any().add(screening);
}

double MP2_energy(int wf, int screening) {
    if (parallel::Context().rank() != 0) return 0;
    return mp2::energy(*bindings::any().find<Wavefunction*>(wf),
		       *bindings::any().find<Integral::Screening*>(screening));
}

