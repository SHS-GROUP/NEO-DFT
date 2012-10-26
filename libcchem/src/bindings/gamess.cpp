#include "parallel.hpp"
#include "runtime.hpp"
#include "exception.hpp"
#include "basis/basis.hpp"

#include "bindings/fortran.h"
#include "bindings/gamess.h"
#include "bindings/bindings.h"
#include "bindings/bindings.hpp"

#include <boost/static_assert.hpp>

namespace gamess {

    struct Task : parallel::Task {
	void reset() {
	    integer_t tag = 0;
	    FORTRAN_CALL(ddi_sync)(&tag);
	    FORTRAN_CALL(ddi_dlbreset)();
	    FORTRAN_CALL(ddi_sync)(&tag);
	}
	size_t next() {
	    integer_t next;
#pragma omp critical(comm)
	    FORTRAN_CALL(ddi_dlbnext)(&next);
	    // std::cout <<  "next task: " << next << std::endl;
	    return next;
	}
    };

    struct Node : parallel::Comm {
	size_t rank() const {
	    integer_t size, rank;
	    FORTRAN_CALL(ddi_smp_nproc)(&size, &rank);
	    // std::cout << rank << std::endl;
	    return rank;
	}
	size_t size() const {
	    integer_t size, rank;
	    FORTRAN_CALL(ddi_smp_nproc)(&size, &rank);
	    // std::cout << size << std::endl;
	    return size;
	}
	void barrier() const {
	    throw CCHEM_EXCEPTION("not implemented");
	}
	void reduce(std::string op, int *data, size_t size) const {
	    throw CCHEM_EXCEPTION("not implemented");
	}
	void reduce(std::string op, double *data, size_t size) const {
	    throw CCHEM_EXCEPTION("not implemented");
	}
	void broadcast(double *data, size_t size, int root) const {
	    throw CCHEM_EXCEPTION("not implemented");
	};
    };

    struct Cluster : parallel::Comm {
	size_t rank() const {
	    integer_t size, rank;
	    FORTRAN_CALL(ddi_nproc)(&size, &rank);
	    size_t cluster = rank;
	    FORTRAN_CALL(ddi_smp_nproc)(&size, &rank);
	    cluster /= size;
	    // std::cout << rank << std::endl;
	    return cluster;
	}
	size_t size() const {
	    integer_t size, rank;
	    FORTRAN_CALL(ddi_nproc)(&size, &rank);
	    size_t cluster = size;
	    FORTRAN_CALL(ddi_smp_nproc)(&size, &rank);
	    assert(cluster%size == 0);
	    cluster /= size;
	    // std::cout << size << std::endl;
	    return size;
	}
	void barrier() const {
	    double tag = rank();
	    reduce("+", &tag, 1);
	    //throw CCHEM_EXCEPTION("not implemented");
	}
	void reduce(std::string op, int *data, size_t size) const {
	    if (op != "+") throw CCHEM_EXCEPTION("not implemented: " + op);
	    integer_t tag = 0, n = size;
	    std::vector<integer_t> data_(data, data+size); 
	    FORTRAN_CALL(ddi_masters_gsumi)(&tag, &data_[0], &n);
	    std::copy(data_.begin(), data_.end(), data);
	}
	void reduce(std::string op, double *data, size_t size) const {
	    if (op != "+") throw CCHEM_EXCEPTION("not implemented: " + op);
	    integer_t tag = 0, n = size;
	    FORTRAN_CALL(ddi_masters_gsumf)(&tag, data, &n);
	}
	void broadcast(double *data, size_t size, int root) const {
	    integer_t tag = 0, n = size, r = root;
	    char type = 'F';
	    FORTRAN_CALL(ddi_masters_bcast)(&tag, &type, data, &n, &r);
	}
    };

}

Parallel::Parallel() 
    : node_(new gamess::Node()),
      cluster_(new gamess::Cluster()),
      task_(new gamess::Task()) {}

size_t Parallel::rank() const {
    integer_t size, rank;
    FORTRAN_CALL(ddi_nproc)(&size, &rank);
    return  rank;
}

size_t Parallel::size() const {
    integer_t size, rank;
    FORTRAN_CALL(ddi_nproc)(&size, &rank);
    return size;
}

void Parallel::barrier() const {
    integer_t tag = 0;
    FORTRAN_CALL(ddi_sync)(&tag);
}

void Parallel::reduce(std::string op, double *data, size_t size) const {
    if (op == "max" && size == 1) {
	std::vector<double> r(this->size(), 0);
	r[this->rank()] = *data;
	reduce("+", r);
	*data = *std::max_element(r.begin(), r.end());
	return;
    }
    if (op != "+") throw CCHEM_EXCEPTION("not implemented: " + op);
    integer_t tag = 0, n = size;
    FORTRAN_CALL(ddi_gsumf)(&tag, data, &n);
}

void Parallel::reduce(std::string op, int *data, size_t size) const {
    if (op != "+") throw CCHEM_EXCEPTION("not implemented: " + op);
    integer_t tag = 0, n = size;
    std::vector<integer_t> data_(data, data+size); 
    FORTRAN_CALL(ddi_gsumi)(&tag, &data_[0], &n);
    std::copy(data_.begin(), data_.end(), data);
}

void Parallel::broadcast(int *data, size_t size, int root) const {
    BOOST_STATIC_ASSERT(sizeof(integer_t)%sizeof(int) == 0);
    CCHEM_ASSERT(size == 1);
    integer_t tag = 0, n = 1, r = root;
    int data_[sizeof(integer_t)/sizeof(int)];
    data_[0] = data[0];
    char type = 'I';
    FORTRAN_CALL(ddi_bcast)(&tag, &type, data_, &n, &r);
    data[0] = data_[0];
}

void Parallel::broadcast(size_t *data, size_t size, int root) const {
    BOOST_STATIC_ASSERT((sizeof(size_t)%sizeof(integer_t) == 0));
    integer_t tag = 0, n = size, r = root;
    char type = 'I';
    FORTRAN_CALL(ddi_bcast)(&tag, &type, data, &n, &r);
}

void Parallel::broadcast(double *data, size_t size, int root) const {
    integer_t tag = 0, n = size, r = root;
    char type = 'F';
    FORTRAN_CALL(ddi_bcast)(&tag, &type, data, &n, &r);    
}

extern "C" {

    /** \brief creates a new molecule from games, common block
	\return new molecule handle
    */
    integer_t molecule_new_gamess() {

	// BOOST_AUTO(cout, Runtime::rt().cout());
	// cout << Runtime::rt() << std::endl;

	int ian[INFOA.nat];

	std::copy(INFOA.ian, INFOA.ian + INFOA.nat, ian);

	//Create a molecule

	int mol = New_molecule(INFOA.nat, ian, (double*) INFOA.c, INFOA.zan);

	//std::cout << Molecule_find_object(mol);

	return mol;
    }

    integer_t molecule_new_gamess_();
#pragma weak molecule_new_gamess_ = molecule_new_gamess

    /** \brief creates a new basis from games, common block end a given molecule
	\param molecule molecule handle
	\return new basis handle
    */
    integer_t basis_new_gamess(integer_t *molecule) {
	int basis = New_basis(*molecule);
	//std::cout << basis << std::endl;
	Basis* pb = bindings::any().find<Basis*>(basis);
	assert(pb && "Basis not found");

	double *C[7] = { NSHEL.cs, NSHEL.cp, NSHEL.cd, NSHEL.cf, NSHEL.cg, NSHEL.ch, NSHEL.ci };

	for (int i = 0; i < NSHEL.nshell; ++i) {
	    int kstart = NSHEL.kstart[i] - 1;
	    int K = NSHEL.kng[i];
	    int L = NSHEL.ktype[i] - 1;
	    int center = NSHEL.katom[i] - 1;

	    if (NSHEL.kmin[i] == 1 && NSHEL.kmax[i] == 4) {
		pb->add(K, &NSHEL.ex[kstart], &C[0][kstart], &C[1][kstart], center);
	    }
	    else {
		pb->add(K, &NSHEL.ex[kstart], &C[L][kstart], L, center);
	    }
	}

	// std::cout << *basisObj;
	return basis;
    }

    integer_t basis_new_gamess_(integer_t *molecule);
#pragma weak basis_new_gamess_ = basis_new_gamess

}
