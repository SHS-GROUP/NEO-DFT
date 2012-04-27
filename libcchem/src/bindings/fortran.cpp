/**
   @file
   @brief libcchem  Fortran bindings
   @details The fortran binding provides lower case subroutines with and without
   trailing underscore for compatibility with different fortran compilers.
*/

#include "bindings/fortran.h"
#include "bindings/bindings.h"
#include "bindings/bindings.hpp"

#include <vector>
#include <string>
#include <algorithm>
#include <boost/algorithm/string/trim.hpp>

namespace {

    template<typename T>
    struct vector {
	template<typename U>
	vector(const U *data, size_t size) : data_(data, data + size) {}
	operator const T*() const { return &data_[0]; }
    private:
	std::vector<T> data_;
    };

    struct string {
	string(const char *data) : data_(data) {}
	string(const char *data, size_t size) :
	    data_(boost::trim_right_copy(std::string(data, size))) {}
	operator const char*() const { return data_.c_str(); }
    private:
	std::string data_;
    };

}

extern "C" {


#define FORTRAN_FUNCTION(F, A, CF, CA) F ## _ A { return CF CA; }
#define FORTRAN_SUBROUTINE(F, A, CF, CA) void F ## _ A { CF CA; }

    FORTRAN_SUBROUTINE(cxx_cout, (const char *type, const Integer *n,
				  const void *data),
		       CXX_cout, (*type, *n, data))

    // FORTRAN_SUBROUTINE(runtime_set_int,
    // 		       (const char *key, const Integer *value, Integer L),
    // 		       Runtime_set_int,
    // 		       (string(key, key+L), *value))
					 

    FORTRAN_SUBROUTINE(delete_molecule, (const Integer *key),
		       Delete_molecule, (*key))

    FORTRAN_SUBROUTINE(delete_wavefunction, (const Integer *key),
		       Delete_wavefunction, (*key))

    FORTRAN_SUBROUTINE(delete_basis, (const Integer *key),
		       Delete_basis, (*key))

    FORTRAN_SUBROUTINE(basis_sort, (const Integer *key),
		       Basis_sort, (*key))

    FORTRAN_SUBROUTINE(basis_recenter, (const Integer *basis,
					const double *x, const Integer *incx,
					const double *y, const Integer *incy,
					const double *z, const Integer *incz),
		       Basis_recenter, (*basis, x, *incx, y, *incy, z, *incz))


    FORTRAN_FUNCTION(Integer new_wavefunction,
		     (const Integer *basis,
		      const Integer *Nc, const Integer *No, const Integer *Nv,
		      const char *trans, const double *C, const Integer *ldC,
		      const double *E),
		     New_wavefunction,
		     (*basis, *Nc, *No, *Nv,
		      !(*trans == 'T' || *trans == 't'), C, *ldC, E))

    FORTRAN_FUNCTION(Integer new_integral_screening,
		     (const Integer *basis, const Integer *M,
		      const double *K, const double *cutoff),
		     New_integral_screening, (*basis, *M, K, *cutoff))

    FORTRAN_SUBROUTINE(hf_fock, (const Integer *basis, const double *D, double *F,
				 const double *cscale, const double *xscale,
				 const double *screen, const double *tolerance),
		       HF_fock, (*basis, D, F, *cscale, *xscale,
				 screen, *tolerance))

    FORTRAN_FUNCTION(double mp2_energy, (const Integer *wavefunction,
					 const Integer *screening),
		     MP2_energy, (*wavefunction, *screening))


    FORTRAN_SUBROUTINE(new_array,
		       (const Integer *N, const Integer *dims,
			const char *name, const Integer L),
		       bindings::new_array,
		       (*N, dims, string(name, L)))

    FORTRAN_SUBROUTINE(array_put,
		       (const double *buffer,
			const Integer *start, const Integer *stop,
			const char *name, const Integer L),
		       bindings::array_put,
		       (buffer, start, stop, string(name, L)))

    FORTRAN_SUBROUTINE(array_get,
		       (double *buffer,
			const Integer *start, const Integer *stop,
			const char *name, const Integer L),
		       bindings::array_get,
		       (buffer, start, stop, string(name, L)))

    FORTRAN_FUNCTION(double array_at,
		     (const Integer *index,const char *name, const Integer L),
		     bindings::array_at, (index, string(name, L)))
    
    FORTRAN_SUBROUTINE(cchem_cc_sd_vvvv,
		       (const Integer *wf),
		       CC_sd_vvvv, (*wf))

    FORTRAN_SUBROUTINE(cchem_cc_triples,
		       (const Integer *no, const Integer *nu,
			const double *eh, const double *ep,
			double *E),
		       CC_triples, (*no, *nu, eh, ep, E))

    // FORTRAN_SUBROUTINE(cchem_cc_t3_ijk,
    // 		       (const Integer *no, const Integer *nu,
    // 			const double *eh, const double *ep,
    // 			const Integer *i, const Integer *j, const Integer *k,
    // 			double *C, double *u1, double *T3),
    // 		       CC_T3_ijk, (*no, *nu, eh, ep, *i, *j, *k, C, u1, T3))


    FORTRAN_SUBROUTINE(dft_fock,
		       (const Integer *wavefunction,
			const double *rcut, const double *ccut,
			const Integer *N,
			const double *dr, const Integer *ldr,
			const double *w, const Integer *ldw,
			double *F,
			double *Xa, double *Xg, double *Ec, double *Rho,
			const char *functional, const int strlen),
		       DFT_fock,
		       (*wavefunction, *rcut, *ccut, *N,
			dr, *ldr, w, *ldw,
			F, Xa, Xg, Ec, Rho,
			string(functional, strlen)))


    // FORTRAN_FUNCTION(Integer eigen_jacoby, (const Integer *N,
    // 					    double *A, double *w,
    // 					    double *V, const Integer *ldV),
    // 		     Eigen_jacoby, (*N, A, w, V, *ldV))

}
