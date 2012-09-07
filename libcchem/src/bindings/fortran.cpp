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

    FORTRAN_SUBROUTINE(cchem_initialize, (), CChem_initialize, ());

    FORTRAN_SUBROUTINE(cchem_runtime_set_int,
    		       (const char *key, const integer_t *value, integer_t L),
    		       CChem_runtime_set_int,
    		       (string(key, L), *value))

    FORTRAN_SUBROUTINE(cchem_runtime_set_double,
    		       (const char *key, const double *value, integer_t L),
    		       CChem_runtime_set_double,
    		       (string(key, L), *value))
					 

    FORTRAN_SUBROUTINE(delete_molecule, (const integer_t *key),
		       Delete_molecule, (*key))

    FORTRAN_SUBROUTINE(delete_wavefunction, (const integer_t *key),
		       Delete_wavefunction, (*key))

    FORTRAN_SUBROUTINE(delete_basis, (const integer_t *key),
		       Delete_basis, (*key))

    FORTRAN_SUBROUTINE(basis_sort, (const integer_t *key),
		       Basis_sort, (*key))

    FORTRAN_SUBROUTINE(basis_recenter, (const integer_t *basis,
					const double *x, const integer_t *incx,
					const double *y, const integer_t *incy,
					const double *z, const integer_t *incz),
		       Basis_recenter, (*basis, x, *incx, y, *incy, z, *incz))


    FORTRAN_FUNCTION(integer_t new_wavefunction, (const integer_t *basis),
		     New_wavefunction, (*basis))

    FORTRAN_SUBROUTINE(wavefunction_set_occupied,
		       (const integer_t *wf,
			const integer_t *start, const integer_t *stop,
			const integer_t *core),
		       Wavefunction_set_occupied,
		       (*wf, *start, *stop, *core))

    FORTRAN_SUBROUTINE(wavefunction_set_virtuals,
		       (const integer_t *wf,
			const integer_t *start, const integer_t *stop),
		       Wavefunction_set_virtuals,
		       (*wf, *start, *stop))

    FORTRAN_SUBROUTINE(wavefunction_set_c,
		       (const integer_t *wf,
			const char *trans,
			const integer_t *m, const integer_t *n,
			const double *C, const integer_t *ld),
		       Wavefunction_set_C, (*wf, *trans, *m, *n, C, *ld))

    FORTRAN_SUBROUTINE(wavefunction_set_f,
		       (const integer_t *wf,
			const double *F, const integer_t *ldF),
		       Wavefunction_set_F, (*wf, F, *ldF))

    FORTRAN_SUBROUTINE(wavefunction_set_e,
		       (const integer_t *wf, const double *e),
		       Wavefunction_set_e, (*wf, e))

    FORTRAN_SUBROUTINE(wavefunction_normalize, (const integer_t *wf),
		       Wavefunction_normalize, (*wf))


    FORTRAN_FUNCTION(integer_t new_integrals_screening,
		     (const integer_t *basis, const integer_t *M,
		      const double *K, const double *cutoff),
		     New_integrals_screening, (*basis, *M, K, *cutoff))

    FORTRAN_FUNCTION(integer_t cchem_hf_fock,
		     (const integer_t *basis,
		      const double *D, double *F,
		      const double *cscale, const double *xscale,
		      const double *screen, const double *tolerance),
		     CChem_hf_fock, (*basis, D, F, *cscale, *xscale,
				     screen, *tolerance))

    FORTRAN_FUNCTION(integer_t cchem_mp2_energy,
		     (const integer_t *wavefunction, double *E),
		     CChem_mp2_energy, (*wavefunction, E))


    // FORTRAN_SUBROUTINE(new_array,
    // 		       (const integer_t *N, const integer_t *dims,
    // 			const char *name, const integer_t L),
    // 		       bindings::new_array,
    // 		       (*N, dims, string(name, L)))

    // FORTRAN_SUBROUTINE(array_put,
    // 		       (const double *buffer,
    // 			const integer_t *start, const integer_t *stop,
    // 			const char *name, const integer_t L),
    // 		       bindings::array_put,
    // 		       (buffer, start, stop, string(name, L)))

    // FORTRAN_SUBROUTINE(array_get,
    // 		       (double *buffer,
    // 			const integer_t *start, const integer_t *stop,
    // 			const char *name, const integer_t L),
    // 		       bindings::array_get,
    // 		       (buffer, start, stop, string(name, L)))

    // FORTRAN_FUNCTION(double array_at,
    // 		     (const integer_t *index,const char *name, const integer_t L),
    // 		     bindings::array_at, (index, string(name, L)))

    
    FORTRAN_FUNCTION(integer_t cchem_cc_energy,
		     (const integer_t *wf, double *E,
		      const char *method, const integer_t L),
		     CChem_cc_energy, (*wf, E, string(method, L)))


    // FORTRAN_SUBROUTINE(cc_integrals, (const integer_t *wf),
    // 		       CC_integrals, (*wf))

    // FORTRAN_SUBROUTINE(cc_sd_vvvv, (const integer_t *wf),
    // 		       CC_sd_vvvv, (*wf))

    // FORTRAN_SUBROUTINE(cc_sd_evaluate, (const integer_t *wf),
    // 		       CC_sd_evaluate, (*wf))

    // FORTRAN_SUBROUTINE(cc_triples,
    // 		       (const integer_t *no, const integer_t *nu,
    // 			const double *eh, const double *ep,
    // 			double *E),
    // 		       CC_triples, (*no, *nu, eh, ep, E))

    // FORTRAN_SUBROUTINE(cchem_cc_t3_ijk,
    // 		       (const integer_t *no, const integer_t *nu,
    // 			const double *eh, const double *ep,
    // 			const integer_t *i, const integer_t *j, const integer_t *k,
    // 			double *C, double *u1, double *T3),
    // 		       CC_T3_ijk, (*no, *nu, eh, ep, *i, *j, *k, C, u1, T3))


    FORTRAN_SUBROUTINE(cchem_dft_fock,
		       (const integer_t *wavefunction,
			const double *rcut, const double *ccut,
			const integer_t *N,
			const double *dr, const integer_t *ldr,
			const double *w, const integer_t *ldw,
			double *F,
			double *Xa, double *Xg, double *Ec, double *Rho,
			const char *functional, const int strlen),
		       CChem_dft_fock,
		       (*wavefunction, *rcut, *ccut, *N,
			dr, *ldr, w, *ldw,
			F, Xa, Xg, Ec, Rho,
			string(functional, strlen)))


    // FORTRAN_FUNCTION(integer_t eigen_jacoby, (const integer_t *N,
    // 					    double *A, double *w,
    // 					    double *V, const integer_t *ldV),
    // 		     Eigen_jacoby, (*N, A, w, V, *ldV))

}
