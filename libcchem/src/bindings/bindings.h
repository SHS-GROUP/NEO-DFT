#ifndef BINDINGS_BINDINGS_H
#define BINDINGS_BINDINGS_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

    void CXX_cout(char type, size_t size, const void *data);

    void CChem_initialize();

    void CChem_runtime_set_int(const char *key, int value);
    void CChem_runtime_set_double(const char *key, double value);

    void Delete_molecule(int key);

    void Delete_basis(int key);

    void Delete_wavefunction(int key);

    void Delete_integrals_screening(int key);

    int New_molecule(int size, int *Z, double *r3, double *Q);

    int New_basis(int molecule);

    void Basis_sort(int basis);

    void Basis_recenter(int basis, 
			const double* x, int incx,
			const double* y, int incy,
			const double* z, int incz);


    int New_wavefunction(int basis);

    void Wavefunction_set_virtuals(int wf, size_t start, size_t stop);

    void Wavefunction_set_occupied(int wf, size_t start, size_t stop, size_t core);

    void Wavefunction_set_C(int wf, char trans, size_t m, size_t n,
			    const double* C, size_t ld );

    void Wavefunction_set_F(int wf, const double* F, size_t ld);

    void Wavefunction_set_e(int wf, const double* e);

    void Wavefunction_normalize(int wf);

    int New_integrals_screening(int basis, size_t N,
			       const double* K, double cutoff);
    
    int CChem_hf_fock(int basis, const double *D, double *F,
		      double cscale, double xscale,
		      const double *screen, double tolerance);

    int CChem_mp2_energy(int wf, double *E);

    int CChem_cc_energy(int wf, double *E, const char *method);

    void CChem_dft_fock(int wavefunction, double rcut, double ccut, int N,
			const double *dr, int ldr, const double *w, int ldw,
			double *F,
			double *Xa, double *Xg, double *Ec, double *Rho,
			const char *functional);


#ifdef __cplusplus
}
#endif

#endif /* BINDINGS_BINDINGS_H */
