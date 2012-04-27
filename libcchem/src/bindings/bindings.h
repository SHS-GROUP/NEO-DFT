#ifndef BINDINGS_BINDINGS_H
#define BINDINGS_BINDINGS_H

#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

    void CXX_cout(char type, size_t size, const void *data);

    //void Runtime_set_int(const char *key, int value);

    void Delete_molecule(int key);

    void Delete_basis(int key);

    void Delete_wavefunction(int key);

    void Delete_integral_screening(int key);

    int New_molecule(int size, int *Z, double *r3, double *Q);

    int New_basis(int molecule);

    void Basis_sort(int basis);

    void Basis_recenter(int basis, 
			const double* x, int incx,
			const double* y, int incy,
			const double* z, int incz);

    int New_wavefunction(int basis,
			 size_t core, size_t active, size_t virtuals,
			 bool axm, const double* C, size_t ldC,
			 const double* E);

    int New_integral_screening(int basis, size_t N,
			       const double* K, double cutoff);

    void New_array(size_t N, const size_t* dims, const char* name);

    void Array_put(const double *buffer, const size_t *start, const size_t *stop,
		   const char* name);
    
    void HF_fock(const int basis, const double *D, double *F,
		 double cscale, double xscale,
		 const double *screen, double tolerance);

    double MP2_energy(int wf, int screening);

    void CC_sd_vvvv(int wf);

    void CC_triples(size_t no, size_t nv,
		    const double *eh, const double *ep,
		    double* C);

    void DFT_fock(int wavefunction, double rcut, double ccut, int N,
		  const double *dr, int ldr, const double *w, int ldw,
		  double *F,
		  double *Xa, double *Xg, double *Ec, double *Rho,
		  const char *functional);


#ifdef __cplusplus
}
#endif

#endif /* BINDINGS_BINDINGS_H */
