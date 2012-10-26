#ifndef GAMESS_H_
#define GAMESS_H_

#include "bindings/fortran.h"
#include <stdlib.h>

/**
   @brief MXATM parameter in GAMESS
*/
#define MXATM 2000

/**
   @brief MXGTOT parameter in GAMESS
*/
#define MXGTOT 20000

/**
   @brief MXSH parameter in GAMESS
*/
#define MXSH 5000

#define INFOA infoa_
#define NSHEL nshel_
#define SHLNOS shlnos_

extern "C" {

    /**
       @brief C interface to GAMESS INFOA common block
    */
    extern struct {
	integer_t nat; //number of atom
	integer_t ich; //integer charge
	integer_t mul; //multiplicity
	integer_t num; //something related to symmetry
	integer_t nqmt; //something related to symmetry
	integer_t ne; //number of electrons
	integer_t na; //number of alpha electrons
	integer_t nb; //number of beta electrons
	double zan[MXATM]; //atomic charges/numbers
	double c[MXATM][3]; //atomic centers
	integer_t ian[MXATM]; //atomic numbers
    }INFOA;

    /**
       @brief C interface to GAMESS NSHEL common block
    */
    extern struct {
	double ex[MXGTOT]; //Exponents
	double cs[MXGTOT]; //Coefficients
	double cp[MXGTOT]; //P. coefficients
	double cd[MXGTOT];//D. coefficients
	double cf[MXGTOT];//F. coefficients
	double cg[MXGTOT];//G. coefficients
	double ch[MXGTOT];//H. coefficients
	double ci[MXGTOT];//I. coefficients
	integer_t kstart[MXSH]; //kstart index (+1) to ex,cs,...
	integer_t katom[MXSH]; //Center + 1
	integer_t ktype[MXSH]; //L + 1
	integer_t kng[MXSH]; //K, number of primitive functions
	integer_t kloc[MXSH];//Shell functions location
	integer_t kmin[MXSH]; //kmin==1 && kmax==4 => type = SP
	integer_t kmax[MXSH];
	integer_t nshell; //Num shells
    } NSHEL;

    extern void FORTRAN_CALL(ddi_dlbreset)();
    extern void FORTRAN_CALL(ddi_sync)(integer_t*);

    extern void FORTRAN_CALL(ddi_dlbnext)(integer_t*);
    extern void FORTRAN_CALL(ddi_nproc)(integer_t*,integer_t*);
    extern void FORTRAN_CALL(ddi_smp_nproc)(integer_t*,integer_t*);

    extern void FORTRAN_CALL(ddi_gsumi)(integer_t*,void*,integer_t*);
    extern void FORTRAN_CALL(ddi_gsumf)(integer_t*,void*,integer_t*);
    extern void FORTRAN_CALL(ddi_bcast)
	(integer_t*,char*,void*,integer_t*,integer_t*);

    extern void FORTRAN_CALL(ddi_masters_gsumi)(integer_t*,void*,integer_t*);
    extern void FORTRAN_CALL(ddi_masters_gsumf)(integer_t*,void*,integer_t*);
    extern void FORTRAN_CALL(ddi_masters_bcast)
	(integer_t*,char*,void*,integer_t*,integer_t*);

}

#endif /* GAMESS_H_ */
