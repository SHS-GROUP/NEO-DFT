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
	Integer nat; //number of atom
	Integer ich; //integer charge
	Integer mul; //multiplicity
	Integer num; //something related to symmetry
	Integer nqmt; //something related to symmetry
	Integer ne; //number of electrons
	Integer na; //number of alpha electrons
	Integer nb; //number of beta electrons
	double zan[MXATM]; //atomic charges/numbers
	double c[MXATM][3]; //atomic centers
	Integer ian[MXATM]; //atomic numbers
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
	Integer kstart[MXSH]; //kstart index (+1) to ex,cs,...
	Integer katom[MXSH]; //Center + 1
	Integer ktype[MXSH]; //L + 1
	Integer kng[MXSH]; //K, number of primitive functions
	Integer kloc[MXSH];//Shell functions location
	Integer kmin[MXSH]; //kmin==1 && kmax==4 => type = SP
	Integer kmax[MXSH];
	Integer nshell; //Num shells
    } NSHEL;


#define DDI_DLBRESET ddi_dlbreset_
#define DDI_DLBNEXT ddi_dlbnext_
#define DDI_NPROC ddi_nproc_
#define DDI_SMP_NPROC ddi_smp_nproc_

    extern void DDI_DLBRESET();
    extern void DDI_DLBNEXT(Integer*);
    extern void DDI_NPROC(Integer*,Integer*);
    extern void DDI_SMP_NPROC(Integer*,Integer*);

}

#endif /* GAMESS_H_ */
