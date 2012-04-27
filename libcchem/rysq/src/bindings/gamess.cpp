#include <rysq.h>
#include "bindings/gamess.h"
#include <stdint.h>
#include <algorithm>
#include <iostream>


extern "C" {
    void rysq_eri1_gamess(const Integer *flags, const Integer *a, const Integer *b,
			  const Integer *c, const Integer *d, double *I);

    void rysq_eri1_gamess_(const Integer *flags, const Integer *a, const Integer *b,
			   const Integer *c, const Integer *d, double *I);

#pragma weak rysq_eri1_gamess_ = rysq_eri1_gamess

}

// static Rysq_shell_t Rysq_new_shell(int index1) {
//     int index0 = index1 - 1;
//     int L = NSHEL.ktype[index0] - 1;
//     int K = NSHEL.kng[index0];
//     int k = NSHEL.kstart[index0]-1;

//     bool sp = (NSHEL.kmin[index0] == 1 && NSHEL.kmax[index0] == 4);
//     const double *C[] = {NSHEL.cs, NSHEL.cp, NSHEL.cd, NSHEL.cf, NSHEL.cg, NSHEL.ch, NSHEL.ci};

//     if (sp) return Rysq_new_shell2(K, &NSHEL.ex[k], &C[0][k], &C[1][k]);
//     return Rysq_new_shell(L, K, &NSHEL.ex[k], &C[L][k]);
// }

static inline double* Gamess_rysq_coord(int a) {
    return (const_cast<double*>(INFOA.c[NSHEL.katom[a-1]-1]));
}

/**
   @brief librysq GAMESS interface
   @warning All input indices are assumed to start at 1
   @param flags, 1 - normalize, 2 - float
   @param a first shell index, start at 1
   @param b second shell index, start at 1
   @param c third shell index, start at 1
   @param d fourth shell index, start at 1
   @param[out] Evaluated eri block 
*/
void rysq_eri1_gamess(const Integer *flags, const Integer *a, const Integer *b,
		      const Integer *c, const Integer *d, double *I) {  
	
    double cutoff = 1.0e-10;
    double scale = SHLNOS.qq4;

//     Rysq_shell_t A = Rysq_new_shell(*a);
//     Rysq_shell_t B = Rysq_new_shell(*b);
//     Rysq_shell_t C = Rysq_new_shell(*c);
//     Rysq_shell_t D = Rysq_new_shell(*d);

    double *ri = Gamess_rysq_coord(*a);
    double *rj = Gamess_rysq_coord(*b);
    double *rk = Gamess_rysq_coord(*c);
    double *rl = Gamess_rysq_coord(*d);
	
//     Rysq_eri_zero(A, B, C, D, I);

    //	Rysq_eri1(RYSQ_NORMALIZE | RYSQ_SINGLE, 1.0e-10, A, Gamess_rysq_coord(*a), B, Gamess_rysq_coord(*b),

    //     int L = A->L + B->L + C->L + D->L;
    //     int nc = A->nc*B->nc*C->nc*D->nc;

    // #ifdef RYSQ_ENABLE_CUDA
    // 	if( L < 2 ) {
    // 	    cuRysq_eri1(*flags, tol, A, ri, B, rj,
    // 			C, rk, D, rl, SHLNOS.qq4, I);
    // 	}

    // 	else if(A.type == RYSQ_TYPE_D && B.type == RYSQ_TYPE_S &&
    // 		C.type == RYSQ_TYPE_S && D.type == RYSQ_TYPE_S) {
    // 	    cuRysq_eri1_r2(*flags, tol, A, ri, B, rj,
    // 			   C, rk, D, rl, SHLNOS.qq4, I);


    // // 	    std::cout<<"*" << I[0]<< " " << I[1]<<  " " << I[2]<<  " " << I[5]<<  std::endl;
    // // 	std::fill(I, I+A.size*B.size*C.size*D.size, 0.0);
    // // 	    Rysq_eri1(*flags, tol, A, ri, B, rj,
    // // 		      C, rk, D, rl, SHLNOS.qq4, I);
    // // 	    std::cout<<" " << I[0]<< " " << I[1]<<  " " << I[2]<<  " " << I[5]<<  std::endl;

    // 	}

    // // 	else if(A.type == RYSQ_TYPE_P && B.type == RYSQ_TYPE_S &&
    // // 		C.type == RYSQ_TYPE_P && D.type == RYSQ_TYPE_S) {
    // // 	    cuRysq_eri1_r2(*flags, tol, A, ri, B, rj,
    // // 			   C, rk, D, rl, SHLNOS.qq4, I);
    // // 	}

    // // 	else if(A.type == RYSQ_TYPE_F && B.type == RYSQ_TYPE_S &&
    // // 		C.type == RYSQ_TYPE_S && D.type == RYSQ_TYPE_S) {
    // // 	    cuRysq_eri1_r2(*flags, tol, A, ri, B, rj,
    // // 			   C, rk, D, rl, SHLNOS.qq4, I);
    // // 	}

    // 	else {
    // 	    Rysq_eri1(*flags, tol, A, ri, B, rj,
    // 		      C, rk, D, rl, SHLNOS.qq4, I);
    // 	}
    // #else

//     Rysq_eri1(*flags, cutoff, A, ri, B, rj, C, rk, D, rl, scale, I);

    //#endif
	
    ERIOUT.iout = *d;
    ERIOUT.jout = *c;
    ERIOUT.kout = *b;
    ERIOUT.lout = *a;

//     ERIOUT.incl = 1;
//     ERIOUT.inck = ERIOUT.incl*Rysq_shell_size(A);
//     ERIOUT.incj = ERIOUT.inck*Rysq_shell_size(B);
//     ERIOUT.inci = ERIOUT.incj*Rysq_shell_size(C);

//     Rysq_delete_shell(A);
//     Rysq_delete_shell(B);
//     Rysq_delete_shell(C);
//     Rysq_delete_shell(D);

}
    
