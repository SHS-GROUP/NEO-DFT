#include "ddi_base.h"

/** Fill array with zeroes */
void DDI_ARR_zero(int dA, int dAi, int dAi2, int dAj, int dAj2) {
  DDI_Patch dAPatch;

  dAPatch.handle = dA;
  #if defined WINTEL
  dAPatch.oper = DDI_ARR_ZERO_OP;
  #else
  dAPatch.oper = DDI_ARR_ZERO;
  #endif
  dAPatch.ilo = dAi;
  dAPatch.ihi = dAi2;
  dAPatch.jlo = dAj;
  dAPatch.jhi = dAj2;
  
  DDI_ARR_scalar_(&dAPatch, (double)0);
}

/** Fill array with alpha */
void DDI_ARR_fill(int dA, int dAi, int dAi2, int dAj, int dAj2, double alpha) {
  DDI_Patch dAPatch;

  dAPatch.handle = dA;
  #if defined WINTEL
  dAPatch.oper = DDI_ARR_FILL_OP;
  #else
  dAPatch.oper = DDI_ARR_FILL;
  #endif
  dAPatch.ilo = dAi;
  dAPatch.ihi = dAi2;
  dAPatch.jlo = dAj;
  dAPatch.jhi = dAj2;
  
  DDI_ARR_scalar_(&dAPatch, alpha);
}

/** Scale array by alpha */
void DDI_ARR_scale(int dA, int dAi, int dAi2, int dAj, int dAj2, double alpha) {
  DDI_Patch dAPatch;
  
  dAPatch.handle = dA;
  #if defined WINTEL
  dAPatch.oper = DDI_ARR_SCALE_OP;
  #else
  dAPatch.oper = DDI_ARR_SCALE;
  #endif
  dAPatch.ilo = dAi;
  dAPatch.ihi = dAi2;
  dAPatch.jlo = dAj;
  dAPatch.jhi = dAj2;
  
  DDI_ARR_scalar_(&dAPatch, alpha);
}

/** Find minimum value and its index in array */
void DDI_ARR_min(int dA, int dAi, int dAi2, int dAj, int dAj2, double *alpha, int *index) {
  DDI_Patch dAPatch;

  dAPatch.handle = dA;
  #if defined WINTEL
  dAPatch.oper = DDI_ARR_MIN_OP;
  #else
  dAPatch.oper = DDI_ARR_MIN;
  #endif
  dAPatch.ilo = dAi;
  dAPatch.ihi = dAi2;
  dAPatch.jlo = dAj;
  dAPatch.jhi = dAj2;

  DDI_ARR_select_(&dAPatch, alpha, index);
}

/** Find maximum value and its index in array */
void DDI_ARR_max(int dA, int dAi, int dAi2, int dAj, int dAj2, double *alpha, int *index) {
  DDI_Patch dAPatch;

  dAPatch.handle = dA;
  #if defined WINTEL
  dAPatch.oper = DDI_ARR_MAX_OP;
  #else
  dAPatch.oper = DDI_ARR_MAX;
  #endif
  dAPatch.ilo = dAi;
  dAPatch.ihi = dAi2;
  dAPatch.jlo = dAj;
  dAPatch.jhi = dAj2;

  DDI_ARR_select_(&dAPatch, alpha, index);
}

/** Compute dot product of two arrays */
void DDI_ARR_dot(int dA, int dAi, int dAi2, int dAj, int dAj2,
	     int dB, int dBi, int dBi2, int dBj, int dBj2, double *x) {
  DDI_Patch dAPatch, dBPatch;

  dAPatch.handle = dA;
  #if defined WINTEL
  dAPatch.oper = DDI_ARR_DOT_OP;
  #else
  dAPatch.oper = DDI_ARR_DOT;
  #endif
  dAPatch.ilo = dAi;
  dAPatch.ihi = dAi2;
  dAPatch.jlo = dAj;
  dAPatch.jhi = dAj2;

  dBPatch.handle = dB;
  #if defined WINTEL
  dBPatch.oper = DDI_ARR_DOT_OP;
  #else
  dBPatch.oper = DDI_ARR_DOT;
  #endif
  dBPatch.ilo = dBi;
  dBPatch.ihi = dBi2;
  dBPatch.jlo = dBj;
  dBPatch.jhi = dBj2;
  
  DDI_ARR_dot_(&dAPatch, &dBPatch, x);
}

/** Compute the sum of two arrays and store result in array C */ 
void DDI_ARR_add(int dA, int dAi, int dAi2, int dAj, int dAj2, double alpha,
	     int dB, int dBi, int dBi2, int dBj, int dBj2, double beta,
	     int dC, int dCi, int dCi2, int dCj, int dCj2) {
  DDI_Patch dAPatch, dBPatch, dCPatch;

  dAPatch.handle = dA;
  #if defined WINTEL
  dAPatch.oper = DDI_ARR_ADD_OP;
  #else
  dAPatch.oper = DDI_ARR_ADD;
  #endif
  dAPatch.ilo = dAi;
  dAPatch.ihi = dAi2;
  dAPatch.jlo = dAj;
  dAPatch.jhi = dAj2;

  dBPatch.handle = dB;
  #if defined WINTEL
  dBPatch.oper = DDI_ARR_ADD_OP;
  #else
  dBPatch.oper = DDI_ARR_ADD;
  #endif
  dBPatch.ilo = dBi;
  dBPatch.ihi = dBi2;
  dBPatch.jlo = dBj;
  dBPatch.jhi = dBj2;

  dCPatch.handle = dC;
  #if defined WINTEL
  dCPatch.oper = DDI_ARR_ADD_OP;
  #else
  dCPatch.oper = DDI_ARR_ADD;
  #endif
  dCPatch.ilo = dCi;
  dCPatch.ihi = dCi2;
  dCPatch.jlo = dCj;
  dCPatch.jhi = dCj2;

  DDI_ARR_add_(&dAPatch, alpha, &dBPatch, beta, &dCPatch);
}

/** Accumulate into array */
void DDI_ARR_acc(int dA, int dAi, int dAi2, int dAj, int dAj2, double alpha, double *buf) {
  DDI_Patch dAPatch;

  dAPatch.handle = dA;
  #if defined WINTEL
  dAPatch.oper = DDI_ARR_ACC_OP;
  #else
  dAPatch.oper = DDI_ARR_ACC;
  #endif
  dAPatch.ilo = dAi;
  dAPatch.ihi = dAi2;
  dAPatch.jlo = dAj;
  dAPatch.jhi = dAj2;

  DDI_ARR_acc_(&dAPatch, alpha, buf);
}
