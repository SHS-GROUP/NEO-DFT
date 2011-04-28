#include "mmath.h"

/** @see mmath.h */
void mfill(double alpha, double *A, int ALda, int Ai, int Aj, int m, int n) {
  int i, j;
  int offsetA = Aj * ALda + Ai;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      A[offsetA + i] = alpha;
    }
    offsetA += ALda;
  }
}

/** @see mmath.h */
void mscale(double alpha, double *A, int ALda, int Ai, int Aj, int m, int n) {
  int i, j;
  int offsetA = Aj * ALda + Ai;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      A[offsetA + i] *= alpha;
    }
    offsetA += ALda;
  }
}

/** @see mmath.h */
void mmin(double *A, int ALda, int Ai, int Aj, int m, int n, double *alpha, int *index) {
  int i, j;
  int offsetA = Aj * ALda + Ai;

  /* initial element */
  *alpha = A[offsetA];
  index[0] = 0;
  index[1] = 0;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      if (A[offsetA + i] < *alpha) {
	*alpha = A[offsetA + i];
	index[0] = i;
	index[1] = j;
      }
    }
    offsetA += ALda;
  }
  index[0] += Ai;
  index[1] += Aj;
}

/** @see mmath.h */
void mmax(double *A, int ALda, int Ai, int Aj, int m, int n, double *alpha, int *index) {
  int i, j;
  int offsetA = Aj * ALda + Ai;

  /* initial element */
  *alpha = A[offsetA];
  index[0] = 0;
  index[1] = 0;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      if (A[offsetA + i] > *alpha) {
	*alpha = A[offsetA + i];
	index[0] = i;
	index[1] = j;
      }
    }
    offsetA += ALda;
  }
  index[0] += Ai;
  index[1] += Aj;
}

/** @see mmath.h */
void mmadd(double alpha, double *A, int ALda, int Ai, int Aj,
	   double beta, double *B, int BLda, int Bi, int Bj,
	   double *C, int CLda, int Ci, int Cj, int m, int n) {
  int i, j;
  int offsetA = Aj * ALda + Ai;
  int offsetB = Bj * BLda + Bi;
  int offsetC = Cj * CLda + Ci;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      C[offsetC + i] = alpha * A[offsetA + i] + beta * B[offsetB + i];
    }
    offsetA += ALda;
    offsetB += BLda;
    offsetC += CLda;
  }
}

/** @see mmath.h */
void mmdot(double *A, int ALda, int Ai, int Aj, double *B, int BLda, int Bi, int Bj,
	   double *x, int m, int n) {
  int i, j;
  int offsetA = Aj * ALda + Ai;
  int offsetB = Bj * BLda + Bi;

  for (j = 0; j < n; ++j) {
    for (i = 0; i < m; ++i) {
      x[i] += A[offsetA + i] * B[offsetB + i];
    }
    offsetA += ALda;
    offsetB += BLda;
  }
}
