#ifndef MMATH_H
#define MMATH_H

/**
   @file
   @brief Matrix math functions.
   @note All functions assume row-wise matrix layout as in Fortran
   and indices starting at 0.
   @warning The functions are not optimized.
*/

/**
   @brief Fill matrix A with a scalar.  A = alpha.
   @param[in] alpha Scalar.
   @param[out] A Matrix A
   @param[in] ALda Matrix A leading dimension.
   @param[in] Ai Matrix A first row.
   @param[in] Aj Matrix A first column.
   @param[in] m Row count.
   @param[in] n Column count.
 */
void mfill(double alpha, double *A, int ALda, int Ai, int Aj, int m, int n);

/**
   @brief Scale matrix A.  A = alpha*A.
   @param[in] alpha Scale factor.
   @param[in,out] A Matrix A.
   @param[in] ALda Matrix A leading dimension.
   @param[in] Ai Matrix A first row.
   @param[in] Aj Matrix A first column.
   @param[in] m Row count.
   @param[in] n Column count.
*/
void mscale(double alpha, double *A, int ALda, int Ai, int Aj, int m, int n);

/**
   @brief Find minimum value and its indices of matrix A.
   @note Out of two equal values, the first along the row takes precedence.
   @param[in] A Matrix A.
   @param[in] ALda Matrix A leading dimension.
   @param[in] Ai Matrix A first.
   @param[in] Aj Matrix A first column.
   @param[in] m Row count.
   @param[in] n Column count.
   @param[out] alpha Minimum value.
   @param[out] index 2-dimensional index array.
 */
void mmin(double *A, int ALda, int Ai, int Aj, int m, int n, double *alpha, int *index);

/**
   @brief Find maximum value and its indices of matrix A.
   @note Out of two equal values, the first along the row takes precedence.
   @param[in] A Matrix A.
   @param[in] ALda Matrix A leading dimension.
   @param[in] Ai Matrix A first row.
   @param[in] Aj Matrix A first column.
   @param[in] m Row count.
   @param[in] n Column count.
   @param[out] alpha Maximum value.
   @param[out] index 2-dimensional index array.
 */
void mmax(double *A, int ALda, int Ai, int Aj, int m, int n, double *alpha, int *index);

/**
   @brief Add two scaled matrices.  C = alpha*A + beta*B.
   @param[in] alpha Matrix A scale factor.
   @param[in] A Matrix A.
   @param[in] ALda Matrix A leading dimension.
   @param[in] Ai Matrix A first row.
   @param[in] Aj Matrix A first column.
   @param[in] beta Matrix B scale factor.
   @param[in] B Matrix B.
   @param[in] BLda Matrix B leading dimension.
   @param[in] Bi Matrix B first row.
   @param[in] Bj Matrix B first column.
   @param[out] C Matrix C.
   @param[in] CLda Matrix C leading dimension.
   @param[in] Ci Matrix C first row.
   @param[in] Cj Matrix C first column.
   @param[in] m Row count.
   @param[in] n Column count.
 */
void mmadd(double alpha, double *A, int ALda, int Ai, int Aj,
	   double beta, double *B, int BLda, int Bi, int Bj,
           double *C, int CLda, int Ci, int Cj, int m, int n);

/**
   @brief Dot product of two matrices.  x(i) = A(i:).B(i:).
   @note Dot product is computed from elements in the same row.
   @param[in] A Matrix A.
   @param[in] ALda Matrix A leading dimension.
   @param[in] Ai Matrix A first row.
   @param[in] Aj Matrix A first column.
   @param[in] B Matrix B.
   @param[in] BLda Matrix B leading dimension.
   @param[in] Bi Matrix B first row.
   @param[in] Bj Matrix B first column.
   @param[out] x The array of length m to store dot products.
   @param[in] m Row count.
   @param[in] n Column count.
*/
void mmdot(double *A, int ALda, int Ai, int Aj, double *B, int BLda, int Bi, int Bj,
           double *x, int m, int n);

#endif /* MMATH_H */

