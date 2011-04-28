#ifndef DDI_ARR_H
#define DDI_ARR_H

/**
   @file
   @brief Math operations on distributed arrays (Not completely documented yet).
   @note Distributed arrays use row-wise layout as in Fortran,
   indices start at 0, and are distributed columns-wise, i.e. broken
   down along column lines.
   @note The scope of the distributed array is the scope in which that array was created.
*/

/**
   @brief Fill the distributed array A segment with zero.
   @note This operation is collective in the scope of the distributed array A.
   @param dA Distributed array A handle.
   @param dAi Distributed array A segment first row.
   @param dAi2 Distributed array A segment last row.
   @param dAj Distributed array A segment first column.
   @param dAj2 Distributed array A segment last column.
*/
void DDI_ARR_zero(int dA, int dAi, int dAi2, int dAj, int dAj2);

/**
   @brief Fill the distributed array A segment with a scalar.
   @note This operation is collective in the scope of the distributed array A.
   @param dA Distributed array A handle.
   @param dAi Distributed array A segment first row.
   @param dAi2 Distributed array A segment last row.
   @param dAj Distributed array A segment first column.
   @param dAj2 Distributed array A segment last column.
   @param alpha Scalar.
*/
void DDI_ARR_fill(int dA, int dAi, int dAi2, int dAj, int dAj2, double alpha);

/**
   @brief Scale the distributed array A segment by a scalar.
   @note This operation is collective in the scope of the distributed array A.
   @param dA Distributed array A handle.
   @param dAi Distributed array A segment first row.
   @param dAi2 Distributed array A segment last row.
   @param dAj Distributed array A segment first column.
   @param dAj2 Distributed array A segment last column.
   @param alpha Scalar.
*/
void DDI_ARR_scale(int dA, int dAi, int dAi2, int dAj, int dAj2, double alpha);

/**
   @brief Find minimum value and its index in the distributed array A segment.
   @details The value and its index are broadcast in the scope of the distributed array A.
   @note This operation is collective in the scope of the distributed array A.
   @param dA Distributed array A handle.
   @param dAi Distributed array A segment first row.
   @param dAi2 Distributed array A segment last row.
   @param dAj Distributed array A segment first column.
   @param dAj2 Distributed array A segment last column.
   @param alpha Pointer to store minimum value.
   @param index Pointer to store 2-dimension index.
*/
void DDI_ARR_min(int dA, int dAi, int dAi2, int dAj, int dAj2, double *alpha, int *index);

/**
   @brief Find maximum value and its index in the distributed array A segment.
   @details The value and its index are broadcast in the scope of the distributed array A.
   @note This operation is collective in the scope of the distributed array A.
   @param dA Distributed array A handle.
   @param dAi Distributed array A segment first row.
   @param dAi2 Distributed array A segment last row.
   @param dAj Distributed array A segment first column.
   @param dAj2 Distributed array A segment last column.
   @param alpha Pointer to store maximum value.
   @param index Pointer to store 2-dimension index.
*/
void DDI_ARR_max(int dA, int dAi, int dAi2, int dAj, int dAj2, double *alpha, int *index);

/**
   @brief Store dot product of distributed arrays A and B segments in local array.
   x(i) = A(i:).B(i:).
   @note Distributed arrays A and B must have the same scope.
   @note Distributed array A and B segments must have the same dimensions
   and the same processor locality.
   @details The dot product is broadcast in the scope of the distributed array A.
   @note This operation is collective in the scope of the distributed arrays A and B.
   @param dA Distributed array A handle.
   @param dAi Distributed array A segment first row.
   @param dAi2 Distributed array A segment last row.
   @param dAj Distributed array A segment first column.
   @param dAj2 Distributed array A segment last column.
   @param dB Distributed array B handle.
   @param dBi Distributed array B segment first row.
   @param dBi2 Distributed array B segment last row.
   @param dBj Distributed array B segment first column.
   @param dBj2 Distributed array B segment last column.
   @param x Pointer to local array of size dAi2 - dAi.
*/
void DDI_ARR_dot(int dA, int dAi, int dAi2, int dAj, int dAj2,
	     int dB, int dBi, int dBi2, int dBj, int dBj2, double *x);

/**
   @brief Add scaled distributed arrays A and B segments and store result in array C segment.
   C = alpha*A + beta*B.
   @note Distributed arrays A, B, and C must have the same scope.
   @note Distributed array A,B, and C segments must have the same dimensions
   and the same processor locality.
   @note This operation is collective in the scope of the distributed arrays A, B, and C.
   @param dA Distributed array A handle.
   @param dAi Distributed array A segment first row.
   @param dAi2 Distributed array A segment last row.
   @param dAj Distributed array A segment first column.
   @param dAj2 Distributed array A segment last column.
   @param alpha Distributed array A segment scale factor.
   @param dB Distributed array B handle.
   @param dBi Distributed array B segment first row.
   @param dBi2 Distributed array B segment last row.
   @param dBj Distributed array B segment first column.
   @param dBj2 Distributed array B segment last column.
   @param beta Distributed array B segment scale factor.
   @param dC Distributed array C handle.
   @param dCi Distributed array C segment first row.
   @param dCi2 Distributed array C segment last row.
   @param dCj Distributed array C segment first column.
   @param dCj2 Distributed array C segment last column.
*/ 
void DDI_ARR_add(int dA, int dAi, int dAi2, int dAj, int dAj2, double alpha,
	     int dB, int dBi, int dBi2, int dBj, int dBj2, double beta,
	     int dC, int dCi, int dCi2, int dCj, int dCj2);

/**
   @brief One-sided accumulate operation.
   @details Scale and add local buffer to distributed array segment.
   @note Local buffer is unchanged by this operation.
   @param dA Distributed array A handle.
   @param dAi Distributed array A segment first row.
   @param dAi2 Distributed array A segment last row.
   @param dAj Distributed array A segment first column.
   @param dAj2 Distributed array A segment last column.
   @param alpha Local buffer scale factor.
   @param buf Pointer to local buffer.
*/
void DDI_ARR_acc(int dA, int dAi, int dAi2, int dAj, int dAj2, double alpha, double *buf);

#endif /* DDI_ARR_H */
#ifndef DDI_ARR_BASE_H
#define DDI_ARR_BASE_H

#if defined USE_SYSV || DDI_MPI2 || defined DDI_ARMCI
#define DDI_ARR_DMA
#endif

/* DDI_ARR operations */
#define DDI_ARR_ZERO              65
#define DDI_ARR_FILL              66
#define DDI_ARR_SCALE             67
#define DDI_ARR_MIN               68
#define DDI_ARR_MAX               69
#define DDI_ARR_DOT               70
#define DDI_ARR_ADD               71
#define DDI_ARR_ACC               72

typedef struct {
  double alpha;
  int index[2];
} DDI_ARR_Element;

typedef struct {
  double alpha;
  DDI_Patch patch;
} DDI_ARR_Patch;

/* Scalar DDI_ARR operations */
void DDI_ARR_scalar_(DDI_Patch *dAPatch, double alpha);
void DDI_ARR_scalar_local(DDI_Patch *dAPatch, double alpha);
void DDI_ARR_scalar_remote(DDI_Patch *dAPatch, double alpha, int rank);
void DDI_ARR_scalar_server(DDI_Patch *dAPatch, int rank);

/* Select DDI_ARR operations */
void DDI_ARR_select_(DDI_Patch *dAPatch, double *alpha, int *index);
void DDI_ARR_select_local(DDI_Patch *dAPatch, DDI_ARR_Element* element);
void DDI_ARR_select_remote(DDI_Patch *dAPatch, DDI_ARR_Element* element, int rank);
void DDI_ARR_select_server(DDI_Patch *dAPatch, int rank);
void DDI_ARR_Element_select(DDI_ARR_Element *a, DDI_ARR_Element *b, int op);

/* C = alpha*A + beta*B operations */
void DDI_ARR_add_(DDI_Patch *dAPatch, double alpha, DDI_Patch *dBPatch, double beta, DDI_Patch *dCPatch);
void DDI_ARR_add_local(DDI_Patch *dAPatch, double alpha, DDI_Patch *dBPatch, double beta, DDI_Patch *dCPatch);
void DDI_ARR_add_remote(DDI_Patch *dAPatch, double alpha, DDI_Patch *dBPatch, double beta, DDI_Patch *dCPatch, int rank);
void DDI_ARR_add_server(DDI_Patch *request, int rank);

/* x(i) = A(i:).B(i:) operations */
void DDI_ARR_dot_(DDI_Patch *dAPatch, DDI_Patch *dBPatch, double *x);
void DDI_ARR_dot_local(DDI_Patch *dAPatch, DDI_Patch *dBPatch, double *x);
void DDI_ARR_dot_remote(DDI_Patch *dAPatch, DDI_Patch *dBPatch, double *x, int rank);
void DDI_ARR_dot_server(DDI_Patch *dAPatch, int rank);

/* Accumulate DDI_ARR operations */
void DDI_ARR_acc_(DDI_Patch *dAPatch, double alpha, double *buf);
/* DDI_ARR_acc uses DDI_Acc */
/*
void DDI_ARR_acc_local(DDI_Patch *dAPatch, double alpha, double *buf);
void DDI_ARR_acc_remote(DDI_Patch *dAPatch, double alpha, double *buf, int rank);
*/
void DDI_ARR_acc_server(DDI_Patch *dAPatch, int rank);

#endif /* DDI_ARR_H */
