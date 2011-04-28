#include "ddi_base.h"

/* Dot product of two matrices
   do i=1,m
     do j=1,n
       s(i) = s(i) + A(i,j)*B(i,j)
     enddo
   enddo
*/  

void DDI_ARR_dot_(DDI_Patch *dAPatch, DDI_Patch *dBPatch, double *x) {
  int i, m;
  DDI_Patch dASubPatch[MAX_SMP_PROCS];
  DDI_Patch dBSubPatch[MAX_SMP_PROCS];
  int dASubPatchRank[MAX_SMP_PROCS], dANSubPatch;
  int dBSubPatchRank[MAX_SMP_PROCS], dBNSubPatch;
  const DDI_Comm *comm = (const DDI_Comm *)Comm_find(DDI_WORKING_COMM);
  int rank = comm->me, node = comm->my;

  DDI_Subpatch_node(dAPatch->handle, dAPatch, &dANSubPatch, dASubPatchRank, dASubPatch, node);
  DDI_Subpatch_node(dBPatch->handle, dBPatch, &dBNSubPatch, dBSubPatchRank, dBSubPatch, node);

  /* number of A and B segments must be the same */
  if (dANSubPatch != dBNSubPatch) {
    Fatal_error(1000);
  }

  m = dAPatch->ihi - dAPatch->ilo + 1;
  /* zero out dot vector */
  for (i = 0; i < m; ++i) {
    *((double*)x+i) = (double)0;
  }

  for (i = 0; i < dANSubPatch; ++i) {
    if (rank == dASubPatchRank[i]) {
      /* ensure segments are aligned on the processor */
      if (dBSubPatchRank[i] != dASubPatchRank[i]) {
	Fatal_error(1000);
      }
#ifdef DDI_ARR_DMA
      DDI_ARR_dot_local(&dASubPatch[i], &dBSubPatch[i], x);
#else
      dASubPatch[i].size = m * sizeof(double);
      DDI_ARR_dot_remote(&dASubPatch[i], &dBSubPatch[i], x, rank);
#endif
    }
  }
  /* sum partial dot products from other processors */
  Comm_gsum(x, m, 0, comm);
}

void DDI_ARR_dot_local(DDI_Patch *dAPatch, DDI_Patch *dBPatch, double *x) {
  DDA_Index *Index = gv(dda_index);
  int dA = dAPatch->handle, dB = dBPatch->handle;
  double *dALocal, *dBLocal;

  int dALdaLocal, dAiLocal, dAjLocal, dAmLocal, dAnLocal;
  int dBLdaLocal, dBiLocal, dBjLocal, dBmLocal, dBnLocal;

  /* A segment dimensions */
  dALdaLocal = Index[dA].ihi - Index[dA].ilo + 1;
  dAiLocal = dAPatch->ilo - Index[dA].ilo;
  dAjLocal = dAPatch->jlo - Index[dA].jlo;
  dAmLocal = dAPatch->ihi - dAPatch->ilo + 1;
  dAnLocal = dAPatch->jhi - dAPatch->jlo + 1;

  /* B segment dimensions */
  dBLdaLocal = Index[dB].ihi - Index[dB].ilo + 1;
  dBiLocal = dBPatch->ilo - Index[dB].ilo;
  dBjLocal = dBPatch->jlo - Index[dB].jlo;
  dBmLocal = dBPatch->ihi - dBPatch->ilo + 1;
  dBnLocal = dBPatch->jhi - dBPatch->jlo + 1;
  
# if defined USE_SYSV
    if(USING_DATA_SERVERS()) {
      DDI_Fence_check(dA);
      DDI_Fence_check(dB);
    }
# endif

    DDI_Acquire(Index, dA, DDI_READ_ACCESS, (void**)&dALocal);
    if (dB == dA) dBLocal = dALocal;
    else DDI_Acquire(Index, dB, DDI_READ_ACCESS, (void**)&dBLocal);

    /* Matrix segments have to be the same. */
    if ((dAmLocal != dBmLocal) || (dAnLocal != dBnLocal)) {
      Fatal_error(1000);
    }
    /* partial dot product */
    mmdot(dALocal, dALdaLocal, dAiLocal, dAjLocal, dBLocal, dBLdaLocal, dBiLocal, dBjLocal,
	  x, dAmLocal, dAnLocal); 
    
    DDI_Release(Index, dA, DDI_READ_ACCESS);
    if (dB != dA) DDI_Release(Index, dB, DDI_READ_ACCESS);
}

void DDI_ARR_dot_remote(DDI_Patch *dAPatch, DDI_Patch *dBPatch, double *x, int rank) {
  const DDI_Comm *comm = (const DDI_Comm *)Comm_find(DDI_WORKING_COMM);

  /* send A and B segments and receive partial dot product */
  DDI_Send_request(dAPatch, &rank, NULL);
  Comm_send(dBPatch, sizeof(DDI_Patch), rank, comm);
  Comm_recv(x, dAPatch->size, rank, comm);
}

void DDI_ARR_dot_server(DDI_Patch *dAPatch, int rank) {
  DDI_Patch dBPatch;
  const DDI_Comm *comm = (const DDI_Comm *)Comm_find(DDI_WORKING_COMM);

  /* receive B segment, compute, and send partial dot product */

/*        original coding, with run-time storage allocation in next line  */
/*
 *   char x[dAPatch->size];
 *   Comm_recv(&dBPatch, sizeof(DDI_Patch), rank, comm);
 *   DDI_ARR_dot_local(dAPatch, &dBPatch, (double*)x);
 *   Comm_send(x, dAPatch->size, rank, comm);
 */

/*        replaced by Mike's re-coding doing explicit alloc/dealloc.   */

  int *x;
  size_t lenx;
  lenx = dAPatch->size;
  x = malloc(lenx);

  Comm_recv(&dBPatch, sizeof(DDI_Patch), rank, comm);
  DDI_ARR_dot_local(dAPatch, &dBPatch, (double*)x);
  Comm_send(x, dAPatch->size, rank, comm);

  (void) free((void*) x);
}

