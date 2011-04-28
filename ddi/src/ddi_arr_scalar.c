#include "ddi_base.h"

void DDI_ARR_scalar_(DDI_Patch *dAPatch, double alpha) {
  int i;
  int dANSubPatch, dASubPatchRank[MAX_SMP_PROCS];
  DDI_Patch dASubPatch[MAX_SMP_PROCS];
  const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
  int rank = comm->me;
  int node = comm->my;

  DDI_Subpatch_node(dAPatch->handle, dAPatch, &dANSubPatch, dASubPatchRank, dASubPatch, node);

  for (i = 0; i < dANSubPatch; ++i) {
    if (dASubPatchRank[i] == rank) {
#ifdef DDI_ARR_DMA
      DDI_ARR_scalar_local(&dASubPatch[i], alpha);
#else
      DDI_ARR_scalar_remote(&dASubPatch[i], alpha, rank);
#endif
    }
  }
  Comm_sync(1000, comm);
}

void DDI_ARR_scalar_local(DDI_Patch *dAPatch, double alpha) {
  DDA_Index *Index = gv(dda_index);
  int dA = dAPatch->handle;
  double *dALocal;
  int op = dAPatch->oper;

  int dALdaLocal, dAiLocal, dAjLocal, dAmLocal, dAnLocal;

  /* A segment dimensions */
  dALdaLocal = Index[dA].ihi - Index[dA].ilo + 1;
  dAiLocal = dAPatch->ilo - Index[dA].ilo;
  dAjLocal = dAPatch->jlo - Index[dA].jlo;
  dAmLocal = dAPatch->ihi - dAPatch->ilo + 1;
  dAnLocal = dAPatch->jhi - dAPatch->jlo + 1;

# if defined USE_SYSV
  if(USING_DATA_SERVERS()) DDI_Fence_check(dA);
# endif

  DDI_Acquire(Index, dA, DDI_READ_ACCESS, (void**)&dALocal);

  switch(op) {
  case DDI_ARR_ZERO:
  case DDI_ARR_FILL:
    mfill(alpha, dALocal, dALdaLocal, dAiLocal, dAjLocal, dAmLocal, dAnLocal);
    break;
  case DDI_ARR_SCALE:
    mscale(alpha, dALocal, dALdaLocal, dAiLocal, dAjLocal, dAmLocal, dAnLocal);
    break;
  } /* switch */

  DDI_Release(Index, dA, DDI_WRITE_ACCESS);
}

void DDI_ARR_scalar_remote(DDI_Patch *dAPatch, double alpha, int rank) {
  const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
  char ack = 39;

  DDI_Send_request(dAPatch, &rank, NULL);
  /* DDI_ARR_ZERO does not need to send alpha */
  if (dAPatch->oper != DDI_ARR_ZERO) Comm_send(&alpha, sizeof(double), rank, comm);
  Comm_recv(&ack, 1, rank, comm);
}

void DDI_ARR_scalar_server(DDI_Patch *dAPatch, int rank) {
  double alpha;
  char ack = 57;
  const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);

  /* DDI_ARR_ZERO does not need to receive alpha */
  if (dAPatch->oper == DDI_ARR_ZERO) alpha = (double)0;
  else Comm_recv(&alpha, sizeof(double), rank, comm);
  DDI_ARR_scalar_local(dAPatch, alpha);  
  Comm_send(&ack, 1, rank, comm);
}

