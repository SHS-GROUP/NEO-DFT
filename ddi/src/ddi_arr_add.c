#include "ddi_base.h"

/** C = a*A + b*B */
void DDI_ARR_add_(DDI_Patch *dAPatch, double alpha, DDI_Patch *dBPatch, double beta, DDI_Patch *dCPatch) {
  int i;
  DDI_Patch dASubPatch[MAX_SMP_PROCS];
  DDI_Patch dBSubPatch[MAX_SMP_PROCS];
  DDI_Patch dCSubPatch[MAX_SMP_PROCS];
  int dASubPatchRank[MAX_SMP_PROCS], dANSubPatch;
  int dBSubPatchRank[MAX_SMP_PROCS], dBNSubPatch;
  int dCSubPatchRank[MAX_SMP_PROCS], dCNSubPatch;
  const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
  int rank = comm->me, node = comm->my;

  DDI_Subpatch_node(dAPatch->handle, dAPatch, &dANSubPatch, dASubPatchRank, dASubPatch, node);
  DDI_Subpatch_node(dBPatch->handle, dBPatch, &dBNSubPatch, dBSubPatchRank, dBSubPatch, node);
  DDI_Subpatch_node(dCPatch->handle, dCPatch, &dCNSubPatch, dCSubPatchRank, dCSubPatch, node);

  /* number of A, B, and C segments must be the same */
  if (dANSubPatch != dBNSubPatch || dBNSubPatch != dCNSubPatch) {
    Fatal_error(1000);
  }

  for (i = 0; i < dANSubPatch; ++i) {
    if (rank == dASubPatchRank[i]) {
      /* ensure segments are aligned on the processor */
      if (dBSubPatchRank[i] != dASubPatchRank[i] || dCSubPatchRank[i] != dASubPatchRank[i]) {
	Fatal_error(1000);
      }
      
#ifdef DDI_ARR_DMA
      DDI_ARR_add_local(&dASubPatch[i], alpha, &dBSubPatch[i], beta, &dCSubPatch[i]);
#else
      DDI_ARR_add_remote(&dASubPatch[i], alpha, &dBSubPatch[i], beta, &dCSubPatch[i], rank);
#endif
    }
  }
  Comm_sync(1000, comm);
}

void DDI_ARR_add_local(DDI_Patch *dAPatch, double alpha, DDI_Patch *dBPatch, double beta, DDI_Patch *dCPatch) {
  DDA_Index *Index = gv(dda_index);
  int dA = dAPatch->handle, dB = dBPatch->handle, dC = dCPatch->handle;
  double *dALocal, *dBLocal, *dCLocal;
  int op = dAPatch->oper;

  int dALdaLocal, dAiLocal, dAjLocal, dAmLocal, dAnLocal;
  int dBLdaLocal, dBiLocal, dBjLocal, dBmLocal, dBnLocal;
  int dCLdaLocal, dCiLocal, dCjLocal, dCmLocal, dCnLocal;

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

  /* C segment dimensions */
  dCLdaLocal = Index[dC].ihi - Index[dC].ilo + 1;
  dCiLocal = dCPatch->ilo - Index[dC].ilo;
  dCjLocal = dCPatch->jlo - Index[dC].jlo;
  dCmLocal = dCPatch->ihi - dCPatch->ilo + 1;
  dCnLocal = dCPatch->jhi - dCPatch->jlo + 1;

# if defined USE_SYSV
    if(USING_DATA_SERVERS()) {
      DDI_Fence_check(dA);
      DDI_Fence_check(dB);
      DDI_Fence_check(dC);
    }
# endif

    DDI_Acquire(Index, dC, DDI_WRITE_ACCESS, (void**)&dCLocal);
    if (dA == dC) dALocal = dCLocal;
    else DDI_Acquire(Index, dA, DDI_READ_ACCESS, (void**)&dALocal);
    if (dB == dC) dBLocal = dCLocal;
    else DDI_Acquire(Index, dB, DDI_READ_ACCESS, (void**)&dBLocal);

    /* Matrix segments have to be the same. */
    if ((dAmLocal != dBmLocal || dBmLocal != dCmLocal) ||
	(dAnLocal != dBnLocal || dBnLocal != dCnLocal)) {
      Fatal_error(1000);
    }
    mmadd(alpha, dALocal, dALdaLocal, dAiLocal, dAjLocal,
	  beta, dBLocal, dBLdaLocal, dBiLocal, dBjLocal,
	  dCLocal, dCLdaLocal, dCiLocal, dCjLocal, dAmLocal, dAnLocal);
    
    DDI_Release(Index, dC, DDI_WRITE_ACCESS);
    if (dA != dC) DDI_Release(Index, dA, DDI_READ_ACCESS);
    if (dB != dC) DDI_Release(Index, dB, DDI_READ_ACCESS);
}

void DDI_ARR_add_remote(DDI_Patch *dAPatch, double alpha, DDI_Patch *dBPatch, double beta, DDI_Patch *dCPatch, int rank) {
  DDI_Patch request;
  DDI_ARR_Patch sendbuf[3];
  const DDI_Comm *comm = (const DDI_Comm *)Comm_find(DDI_WORKING_COMM);
  char ack = 39;

  /* send request  */
  request.oper = DDI_ARR_ADD;
  DDI_Send_request(&request, &rank, NULL);
  /* send A, B, and C segment, and wait for completion  */
  sendbuf[0].alpha = alpha;
  sendbuf[1].alpha = beta;
  sendbuf[0].patch = *dAPatch;
  sendbuf[1].patch = *dBPatch;
  sendbuf[2].patch = *dCPatch;
  Comm_send(sendbuf, sizeof(DDI_ARR_Patch)*3, rank, comm);
  Comm_recv(&ack, 1, rank, comm);
}

void DDI_ARR_add_server(DDI_Patch *request, int rank) {
  DDI_ARR_Patch recvbuf[3];
  DDI_Patch *dAPatch, *dBPatch, *dCPatch;
  double alpha, beta;
  const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
  char ack = 39;

  /* receive A, B, and C segments */
  Comm_recv(recvbuf, sizeof(DDI_ARR_Patch)*3, rank, comm);
  alpha = recvbuf[0].alpha;
  beta = recvbuf[1].alpha;
  dAPatch = &(recvbuf[0].patch);
  dBPatch = &(recvbuf[1].patch);
  dCPatch = &(recvbuf[2].patch);
  DDI_ARR_add_local(dAPatch, alpha, dBPatch, beta, dCPatch);
  Comm_send(&ack, 1, rank, comm);
}

