#include "ddi_base.h"

void DDI_ARR_select_(DDI_Patch *dAPatch, double *alpha, int *index) {
  int i;
  int op = dAPatch->oper;
  DDI_Patch dASubPatch[MAX_PROCESSORS];
  int dASubPatchRank[MAX_PROCESSORS], dANSubPatch;
  DDI_ARR_Element recvbuf, *element = NULL;
  const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
  int rank = comm->me;
  int node = comm->my;
  int smpRank = comm->me_local;
  int smpRoot;

  DDI_Subpatch_node(dAPatch->handle, dAPatch, &dANSubPatch, dASubPatchRank, dASubPatch, node);

  for (i = 0; i < dANSubPatch; ++i) {
    if (dASubPatchRank[i] == rank) {
      if (element == NULL) element = (DDI_ARR_Element*)malloc(sizeof(DDI_ARR_Element));
#ifdef DDI_ARR_DMA
      DDI_ARR_select_local(&dASubPatch[i], element);
#else
      DDI_ARR_select_remote(&dASubPatch[i], element, rank);
#endif
    }
  }

  /* Now selected elements need to be gathered from other nodes.
     Point-to-point is used because patch may not be located on all processors,
     but perhaps consider change to collective operation.
     IRecv or Recvany are not used because they seem to give problems.
  */

  /* select intra-node element */
  smpRoot = comm->node_master[node];
  if (smpRank == 0) {
    for (i = 0; i < dANSubPatch; ++i) {
      if (dASubPatchRank[i] == rank) continue;
      Comm_recv(&recvbuf, sizeof(DDI_ARR_Element), dASubPatchRank[i], comm);
      DDI_ARR_Element_select(element, &recvbuf, op);
    }
  }
  else if (element != NULL) {
    Comm_send(element, sizeof(DDI_ARR_Element), smpRoot, comm);
  }

  /* select inter-node element */
  DDI_Subpatch(dAPatch->handle, dAPatch, &dANSubPatch, dASubPatchRank, dASubPatch);
  if (rank == 0) {
    for (i = 0; i < dANSubPatch; ++i) {
      if (dASubPatchRank[i] == rank) continue;
      smpRoot = comm->node_master[dASubPatchRank[i]];
      Comm_recv(&recvbuf, sizeof(DDI_ARR_Element), smpRoot, comm);
      DDI_ARR_Element_select(element, &recvbuf, op);
    }
  }
  else if (element != NULL && rank == smpRoot) {
    Comm_send(element, sizeof(DDI_ARR_Element), 0, comm);
  }

  /* broadcast selected element */
  if (element == NULL) element = (DDI_ARR_Element*)malloc(sizeof(DDI_ARR_Element));
  Comm_bcast(element, sizeof(DDI_ARR_Element), 0, comm);

  *alpha = element->alpha;
  index[0] = element->index[0];
  index[1] = element->index[1];
  free(element);
}


void DDI_ARR_select_local(DDI_Patch *dAPatch, DDI_ARR_Element *element) {
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
  case DDI_ARR_MIN:
    mmin(dALocal, dALdaLocal, dAiLocal, dAjLocal, dAmLocal, dAnLocal, &(element->alpha), element->index);
    break;
    
  case DDI_ARR_MAX:
    mmax(dALocal, dALdaLocal, dAiLocal, dAjLocal, dAmLocal, dAnLocal, &(element->alpha), element->index);
    break;
  } /* switch */
  
  /* adjust element indices to global values */
  element->index[0] += Index[dA].ilo;
  element->index[1] += Index[dA].jlo;
  
  DDI_Release(Index, dA, DDI_READ_ACCESS);
}

void DDI_ARR_select_remote(DDI_Patch *dAPatch, DDI_ARR_Element *element, int rank) {
  const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);

  DDI_Send_request(dAPatch, &rank, NULL);
  Comm_recv(element, sizeof(DDI_ARR_Element), rank, comm);
}

void DDI_ARR_select_server(DDI_Patch *dAPatch, int rank) {
  DDI_ARR_Element element;
  const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);

  DDI_ARR_select_local(dAPatch, &element);
  Comm_send(&element, sizeof(DDI_ARR_Element), rank, comm);
}

void DDI_ARR_Element_select(DDI_ARR_Element *a, DDI_ARR_Element *b, int op) {
  if (a == NULL) {
    a = (DDI_ARR_Element*)malloc(sizeof(DDI_ARR_Element));
    *a = *b;
  }
  else {
    switch(op) {
    case DDI_ARR_MIN:
      if (b->alpha < a->alpha) *a = *b;
      break;
    case DDI_ARR_MAX:
      if (b->alpha > a->alpha) *a = *b;
      break;
    default:
      Fatal_error(1000);
      break;
    }
  }
}
