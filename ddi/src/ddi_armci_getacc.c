#include "ddi_base.h"

/** @see ddi_armci.h */
void DDI_ARMCI_GetAcc(DDI_Patch *patch, void *buf) {
  int handle = patch->handle;
  int pid = patch->jlo; /* patch has dimensions (0:0,pid:pid) */
  int armciPid;
  void *prem;
  long local, value = (long)(*(double*)buf);

  /* Assumption is made here that DDI_GetAcc is only called by DDI_PidDLB_next.
     If called with parameters differing from that of DDI_PidDLB_next, result will
     be wrong. */
  DDI_ARMCI_Acquire(handle, pid, DDI_WRITE_ACCESS, (void**)&prem, &armciPid);
  ARMCI_Rmw(ARMCI_FETCH_AND_ADD_LONG, (void*)(&local), prem, value, armciPid);
  DDI_ARMCI_Release(handle, pid, DDI_WRITE_ACCESS);
  *(double*)buf = (double)local;
}
