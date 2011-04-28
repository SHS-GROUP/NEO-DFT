#include "ddi_base.h"

/** @see ddi_armci.h */
void DDI_ARMCI_Memory_init(size_t size) {
  int code;
  const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_COMM_WORLD);
  
  code = ARMCI_Malloc(gv(armci_mem_addr), size);
  if (code != 0) {
    fprintf(DDI_STDERR, "%s: ARMCI_Malloc(%p, %z) returned %i\n",
	    DDI_Id(), gv(armci_mem_addr), size, code);
    DDI_Error(DDI_ARMCI_MEMORY_INIT_ERROR, DDI_ARMCI_MEMORY_INIT_ERROR_MESSAGE);
  }
  gv(dda_index) = (DDA_Index*)gv(armci_mem_addr)[comm->me];

  code = ARMCI_Create_mutexes(MAX_DD_ARRAYS);
  if (code != 0) {
    fprintf(DDI_STDERR, "%s: ARMCI_Create_mutexes(%d) returned %i\n",
	    DDI_Id(), MAX_DD_ARRAYS, code);
    DDI_Error(DDI_ARMCI_MEMORY_INIT_ERROR, DDI_ARMCI_MEMORY_INIT_ERROR_MESSAGE);
  }
}

/** @see ddi_armci.h */
void DDI_ARMCI_Counters_init() {
  int i;
  const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_COMM_WORLD);

  /* The offsets should be the same on all processes. */
  size_t dlb_counter_offset = (size_t)(gv(dlb_counter)) - (size_t)(gv(dda_index));
  size_t gdlb_counter_offset = (size_t)(gv(gdlb_counter)) - (size_t)(gv(dda_index));
  for (i = 0; i < comm->np; ++i) {
    gv(armci_dlb_counter)[i] = gv(armci_mem_addr)[i] + dlb_counter_offset;
    gv(armci_gdlb_counter)[i] = gv(armci_mem_addr)[i] + gdlb_counter_offset;
  }

  ARMCI_PutValueInt(0, gv(dlb_counter), comm->me);
  ARMCI_PutValueInt(0, gv(gdlb_counter), comm->me);
}

/** @see ddi_armci.h */
void DDI_ARMCI_DLBReset() {
    const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
    int pid = comm->me;
    int armciPid = comm->global_pid[pid];
    if (pid == 0) ARMCI_PutValueInt(0, gv(dlb_counter), armciPid);
}

/** @see ddi_armci.h */
void DDI_ARMCI_DLBNext(size_t *counter) {
    long buf;
    const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
    int armciPid = comm->global_pid[0];
    /* long is used in case signed int is too small */
    ARMCI_Rmw(ARMCI_FETCH_AND_ADD_LONG, (void*)(&buf), gv(armci_dlb_counter)[armciPid], 1, armciPid);
    *counter = (size_t)buf;
}

/** @see ddi_armci.h */
void DDI_ARMCI_GDLBReset() {
    const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_COMM_WORLD);
    int armciPid = comm->me;
    if (armciPid == 0) ARMCI_PutValueInt(0, gv(gdlb_counter), 0);
}

/** @see ddi_armci.h */
void DDI_ARMCI_GDLBNext(size_t *counter) {
    int tmp;
    const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
    if (comm->me == 0) ARMCI_Rmw(ARMCI_FETCH_AND_ADD, &tmp, gv(armci_gdlb_counter)[0], 1, 0);
    MPI_Bcast(&tmp, sizeof(int), MPI_BYTE, 0, comm->compute_comm);
    *counter = (size_t)tmp;
}

/** @see ddi_armci.h */
void DDI_ARMCI_Index_create(DDA_Index *index, int handle) {
  DDA_Remote_Index *remoteIndex = gv(dda_remote_index)[handle];
  const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);

  /* remoteIndex[pid] is indexed by DDI_WORKING_COMM process id. */
  remoteIndex[comm->me].offset = index[handle].offset;
  remoteIndex[comm->me].mutex = handle;
  /* remoteIndex[pid].offset may differ from process to process */
  MPI_Allgather(&remoteIndex[comm->me], sizeof(DDA_Remote_Index), MPI_BYTE,
		remoteIndex, sizeof(DDA_Remote_Index), MPI_BYTE, comm->compute_comm);
}

/** @see ddi_armci.h */
void DDI_ARMCI_Acquire(int handle, int pid, int ltype, void **array, int *armciPid) {
  DDA_Remote_Index *remoteIndex = gv(dda_remote_index)[handle];
  char *buf = NULL;
  int commId = gv(dda_index)[handle].commId;
  const DDI_Comm *comm = (const DDI_Comm*)Comm_find(commId);
  int armciPid_;
  int mutex;

  if (pid < 0) pid = comm->me;
  armciPid_ = comm->global_pid[pid];
  if (armciPid != NULL) *armciPid = armciPid_;
  mutex = remoteIndex[pid].mutex;

  *array = (void*)((char*)(gv(armci_mem_addr)[armciPid_]) + remoteIndex[pid].offset);

#if defined DDI_ARMCI_LOCK
  DDI_ARMCI_Lock(mutex, armciPid_);
#endif
}

/** @see ddi_armci.h */
void DDI_ARMCI_Release(int handle, int pid, int ltype) {
  DDA_Remote_Index *remoteIndex = gv(dda_remote_index)[handle];
  int commid = gv(dda_index)[handle].commId;
  const DDI_Comm *comm = (const DDI_Comm*)Comm_find(commid);
  int armciPid_;
  int mutex;

  if (pid < 0) pid = comm->me;
  armciPid_ = comm->global_pid[pid];
  mutex = remoteIndex[pid].mutex;
  
#if defined DDI_ARMCI_LOCK
  DDI_ARMCI_Unlock(mutex, armciPid_);
#endif
}

/** @see ddi_armci.h */
inline void DDI_ARMCI_Lock(int mutex, int pid) {
    ARMCI_Lock(mutex, pid);
}

/** @see ddi_armci.h */
inline void DDI_ARMCI_Unlock(int mutex, int pid) {
    ARMCI_Unlock(mutex, pid);
}

/** @see ddi_armci.h */
void DDI_ARMCI_Barrier(MPI_Comm comm) {
  ARMCI_WaitAll();
  MPI_Barrier(comm);
}

/** @see ddi_armci.h */
void DDI_ARMCI_Memory_finalize() {
  int code;
 
  code = ARMCI_Free(gv(dda_index));
  if (code != 0) {
    fprintf(DDI_STDERR, "%s: ARMCI_Free(%p) returned %i\n", DDI_Id(), gv(dda_index), code);
    DDI_Error(DDI_ARMCI_MEMORY_FINALIZE_ERROR, DDI_ARMCI_MEMORY_FINALIZE_ERROR_MESSAGE);
  }
  
  code = ARMCI_Destroy_mutexes();
  if (code != 0) {
    fprintf(DDI_STDERR, "%s: ARMCI_Destroy_mutexes() returned %i\n", DDI_Id(), code);
    DDI_Error(DDI_ARMCI_MEMORY_FINALIZE_ERROR, DDI_ARMCI_MEMORY_FINALIZE_ERROR_MESSAGE);
  }
}

/** @see ddi_armci.h */
void DDI_ARMCI_Finalize() {
  ARMCI_Finalize();
}
