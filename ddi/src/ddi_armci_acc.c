#include "ddi_base.h"

/** @see ddi_armci.h */
void DDI_ARMCI_Acc(DDI_Patch *patch, double alpha, void *buf) {
    int handle = patch->handle;
    int commId = gv(dda_index)[handle].commId;
    int i, nops = 0;
    int nsubp, ranks[MAX_NODES];
    DDI_Patch subp[MAX_NODES];
    char *working_buffer = (char*)buf;
    
    DDI_Subpatch_comm(handle, patch, &nsubp, ranks, subp, commId);
    
    for(i=0; i<nsubp; i++) {
#if FULL_SMP
	nops += DDI_ARMCI_Acc_domain_SMP(&subp[i], alpha, working_buffer, ranks[i]);
#else
	nops += DDI_ARMCI_Acc_proc(&subp[i], alpha, working_buffer, ranks[i]);
#endif
	working_buffer += subp[i].size;
    }
    
#if defined DDI_ARMCI_IMPLICIT_NBACC && defined DDI_ARMCI_IMPLICIT_WAIT
    /* wait for implicit non-blocking operations */
    ARMCI_WaitAll();
#endif

    return;
}

/** @see ddi_armci.h */
inline int DDI_ARMCI_Acc_domain_SMP(DDI_Patch *patch, double alpha, void *buf, int domain) {
    int handle = patch->handle;
    int i, nsubp, nops = 0;
    int ranks[MAX_SMP_PROCS];
    DDI_Patch subp[MAX_SMP_PROCS];
    char *src = (char*)buf;
    
    DDI_Subpatch_node(handle, patch, &nsubp, ranks, subp, domain);
    
    for (i = 0; i < nsubp; ++i) {
	nops += DDI_ARMCI_Acc_proc(&subp[i], alpha, src, ranks[i]);
	src += subp[i].size;
    }
    return nops;
}

/** @see ddi_armci.h */
inline int DDI_ARMCI_Acc_proc(DDI_Patch *patch, double alpha, void *buf, int pid) {
    int handle = patch->handle;
    int nops = 1;

    int trows,tcols,nrows,ncols;
    size_t offset;
    char *dst,*src = (char*)buf;
    int src_stride_arr[2],dst_stride_arr[2],count[2];
    int stride_levels = 1;
    int armciPid;
    
    trows = gv(dda_index)[handle].nrows;
    tcols = gv(pcmap)[handle][pid+1] - gv(pcmap)[handle][pid];
    nrows = patch->ihi - patch->ilo + 1;
    ncols = patch->jhi - patch->jlo + 1;
    
    offset = (patch->jlo - gv(pcmap)[handle][pid])*trows + patch->ilo;
    offset *= sizeof(double);
    
    DDI_ARMCI_Acquire(handle, pid, DDI_WRITE_ACCESS, (void**)&dst, &armciPid);
    dst += offset;
    
    if (nrows == trows) {
	src_stride_arr[0] = sizeof(double)*nrows;
	dst_stride_arr[0] = sizeof(double)*trows;
	count[0] = patch->size;
	stride_levels = 0;
	
#if defined DDI_ARMCI_IMPLICIT_NBACC
	/* Apparantely ARMCI_Acc is not always present */
	/*ARMCI_Acc(ARMCI_ACC_DBL, alpha, (void*)src, (void*)dst, subp[i].size, armciPid, NULL);*/
	ARMCI_NbAccS(ARMCI_ACC_DBL, &alpha,
		     (void*)src, src_stride_arr,
		     (void*)dst, dst_stride_arr,
		     count, stride_levels, armciPid, NULL);
#else
	/* Apparantely ARMCI_Acc is not always present */
	/*ARMCI_Acc(ARMCI_ACC_DBL, 1.0, (void*)src, (void*)dst, subp[i].size, armciPid);*/
	ARMCI_AccS(ARMCI_ACC_DBL, &alpha,
		   (void*)src, src_stride_arr,
		   (void*)dst, dst_stride_arr,
		   count, stride_levels, armciPid);
#endif
    }
    else {
        /* i dimensions */
	src_stride_arr[0] = sizeof(double)*nrows;
	dst_stride_arr[0] = sizeof(double)*trows;
	/* j dimensions */
	src_stride_arr[1] = src_stride_arr[0]*ncols;
	dst_stride_arr[1] = dst_stride_arr[0]*tcols;
	/* block size count, first dimension must be in bytes */
	count[0] = sizeof(double)*nrows;
	count[1] = ncols;
	stride_levels = 1;
	
#if defined DDI_ARMCI_IMPLICIT_NBACC
	ARMCI_NbAccS(ARMCI_ACC_DBL, &alpha,
		     (void*)src, src_stride_arr,
		     (void*)dst, dst_stride_arr,
		     count, stride_levels, armciPid, NULL);
#else
	ARMCI_AccS(ARMCI_ACC_DBL, &alpha,
		   (void*)src, src_stride_arr,
		   (void*)dst, dst_stride_arr,
		   count, stride_levels, armciPid);
#endif
    }
    
    DDI_ARMCI_Release(handle, pid, DDI_WRITE_ACCESS);

    return nops;
}
