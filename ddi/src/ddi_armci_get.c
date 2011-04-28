#include "ddi_base.h"

/** @see ddi_armci.h */
void DDI_ARMCI_Get(DDI_Patch *patch, void *buf) {
    int handle = patch->handle;
    int commid = gv(dda_index)[handle].commId;
    int i, nops = 0;
    int nsubp, ranks[MAX_NODES];
    DDI_Patch subp[MAX_NODES];
    char *working_buffer = (char*)buf;

    DDI_Subpatch_comm(handle, patch, &nsubp, ranks, subp, commid);

    for(i=0; i<nsubp; i++) {
#if FULL_SMP
        nops += DDI_ARMCI_Get_domain_SMP(&subp[i], working_buffer, ranks[i]);
#else
        nops += DDI_ARMCI_Get_proc(&subp[i], working_buffer, ranks[i]);
#endif
        working_buffer += subp[i].size;
    }
    
#if defined DDI_ARMCI_IMPLICIT_NBGET && defined DDI_ARMCI_IMPLICIT_WAIT
    /* wait for implicit non-blocking operations */
    ARMCI_WaitAll();
#endif

    return;
}

/** @see ddi_armci.h */
inline int DDI_ARMCI_Get_domain_SMP(DDI_Patch *patch, void *buf, int domain) {
    int handle = patch->handle;
    int i, nsubp, nops = 0;
    int ranks[MAX_SMP_PROCS];
    DDI_Patch subp[MAX_SMP_PROCS];
    char *dst = (char*)buf;
    
    DDI_Subpatch_node(handle, patch, &nsubp, ranks, subp, domain);
    
    for (i = 0; i < nsubp; ++i) {
	nops += DDI_ARMCI_Get_proc(&subp[i], dst, ranks[i]);
	dst += subp[i].size;
    }
    
    return nops;
}

/** @see ddi_armci.h */
inline int DDI_ARMCI_Get_proc(DDI_Patch *patch, void *buf, int pid) {
    int handle = patch->handle;
    int nops = 1;

    int trows,tcols,nrows,ncols;
    size_t offset;
    char *src,*dst = (char*)buf;
    int src_stride_arr[2],dst_stride_arr[2],count[2];
    int stride_levels = 1;
    int armciPid;
    
    trows = gv(dda_index)[handle].nrows;
    tcols = gv(pcmap)[handle][pid+1] - gv(pcmap)[handle][pid];
    nrows = patch->ihi - patch->ilo + 1;
    ncols = patch->jhi - patch->jlo + 1;
    
    offset = (patch->jlo - gv(pcmap)[handle][pid])*trows + patch->ilo;
    offset *= sizeof(double);
    
    DDI_ARMCI_Acquire(handle, pid, DDI_WRITE_ACCESS, (void**)&src, &armciPid);
    src += offset;
    
    if (nrows == trows) {
#if defined DDI_ARMCI_IMPLICIT_NBGET
	ARMCI_NbGet((void*)src, (void*)dst, patch->size, armciPid, NULL);
#else
	ARMCI_Get((void*)src, (void*)dst, patch->size, armciPid);
#endif
    }
    else {
        /* i dimensions */
	src_stride_arr[0] = sizeof(double)*trows;
	dst_stride_arr[0] = sizeof(double)*nrows;
	/* j dimensions */
	src_stride_arr[1] = src_stride_arr[0]*tcols;
	dst_stride_arr[1] = dst_stride_arr[0]*ncols;
	/* block size count, first dimension must be in bytes */
	count[0] = sizeof(double)*nrows;
	count[1] = ncols;
	stride_levels = 1;
	
#if defined DDI_ARMCI_IMPLICIT_NBGET
	ARMCI_NbGetS((void*)src, src_stride_arr, (void*)dst, dst_stride_arr, count, stride_levels, armciPid, NULL);
#else
	ARMCI_GetS((void*)src, src_stride_arr, (void*)dst, dst_stride_arr, count, stride_levels, armciPid);
#endif
    }
    
    DDI_ARMCI_Release(handle, pid, DDI_WRITE_ACCESS);

    return nops;
}
