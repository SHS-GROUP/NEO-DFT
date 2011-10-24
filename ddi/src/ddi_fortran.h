/* ------------------------------------------------------------- *\
   External DDI Functions called by FORTRAN
   ========================================

   All DDI subroutines are written in C.  DDI subroutines called
   from FORTRAN go through a C wrapper.  All DDI subroutines are
   prefixed with DDI_ and all FORTRAN wrapper subroutines are
   prefixed with F77_.
   
   The following C definitions substitute the FORTRAN wrapper
   subroutine name with the correct machine dependent FORTRAN
   external name.
   
   Note: FORTRAN externals are generally all lowercase, but may
   be uppercase.  Define F77_UPPERCASE if you machine's FORTRAN
   using uppercase externals.
   
   Author: Ryan M. Olson
   Date: December 22, 2002
   CVS $Id: ddi_fortran.h,v 1.2 2007/05/26 20:48:22 andrey Exp $
\* ------------------------------------------------------------- */
 # include "f77_extern.h"

 # if defined INT_SIZE
   typedef INT_SIZE int_f77;
 # else
 # error "INT_SIZE must be defined."
 # endif

 #  if defined F77_UPPERCASE
 #     define F77_PBeg                    F77_Extern(DDI_PBEG)
 #     define F77_Init                    F77_Extern(DDI_INIT)
 #     define F77_PEnd                    F77_Extern(DDI_PEND)
 #     define F77_Finalize                F77_Extern(DDI_FINALIZE)
 #     define F77_NProc                   F77_Extern(DDI_NPROC)
 #     define F77_NNode                   F77_Extern(DDI_NNODE)
 #     define F77_Memory                  F77_Extern(DDI_MEMORY)
 #     define F77_Create                  F77_Extern(DDI_CREATE)
 #     define F77_Create_custom           F77_Extern(DDI_CREATE_CUSTOM)
 #     define F77_Destroy                 F77_Extern(DDI_DESTROY)
 #     define F77_Zero                    F77_Extern(DDI_ZERO)
 #     define F77_Distrib                 F77_Extern(DDI_DISTRIB)
 #     define F77_NDistrib                F77_Extern(DDI_NDISTRIB)
 #     define F77_DLBReset                F77_Extern(DDI_DLBRESET)
 #     define F77_DLBNext                 F77_Extern(DDI_DLBNEXT)
 #     define F77_Send                    F77_Extern(DDI_SEND)
 #     define F77_Recv                    F77_Extern(DDI_RECV)
 #     define F77_Get                     F77_Extern(DDI_GET)
 #     define F77_Put                     F77_Extern(DDI_PUT)
 #     define F77_Acc                     F77_Extern(DDI_ACC)
 #     define F77_ISend                   F77_Extern(DDI_ISEND)
 #     define F77_IRecv                   F77_Extern(DDI_IRECV)
 #     define F77_IGet                    F77_Extern(DDI_IGET)
 #     define F77_IPut                    F77_Extern(DDI_IPUT)
 #     define F77_IAcc                    F77_Extern(DDI_IACC)
 #     define F77_Wait                    F77_Extern(DDI_WAIT)
 #     define F77_Sync                    F77_Extern(DDI_SYNC)
 #     define F77_BCast                   F77_Extern(DDI_BCAST)
 #     define F77_GSumF                   F77_Extern(DDI_GSUMF)
 #     define F77_GSumI                   F77_Extern(DDI_GSUMI)
 #     define F77_Rcvany                  F77_Extern(DDI_RCVANY)
 #     define F77_Output                  F77_Extern(DDI_OUTPUT)
 #     define F77_Level                   F77_Extern(DDI_LEVEL)
 #     define F77_Debug                   F77_Extern(DDI_DEBUG)
 #     define F77_Timer_reset             F77_Extern(DDI_TIMER_RESET)
 #     define F77_Timer_output            F77_Extern(DDI_TIMER_OUTPUT)
 #     define F77_SMP_NProc               F77_Extern(DDI_SMP_NPROC)
 #     define F77_Group_create            F77_Extern(DDI_GROUP_CREATE)
 #     define F77_Group_create_custom     F77_Extern(DDI_GROUP_CREATE_CUSTOM)
 #     define F77_Scope                   F77_Extern(DDI_SCOPE)
 #     define F77_AScope                  F77_Extern(DDI_ASCOPE)
 #     define F77_NGroup                  F77_Extern(DDI_NGROUP)
 #     define F77_GDLBReset               F77_Extern(DDI_GDLBRESET)
 #     define F77_GDLBNext                F77_Extern(DDI_GDLBNEXT)
 #     define F77_ProcDLB_create          F77_Extern(DDI_PROCDLB_CREATE)
 #     define F77_ProcDLB_next            F77_Extern(DDI_PROCDLB_NEXT)
 #     define F77_ProcDLB_reset           F77_Extern(DDI_PROCDLB_RESET)
 #     define F77_ProcDLB_destroy         F77_Extern(DDI_PROCDLB_DESTROY)
/* DDI_ARR subroutines */
 #     define F77_ARR_zero                F77_Extern(DDI_ARR_ZERO)
 #     define F77_ARR_fill                F77_Extern(DDI_ARR_FILL)
 #     define F77_ARR_scale               F77_Extern(DDI_ARR_SCALE)
 #     define F77_ARR_min                 F77_Extern(DDI_ARR_MIN)
 #     define F77_ARR_max                 F77_Extern(DDI_ARR_MAX)
 #     define F77_ARR_dot                 F77_Extern(DDI_ARR_DOT)
 #     define F77_ARR_add                 F77_Extern(DDI_ARR_ADD)
 #     define F77_ARR_acc                 F77_Extern(DDI_ARR_ACC)
 # if defined WINTEL
 #     define F77_SMP_Create              F77_Extern(DDI_SMP_CREATE)
 #     define F77_SMP_Destroy             F77_Extern(DDI_SMP_DESTROY)
 #     define F77_SMP_Offset              F77_Extern(DDI_SMP_OFFSET)
 #     define F77_SMP_sync                F77_Extern(DDI_SMP_SYNC)
 #     define F77_BCast_smp               F77_Extern(DDI_SMP_BCAST)
 #     define F77_GSumF_smp               F77_Extern(DDI_SMP_GSUMF)
 #     define F77_GSumF_node              F77_Extern(DDI_MASTERS_GSUMF)
 #     define F77_BCast_node              F77_Extern(DDI_MASTERS_BCAST)
 #     define F77_Get_comm                F77_Extern(DDI_GET_COMM)
 #     define F77_Put_comm                F77_Extern(DDI_PUT_COMM)
 # endif
 #  else
 #     define F77_PBeg                    F77_Extern(ddi_pbeg)
 #     define F77_Init                    F77_Extern(ddi_init)
 #     define F77_PEnd                    F77_Extern(ddi_pend)
 #     define F77_Finalize                F77_Extern(ddi_finalize)
 #     define F77_NProc                   F77_Extern(ddi_nproc)
 #     define F77_NNode                   F77_Extern(ddi_nnode)
 #     define F77_Memory                  F77_Extern(ddi_memory)
 #     define F77_Create                  F77_Extern(ddi_create)
 #     define F77_Create_custom           F77_Extern(ddi_create_custom)
 #     define F77_Destroy                 F77_Extern(ddi_destroy)
 #     define F77_Zero                    F77_Extern(ddi_zero)
 #     define F77_Distrib                 F77_Extern(ddi_distrib)
 #     define F77_NDistrib                F77_Extern(ddi_ndistrib)
 #     define F77_DLBReset                F77_Extern(ddi_dlbreset)
 #     define F77_DLBNext                 F77_Extern(ddi_dlbnext)
 #     define F77_Send                    F77_Extern(ddi_send)
 #     define F77_Recv                    F77_Extern(ddi_recv)
 #     define F77_Get                     F77_Extern(ddi_get)
 #     define F77_Put                     F77_Extern(ddi_put)
 #     define F77_Acc                     F77_Extern(ddi_acc)
 #     define F77_ISend                   F77_Extern(ddi_isend)
 #     define F77_IRecv                   F77_Extern(ddi_irecv)
 #     define F77_IGet                    F77_Extern(ddi_iget)
 #     define F77_IPut                    F77_Extern(ddi_iput)
 #     define F77_IAcc                    F77_Extern(ddi_iacc)
 #     define F77_Wait                    F77_Extern(ddi_wait)
 #     define F77_Sync                    F77_Extern(ddi_sync)
 #     define F77_BCast                   F77_Extern(ddi_bcast)
 #     define F77_GSumF                   F77_Extern(ddi_gsumf)
 #     define F77_GSumI                   F77_Extern(ddi_gsumi)
 #     define F77_Rcvany                  F77_Extern(ddi_rcvany)
 #     define F77_Output                  F77_Extern(ddi_output)
 #     define F77_Level                   F77_Extern(ddi_level)
 #     define F77_Debug                   F77_Extern(ddi_debug)
 #     define F77_Timer_reset             F77_Extern(ddi_timer_reset)
 #     define F77_Timer_output            F77_Extern(ddi_timer_output)
 #     define F77_SMP_NProc               F77_Extern(ddi_smp_nproc)
 #     define F77_Group_create            F77_Extern(ddi_group_create)
 #     define F77_Group_create_custom     F77_Extern(ddi_group_create_custom)
 #     define F77_Scope                   F77_Extern(ddi_scope)
 #     define F77_AScope                  F77_Extern(ddi_ascope)
 #     define F77_NGroup                  F77_Extern(ddi_ngroup)
 #     define F77_GDLBReset               F77_Extern(ddi_gdlbreset)
 #     define F77_GDLBNext                F77_Extern(ddi_gdlbnext)
 #     define F77_ProcDLB_create          F77_Extern(ddi_procdlb_create)
 #     define F77_ProcDLB_next            F77_Extern(ddi_procdlb_next)
 #     define F77_ProcDLB_reset           F77_Extern(ddi_procdlb_reset)
 #     define F77_ProcDLB_destroy         F77_Extern(ddi_procdlb_destroy)
 #     define F77_SMP_sync                F77_Extern(ddi_smp_sync)
 #     define F77_BCast_smp               F77_Extern(ddi_smp_bcast)
 #     define F77_GSumI_smp               F77_Extern(ddi_smp_gsumi)
 #     define F77_GSumF_smp               F77_Extern(ddi_smp_gsumf)
 #     define F77_BCast_node              F77_Extern(ddi_masters_bcast)
 #     define F77_GSumI_node              F77_Extern(ddi_masters_gsumi)
 #     define F77_GSumF_node              F77_Extern(ddi_masters_gsumf)
 #     define F77_SMP_Create              F77_Extern(ddi_smp_create)
 #     define F77_SMP_Destroy             F77_Extern(ddi_smp_destroy)
 #     define F77_SMP_Offset              F77_Extern(ddi_smp_offset)
 #     define F77_Addr_test               F77_Extern(ddi_addr_test)
 #     define F77_DB_Create               F77_Extern(ddi_db_create)
 #     define F77_DB_Read                 F77_Extern(ddi_db_read)
 #     define F77_DB_Write                F77_Extern(ddi_db_write)
 #     define F77_Get_comm                F77_Extern(ddi_get_comm)
 #     define F77_Put_comm                F77_Extern(ddi_put_comm)
 #     define F77_Acc_comm                F77_Extern(ddi_acc_comm)
/* DDI_ARR subroutines */
 #     define F77_ARR_zero                F77_Extern(ddi_arr_zero)
 #     define F77_ARR_fill                F77_Extern(ddi_arr_fill)
 #     define F77_ARR_scale               F77_Extern(ddi_arr_scale)
 #     define F77_ARR_min                 F77_Extern(ddi_arr_min)
 #     define F77_ARR_max                 F77_Extern(ddi_arr_max)
 #     define F77_ARR_dot                 F77_Extern(ddi_arr_dot)
 #     define F77_ARR_add                 F77_Extern(ddi_arr_add)
 #     define F77_ARR_acc                 F77_Extern(ddi_arr_acc)
 #  endif
