/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Prototypes and Data Structures
 * 
 * Author: Ryan M. Olson
 * 10 Jun 2009 - RMO - give FENCE a value
 * 13 May 2010 - SS  - define a 64 bit integer type, for Windows' benefit
 * xx Jun 2010 - RMO - DDI_Patch (control message) extended to include alpha 
\* -------------------------------------------------------------------- */
 # include <stdio.h>
 # include <stdlib.h>
 # include "debug.h"

#define DDI_VERSION 3
#define DDI_SUBVERSION 0

/* ----------------- *\
   DDI Communicators
\* ----------------- */
 # define DDI_COMM_WORLD 0
 # define DDI_WORKING_COMM gv(ddi_working_comm)

/* -------------- *\
   DDI Data Types
\* -------------- */
# if defined WINDOWS64
 # define DDI_INT64  __int64
# else
 # define DDI_INT64  long
#endif

 # define DDI_DOUBLE 0
 # define DDI_INT    sizeof(int)
 # define DDI_LONG   sizeof(DDI_INT64)

/* ------------------- *\
   DDI Data Structures
\* ------------------- */
   typedef struct {
     int oper;
     int handle;
     int ilo;
     int ihi;
     int jlo;
     int jhi;
     size_t size;
     /*  possible scale factor, e.g. for DAXPY-type accumulates  */
     double alpha;

   # if defined DDI_LAPI
     int    cp_lapi_id;
     void  *cp_lapi_cntr;
     void  *cp_buffer_addr;
     size_t ds_buffer_size;
   # endif

   } DDI_Patch;

   typedef struct {
     int oper;
     int handle;
     int ilo;
     int ihi;
     int jlo;
     int jhi;

       //       size_t start;
       //       size_t end;
       size_t size; //message size in bytes
       size_t nelem; //total number of elements to be accumulated
       size_t nlocal; //number of local elements for subp[i]
     /*  possible scale factor, e.g. for DAXPY-type accumulates  */
     double alpha;

   } DDI_Scattered;


/* ----------------------- *\
   DDI Function Prototypes
\* ----------------------- */
   void DDI_Init(int,char**);
   void DDI_Memory(size_t);
   void DDI_Finalize();
   void DDI_Sync(int);
   void DDI_BCast(void*,size_t,int);
   void DDI_GSum(void*,size_t,int);
   void DDI_Create(int,int,int*);
   void DDI_Create_custom(int,int,int*,int*);
   void DDI_Destroy(int);
   void DDI_Zero(int);
   void DDI_Output(int);
   void DDI_Get(int,int,int,int,int,void*);
   void DDI_Put(int,int,int,int,int,void*);
   void DDI_Acc(int,int,int,int,int,void*);
   void DDI_GetAcc(int,int,int,int,int,void*);
   void DDI_GetP(int,DDI_Patch*,void*);
   void DDI_PutP(int,DDI_Patch*,void*);
   void DDI_AccP(int,DDI_Patch*,double,void*);   
   void DDI_Scatter_AccS(int,DDI_Scattered*,double,long*,void*);   
   void DDI_GetAccP(int,DDI_Patch*,void*);   
   void DDI_DLBReset();
   void DDI_DLBNext(size_t*);
   void DDI_Send(void*,size_t,int);
   void DDI_Recv(void*,size_t,int);
   void DDI_Recvany(void*,size_t,int*);

   void DDI_GetP_comm(int,DDI_Patch*,void*,int);
   void DDI_PutP_comm(int,DDI_Patch*,void*,int);

/* ----------------- *\
   New SMP Functions
\* ----------------- */
   void DDI_SMP_NProc(int*,int*);
   void DDI_SMP_Create(size_t,int*);
   void DDI_SMP_Destroy(int);
   void *DDI_SMP_Array(int);

/* ---------------------------- *\
   New non-blocking subroutines
\* ---------------------------- */
   void DDI_ISend(void *,size_t,int,int*);
   void DDI_IRecv(void *,size_t,int,int*);
   void DDI_Wait(int);

/* ---------------------- *\
   New Group DDI routines
\* ---------------------- */
   void DDI_Scope(int);
   void DDI_AScope(int);
   void DDI_Group_create(int,int*,int*,int*);
   void DDI_Group_create_custom(int,int*,int*,int*,int*);
   void DDI_GDLBReset();
   void DDI_GDLBNext(size_t*);

/* -------------------------------------- *\
   New dynamic load-balancing subroutines
\* -------------------------------------- */
   void DDI_ProcDLB_create(int*);
   void DDI_ProcDLB_destroy(int);
   void DDI_ProcDLB_reset(int);
   void DDI_ProcDLB_next(int,int,int*);

/* ---------------------- *\
   New timing subroutines
\* ---------------------- */
   void DDI_Timer_reset();
   void DDI_Timer_output();

/* -------------------------- *\
   Possible Inlined Functions
\* -------------------------- */
   void DDI_NProc(int*,int*);
   void DDI_NNode(int*,int*);
   void DDI_NGroup(int*,int*);
   void DDI_Distrib(int,int,int*,int*,int*,int*);
   void DDI_DistribP(int,int,DDI_Patch*);
   void DDI_NDistrib(int,int,int*,int*,int*,int*);
   void DDI_NDistribP(int,int,DDI_Patch*);

   void DDI_Scattered_DistribP(int,int,DDI_Scattered*);
   void DDI_Scattered_NDistribP(int,int,DDI_Scattered*);

/* -------------------------- *\
   DDI Communicator Routines
\* -------------------------- */
   void DDI_Comm_divide(int,int,int*);
   void DDI_Comm_create(int,int*,int,int,int,int*);
   void DDI_Comm_destroy(int);

   void DDI_NProc_comm(int,int*,int*);
   void DDI_NNode_comm(int,int*,int*);
   void DDI_Sync_comm(int,int);
   void DDI_GSum_comm(void*,size_t,int,int);
   void DDI_Send_comm(void*,size_t,int,int);
   void DDI_Recv_comm(void*,size_t,int,int);
   void DDI_BCast_comm(void*,size_t,int,int);

   void DDI_Sync_smp(int);
   void DDI_GSum_smp(void*,size_t,int);
   void DDI_BCast_smp(void*,size_t,int);

   void DDI_Sync_node(int);
   void DDI_GSum_node(void*,size_t,int);
   void DDI_BCast_node(void*,size_t,int);

/* -------------- *\
   DDI Parameters
\* -------------- */
 # define DDI_WORLD   0
 # define DDI_GROUP   1
 # define DDI_MASTERS 2


/* -------------- *\
   DDI Operations
\* -------------- */

/* Certain variables have an _OP suffix to avoid
 * namespace collision with fortran method names
 * when using ifort on Windows */

# if defined WINTEL 
 
 # define DDI_MEMORY_OP      0
 # define DDI_CREATE_OP      1
 # define DDI_DESTROY_OP     2
 # define DDI_GET_OP         3
 # define DDI_PUT_OP         4
 # define DDI_ACC_OP         5

# else

 # define DDI_MEMORY         0
 # define DDI_CREATE         1
 # define DDI_DESTROY        2
 # define DDI_GET            3
 # define DDI_PUT            4
 # define DDI_ACC            5

# endif 
 
 # define DDI_GETACC         6
 
# if defined WINTEL 
 
 # define DDI_DLBRESET_OP    7
 # define DDI_DLBNEXT_OP     8
 # define DDI_GDLBRESET_OP   9
 # define DDI_GDLBNEXT_OP   10

# else
 
 # define DDI_DLBRESET       7
 # define DDI_DLBNEXT        8
 # define DDI_GDLBRESET      9
 # define DDI_GDLBNEXT      10

# endif
 
 # define DDI_QUIT          11

# if defined WINTEL 
 
 # define DDI_ZERO_OP       12 
 
# else

 # define DDI_ZERO          12 

# endif
 
 # define DDI_FENCE         13

# if defined WINTEL

 # define DDI_SCATTER_ACC_OP 14

# else

 # define DDI_SCATTER_ACC   14

# endif

 # define DDI_DEBUGFLAG     20
 # define DB_CREATE_ENTRY   30
 # define DB_READ_ENTRY     31
 # define DB_WRITE_ENTRY    32

