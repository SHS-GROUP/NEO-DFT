/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * FORTRAN interface.  The subroutines in this file named F77_XXXXX are
 * converted to the proper FORTRAN external by the F77_Extern macro and
 * the defines in the ddi_fortran.h header file.
 * 
 * Author: Ryan M. Olson
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"
 # include "ddi_arr.h"
 # include "ddi_fortran.h"

/* ----------------------------- *\
   System Specific FORTRAN fixes
\* ----------------------------- */
 # if (defined CRAY)
 # define iargc_  IARGC
 # define getarg_ GETARG
 # endif 

/* ----------------------------------------- *\
   Extern subroutines from FORTRAN libraries
\* ----------------------------------------- */
   extern int iargc_();
   extern char *getarg_(int *,char *,int);


/* ---------------------------- *\
   FORTRAN Wrapper for DDI_Init
\* ---------------------------- */
   void F77_Init() {
      int i,j,lenmax=256,argc=iargc_();
      char **argv = NULL;
      char arg[256];

      STD_DEBUG((stdout," DDI: Entering F77 DDI_Init.\n"))

   /* -------------------------- *\
      Get command line arguments
   \* -------------------------- */
      if(argc) {
         argc++;
         argv = Malloc(argc*sizeof(char*));
         for(i=0; i<argc; i++) {
           for(j=0; j<256; j++) arg[j]=' ';

         # if defined CRAY
             getarg_(&i,arg,&lenmax);
         # else
             getarg_(&i,arg,lenmax);
         # endif

           for(j=0; j<256 && arg[j] != ' '; j++);
           arg[j] = 0;
           argv[i] = (char *) strdup(arg);
         }
      }

      MAX_DEBUG((stdout," DDI: Calling DDI_Init.\n"))

   /* -------------- *\
      Initialize DDI
   \* -------------- */
      DDI_Init(argc,argv);

   }


/* ---------------------------- *\
   FORTRAN Wrapper for DDI_PBeg
\* ---------------------------- */
/*   this deprecated name is just mapped onto F77_Init   */
   void F77_PBeg(int_f77 *nwdvar) {
      F77_Init();
   } 


/* -------------------------------- *\
   FORTRAN Wrapper for DDI_Finalize
\* -------------------------------- */
   void F77_Finalize() {
      DDI_Finalize();
   }


/* ---------------------------- *\
   FORTRAN Wrapper for DDI_PEnd
\* ---------------------------- */
/*   this deprecated name is just mapped onto DDI_Finalize   */
   void F77_PEnd(int_f77 *status) {
      if(*status == 1) Fatal_error(911);
      DDI_Finalize();
   }


/* ----------------------------- *\
   FORTRAN Wrapper for DDI_NProc
\* ----------------------------- */
   void F77_NProc(int_f77* np, int_f77* me) {
      int ddinp,ddime;
      DDI_NProc(&ddinp,&ddime);
      *np = (int_f77) ddinp;
      *me = (int_f77) ddime;
   }


/* ----------------------------- *\
   FORTRAN Wrapper for DDI_NNode
\* ----------------------------- */
   void F77_NNode(int_f77* np, int_f77* me) {
      int ddinp,ddime;
      DDI_NNode(&ddinp,&ddime);
      *np = (int_f77) ddinp;
      *me = (int_f77) ddime;
   }


/* ---------------------------- *\
   FORTRAN Wrapper for DDI_Send
\* ---------------------------- */
   void F77_Send(void *buff,int_f77 *size,int_f77 *to) {
      size_t isize = (size_t) *size;
      int ito      = (int) *to;
      DDI_Send(buff,isize,ito);
   }


/* ---------------------------- *\
   FORTRAN Wrapper for DDI_Recv
\* ---------------------------- */
   void F77_Recv(void *buff,int_f77 *size,int_f77* from) {
      size_t isize = (size_t) *size;
      int ifrom    = (int) *from;
      DDI_Recv(buff,isize,ifrom);
   }
 

/* ------------------------------ *\
   FORTRAN Wrapper for DDI_Memory
\* ------------------------------ */
/*     1st and 3rd argument were used by DDI version 1, 
 *     and are kept for backwards compatibility.
 *     Now, only the middle argument is needed, and should
 *     be in units of MWords, meaning 1d+6 * 8 bytes.
 */
   void F77_Memory(int_f77 *mem_rep, int_f77 *mem_ddi, int_f77 *exetyp) {
      size_t size = (size_t) *mem_ddi;
      DDI_Memory(size);
   }


/* ------------------------------ *\
   FORTRAN Wrapper for DDI_Output
\* ------------------------------ */
/*
          0 means suppress output, presently limited to ddi_create messages
          1 means print output
*/
   void F77_Output(int_f77 *value) {
      int ivalue = (int) *value;
      DDI_Output(ivalue);
   }


/* ----------------------------- *\
   FORTRAN Wrapper for DDI_Level
\* ----------------------------- */
   void F77_Level(int_f77 *implementation) {
/*
          0 means serial
          1 means distributed memory, but no subgroups
          2 means distributed memory, and subgroups (full version)
*/
      *implementation = 2;
   }


/* ------------------------------ *\
   FORTRAN Wrapper for DDI_Create
\* ------------------------------ */
   void F77_Create(int_f77 *idim, int_f77 *jdim, int_f77 *handle) {
      int ihandle;
      DDI_Create( (int) *idim, (int) *jdim, &ihandle );
      *handle = (int_f77) ihandle;
   }

   
/* ------------------------------- *\
   FORTRAN Wrapper for DDI_Destroy
\* ------------------------------- */
   void F77_Destroy(int_f77 *handle) {
      int ihandle = (int) *handle;
      DDI_Destroy(ihandle);
   }

   
/* ---------------------------- *\
   FORTRAN Wrapper for DDI_Zero
\* ---------------------------- */
/*   zero out the entire array   */
   void F77_Zero(int_f77 *handle) {
      int ihandle = (int) *handle;
      DDI_Array_zero(ihandle);
   }


/* ------------------------------- *\
   FORTRAN Wrapper for DDI_DISTRIB
\* ------------------------------- */
   void F77_Distrib(int_f77 *handle, int_f77 *rank,
                    int_f77 *ilo, int_f77 *ihi,
                    int_f77 *jlo, int_f77 *jhi) {
      DDI_Patch Patch;
      DDI_DistribP((int) *handle,(int) *rank,&Patch);
      
      *ilo = (int_f77) Patch.ilo + 1;
      *ihi = (int_f77) Patch.ihi + 1;
      *jlo = (int_f77) Patch.jlo + 1;
      *jhi = (int_f77) Patch.jhi + 1;
   }


/* -------------------------------- *\
   FORTRAN Wrapper for DDI_NDISTRIB
\* -------------------------------- */
   void F77_NDistrib(int_f77 *handle, int_f77 *rank,
                     int_f77 *ilo, int_f77 *ihi,
                     int_f77 *jlo, int_f77 *jhi) {
      DDI_Patch Patch;
      DDI_NDistribP((int) *handle,(int) *rank,&Patch);
      
      *ilo = (int_f77) Patch.ilo + 1;
      *ihi = (int_f77) Patch.ihi + 1;
      *jlo = (int_f77) Patch.jlo + 1;
      *jhi = (int_f77) Patch.jhi + 1;
   }


/* --------------------------- *\
   FORTRAN Wrapper for DDI_Get
\* --------------------------- */
   void F77_Get(int_f77 *handle, int_f77 *ilo, int_f77 *ihi,
                                 int_f77 *jlo, int_f77 *jhi, void *buff) {
      DDI_Patch Patch;
      
      Patch.oper   = DDI_GET;
      Patch.handle = (int) *handle;
      Patch.ilo    = (int) *ilo - 1;
      Patch.ihi    = (int) *ihi - 1;
      Patch.jlo    = (int) *jlo - 1;
      Patch.jhi    = (int) *jhi - 1;
      
      DDI_GetP(Patch.handle,&Patch,buff);
   }

/* -------------------------------- *\
   FORTRAN Wrapper for DDI_Get_comm
\* -------------------------------- */
   void F77_Get_comm(int_f77 *handle, int_f77 *ilo, int_f77 *ihi,
                                      int_f77 *jlo, int_f77 *jhi,
                                      void *buff, int_f77 *commid) {
      DDI_Patch Patch;
      
      Patch.oper   = DDI_GET;
      Patch.handle = (int) *handle;
      Patch.ilo    = (int) *ilo - 1;
      Patch.ihi    = (int) *ihi - 1;
      Patch.jlo    = (int) *jlo - 1;
      Patch.jhi    = (int) *jhi - 1;
      
      DDI_GetP_comm(Patch.handle,&Patch,buff,*commid);
   }

/* --------------------------- *\
   FORTRAN Wrapper for DDI_Put
\* --------------------------- */
   void F77_Put(int_f77 *handle, int_f77 *ilo, int_f77 *ihi,
                                 int_f77 *jlo, int_f77 *jhi, void *buff) {
      DDI_Patch Patch;
      
      Patch.oper   = DDI_PUT;
      Patch.handle = (int) *handle;
      Patch.ilo    = (int) *ilo - 1;
      Patch.ihi    = (int) *ihi - 1;
      Patch.jlo    = (int) *jlo - 1;
      Patch.jhi    = (int) *jhi - 1;
      
      DDI_PutP(Patch.handle,&Patch,buff);
   }

/* -------------------------------- *\
   FORTRAN Wrapper for DDI_Put_comm
\* -------------------------------- */
   void F77_Put_comm(int_f77 *handle, int_f77 *ilo, int_f77 *ihi,
                                      int_f77 *jlo, int_f77 *jhi,
                                      void *buff, int_f77 *commid) {
      DDI_Patch Patch;
      
      Patch.oper   = DDI_PUT;
      Patch.handle = (int) *handle;
      Patch.ilo    = (int) *ilo - 1;
      Patch.ihi    = (int) *ihi - 1;
      Patch.jlo    = (int) *jlo - 1;
      Patch.jhi    = (int) *jhi - 1;
      
      DDI_PutP_comm(Patch.handle,&Patch,buff,*commid);
   }

/* --------------------------- *\
   FORTRAN Wrapper for DDI_Acc
\* --------------------------- */
   void F77_Acc(int_f77 *handle, int_f77 *ilo, int_f77 *ihi,
                                 int_f77 *jlo, int_f77 *jhi, void *buff) {
      DDI_Patch Patch;
      double alpha = (double)1.0;
      
      Patch.oper   = DDI_ACC;
      Patch.handle = (int) *handle;
      Patch.ilo    = (int) *ilo - 1;
      Patch.ihi    = (int) *ihi - 1;
      Patch.jlo    = (int) *jlo - 1;
      Patch.jhi    = (int) *jhi - 1;
      
      DDI_AccP(Patch.handle,&Patch,alpha,buff);
   }


/* -------------------------------- *\
   FORTRAN Wrapper for DDI_DLBReset
\* -------------------------------- */
   void F77_DLBReset() {
      DDI_DLBReset();
   }

      
/* ------------------------------- *\
   FORTRAN Wrapper for DDI_DLBNext
\* ------------------------------- */
   void F77_DLBNext(int_f77 *counter) {
      size_t icounter;
      DDI_DLBNext(&icounter);
      *counter = (int_f77) icounter;
   }


/* ----------------------------- *\
   FORTRAN Wrapper for DDI_ISend
\* ----------------------------- */
   void F77_ISend(void *buffer, int_f77 *size, int_f77 *to, int_f77 *req) {
      int ireq = 0; 
      int ito = (int) *to;
      size_t isize = (size_t) *size;
      DDI_ISend(buffer,isize,ito,&ireq);
      *req = (int_f77) ireq;
   }


/* ----------------------------- *\
   FORTRAN Wrapper for DDI_IRecv
\* ----------------------------- */
   void F77_IRecv(void *buffer, int_f77 *size, int_f77 *from, int_f77 *req) {
      int ireq = 0;
      int ifrom = (int) *from;
      size_t isize = (size_t) *size;
      DDI_IRecv(buffer,isize,ifrom,&ireq);
      *req = (int_f77) ireq;
   }


/* ---------------------------- *\
   FORTRAN Wrapper for DDI_Wait
\* ---------------------------- */
   void F77_Wait(int_f77 *req) {
      int ireq = (int) *req;
      DDI_Wait(ireq);
   }


/* ----------------------------- *\
   FORTRAN Wrapper for DDI_BCast
\* ----------------------------- */
   void F77_BCast(int_f77* msgtag, char *type,
                  void *buffer, int_f77 *len, int_f77 *root) {
      char t = *type;
      size_t size;

      if(t == 'i' || t == 'I')      size = (*len)*sizeof(int_f77);
      else if(t == 'f' || t == 'F') size = (*len)*sizeof(double);
      else {
         fprintf(stdout,"%s: Unknown type passed in DDI_BCast\n",DDI_Id());
         Fatal_error(911);
      }

      DDI_BCast(buffer,size,(int)*root);
   }


/* ----------------------------- *\
   FORTRAN Wrapper for DDI_GSumI
\* ----------------------------- */
   void F77_GSumI(int_f77 *msgtag, void *buffer, int_f77 *len) {
      DDI_GSum(buffer,*len,sizeof(int_f77));
   } 


/* ----------------------------- *\
   FORTRAN Wrapper for DDI_GSumF
\* ----------------------------- */
   void F77_GSumF(int_f77 *msgtag, void *buffer, int_f77 *len) {
      DDI_GSum(buffer,*len,0);
   }


/* ---------------------------- *\
   FORTRAN Wrapper for DDI_Sync
\* ---------------------------- */
   void F77_Sync(int_f77 *msgtag) {
      int tag = (int) *msgtag;
      DDI_Sync(tag);
   }


/* ----------------------------- *\
   FORTRAN Wrapper for DDI_Debug
\* ----------------------------- */
   void F77_Debug(int_f77 *flag) {
      int value = (int) *flag;
      DebugOutput(value);
   }


/* ------------------------------------ *\
   FORTRAN Wrapper for DDI_Group_Create
\* ------------------------------------ */
   void F77_Group_create(int_f77 *ngroups, int_f77 *world,
                         int_f77 *group, int_f77 *master) {
      int handle,type=0;
      int ngrp = (int) *ngroups;
      int iw,ig,im;
      DDI_Group_create(ngrp,&iw,&ig,&im);
      *world  = iw;
      *group  = ig;
      *master = im;
   }


/* ----------------------------- *\
   FORTRAN Wrapper for DDI_Scope
\* ----------------------------- */
   void F77_Scope(int_f77 *scope) {
      int iscope = (int) *scope;
      DDI_Scope(iscope);
   }


/* ------------------------------ *\
   FORTRAN Wrapper for DDI_AScope
\* ------------------------------ */
   void F77_AScope(int_f77 *scope) {
      int iscope = (int) *scope;
      DDI_AScope(iscope);
   }


/* --------------------------------- *\
   FORTRAN Wrapper for DDI_GDLBReset
\* --------------------------------- */
   void F77_GDLBReset() {
      DDI_GDLBReset();
   }


/* -------------------------------- *\
   FORTRAN Wrapper for DDI_GDLBNext
\* -------------------------------- */
   void F77_GDLBNext(int_f77 *counter) {
      size_t value;
      DDI_GDLBNext(&value);
      *counter = (int_f77) value;
   }


/* ----------------------------- *\
   FORTRAN Wrapper for DDI_Ngroup
\* ----------------------------- */
   void F77_NGroup(int_f77 *ngroups, int_f77 *mygroup) {
      int ngr, mygr;
   
      DDI_NGroup(&ngr,&mygr);
      *ngroups = (int_f77) ngr;
      *mygroup = (int_f77) mygr;
   }


/* ------------------------------------- *\
   FORTRAN Wrapper for DDI_Create_custom
\* ------------------------------------- */
   void F77_Create_custom(int_f77 *idim, int_f77 *jdim,
                          int_f77 *jcols, int_f77 *handle) {
      int ihandle,i,np,me;
      int ijcols[MAX_PROCESSORS];
      DDI_NProc(&np,&me);
      for(i=0; i<np; i++) ijcols[i] = (int) *jcols;
      DDI_Create_custom( (int) *idim, (int) *jdim, ijcols, &ihandle );
      *handle = (int_f77) ihandle;
   }


/* ------------------------------------------- *\
   FORTRAN Wrapper for DDI_Group_create_custom
\* ------------------------------------------- */
   void F77_Group_create_custom(int_f77 *ngroups, int_f77 *group_size,
                                int_f77 *world, int_f77 *group,
                                int_f77 *master) {
      int i,handle,type=0;
      int ngrp = (int) *ngroups;
      int grpsize[MAX_NODES];
      int iw,ig,im;
      for (i=0;i<ngrp;i++) grpsize[i]=(int)group_size[i];
      DDI_Group_create_custom(ngrp,grpsize,&iw,&ig,&im);
      *world  = iw;
      *group  = ig;
      *master = im;
   }


/* ----------------------------------- *\
   FORTRAN Wrapper for DDI_Timer_reset
\* ----------------------------------- */
   void F77_Timer_reset() {
      DDI_Timer_reset();
   }


/* ------------------------------------ *\
   FORTRAN Wrapper for DDI_Timer_output
\* ------------------------------------ */
   void F77_Timer_output() {
      DDI_Timer_output();
   }


/* -------------------------------------- *\
   FORTRAN Wrapper for DDI_ProcDLB_create
\* -------------------------------------- */
   void F77_ProcDLB_create(int_f77 *handle) {
      int hnd;
      DDI_ProcDLB_create(&hnd);
      *handle = (int_f77) hnd;
   }


/* --------------------------------------- *\
   FORTRAN Wrapper for DDI_ProcDLB_destroy
\* --------------------------------------- */
   void F77_ProcDLB_destroy(int_f77 *handle) {
      int hnd = (int) *handle;
      DDI_ProcDLB_destroy(hnd);
   }


/* ------------------------------------- *\
   FORTRAN Wrapper for DDI_ProcDLB_reset
\* ------------------------------------- */
   void F77_ProcDLB_reset(int_f77 *handle) {
      int hnd = (int) *handle;
      DDI_ProcDLB_reset(hnd);
   }


/* ------------------------------------ *\
   FORTRAN Wrapper for DDI_ProcDLB_next
\* ------------------------------------ */
   void F77_ProcDLB_next(int_f77 *handle, int_f77 *proc, int_f77 *counter) {
      int hnd = (int) *handle;
      int cpu = (int) *proc;
      int val = 0;
      DDI_ProcDLB_next(hnd,cpu,&val);
      *counter = (int_f77) val;
   }


/* --------------------------------- *\
   FORTRAN Wrapper for DDI_SMP_NProc
\* --------------------------------- */
   void F77_SMP_NProc(int_f77* np, int_f77* me) {
      int ddinp,ddime;
      DDI_SMP_NProc(&ddinp,&ddime);
      *np = (int_f77) ddinp;
      *me = (int_f77) ddime;
   }


/* --------------------------------- *\
   FORTRAN Wrapper for DDI_SMP_BCast
\* --------------------------------- */
   void F77_SMP_sync() {
      DDI_Comm *comm = Comm_find(DDI_COMM_WORLD);
      Comm_sync_smp(comm); 
   }


/* --------------------------------- *\
   FORTRAN Wrapper for DDI_SMP_BCast
\* --------------------------------- */
   void F77_BCast_smp(int_f77* msgtag, char *type, void *buffer,
                      int_f77 *len, int_f77 *root) {
      char t = *type;
      size_t size;

      if(t == 'i' || t == 'I')      size = (*len)*sizeof(int_f77);
      else if(t == 'f' || t == 'F') size = (*len)*sizeof(double);
      else {
         fprintf(stdout,"%s: Unknown type passed in DDI_BCast\n",DDI_Id());
         Fatal_error(911);
      }

      DDI_BCast_smp(buffer,size,*root);
   }


/* --------------------------------- *\
   FORTRAN Wrapper for DDI_SMP_GSumI
\* --------------------------------- */
   void F77_GSumI_smp(int_f77 *msgtag, void *buffer, int_f77 *len) {
      DDI_GSum_smp(buffer,(size_t)*len,sizeof(int_f77));
   } 


/* --------------------------------- *\
   FORTRAN Wrapper for DDI_SMP_GSumF
\* --------------------------------- */
   void F77_GSumF_smp(int_f77 *msgtag, void *buffer, int_f77 *len) {
      DDI_GSum_smp(buffer,(size_t)*len,0);
   }


/* ------------------------------------- *\
   FORTRAN Wrapper for DDI_Masters_gsumi
\* ------------------------------------- */
   void F77_GSumI_node(int_f77 *msgtag, void *buffer, int_f77 *len) {
      DDI_GSum_node(buffer,(size_t) *len,sizeof(int_f77) /*type*/);
   }


/* ------------------------------------- *\
   FORTRAN Wrapper for DDI_Masters_gsumf
\* ------------------------------------- */
   void F77_GSumF_node(int_f77 *msgtag, void *buffer, int_f77 *len) {
      DDI_GSum_node(buffer,(size_t) *len,0 /*type*/);
   }


/* ------------------------------------- *\
   FORTRAN Wrapper for DDI_Masters_bcast
\* ------------------------------------- */
   void F77_BCast_node(int_f77* msgtag, char *type, void *buffer,
                       int_f77 *len, int_f77 *root) {
      char t = *type;
      size_t size;

      if(t == 'i' || t == 'I')      size = (*len)*sizeof(int_f77);
      else if(t == 'f' || t == 'F') size = (*len)*sizeof(double);
      else {
         fprintf(stdout,"%s: Unknown type passed in DDI_BCast\n",DDI_Id());
         Fatal_error(911);
      }

      DDI_BCast_node(buffer,size,(int)*root);
   }

/* ---------------------- *\
   DDI_SMP_Create wrapper
\* ---------------------- */
   void F77_SMP_Create(int_f77 *size, int_f77 *handle) {
      int myhndl = 0;
      size_t mysize = (size_t) *size;
      mysize *= sizeof(double);
      *handle = (int_f77) SMP_create(mysize);
   }


/* ----------------------- *\
   DDI_SMP_Destroy wrapper
\* ----------------------- */
   void F77_SMP_Destroy(int_f77 *handle) {
      DDI_SMP_Destroy((int) *handle);
   }


/* ---------------------- *\
   DDI_SMP_Offset wrapper
\* ---------------------- */
   void F77_SMP_Offset(int_f77 *handle, void *addr, int_f77 *offset) {

      unsigned INT_SIZE loc_addr_val = (unsigned INT_SIZE) addr;
      unsigned INT_SIZE smp_addr_val = (unsigned INT_SIZE) DDI_SMP_Array(*handle);

      *offset = smp_addr_val - loc_addr_val;
/*
      fprintf(stdout,"%s: SMP HANDLE %i = %x\n",DDI_Id(),
              *handle,DDI_SMP_Addr(*handle));
      fflush(stdout);
*/
/*
      if(*offset % 8) {
         printf("Warning, SMP OFFSET not one WORD boundaries!!!\n");
      }
*/
   }


/* --------------------- *\
   DDI_Addr_test wrapper
\* --------------------- */
   void F77_Addr_test(void *addr) {
      fprintf(stdout,"%s: addr=%x\n",DDI_Id(),addr);
      fflush(stdout);
   }


/*    FORTRAN bindings for various array operations     */

void F77_ARR_zero(int_f77 *dA, int_f77 *dAi, int_f77 *dAi2,
                               int_f77 *dAj, int_f77 *dAj2)
{
  DDI_ARR_zero((int)*dA, (int)*dAi - 1, (int)*dAi2 - 1,
                         (int)*dAj - 1, (int)*dAj2 - 1);
}

void F77_ARR_fill(int_f77 *dA, int_f77 *dAi, int_f77 *dAi2,
                               int_f77 *dAj, int_f77 *dAj2, void *alpha)
{
  DDI_ARR_fill((int)*dA, (int)*dAi - 1, (int)*dAi2 - 1,
                         (int)*dAj - 1, (int)*dAj2 - 1, *(double*)alpha);
}

void F77_ARR_scale(int_f77 *dA, int_f77 *dAi, int_f77 *dAi2,
                                int_f77 *dAj, int_f77 *dAj2, void *alpha)
{
  DDI_ARR_scale((int)*dA, (int)*dAi - 1, (int)*dAi2 - 1,
                          (int)*dAj - 1, (int)*dAj2 - 1, *(double*)alpha);
}

void F77_ARR_min(int_f77 *dA, int_f77 *dAi, int_f77 *dAi2,
                              int_f77 *dAj, int_f77 *dAj2,
		 void *alpha, int_f77 *i, int_f77 *j)
{
  int index[2];
  DDI_ARR_min((int)*dA, (int)*dAi - 1, (int)*dAi2 - 1,
                        (int)*dAj - 1, (int)*dAj2 - 1, (double*)alpha, index);
  *i = (int_f77)index[0] + 1;
  *j = (int_f77)index[1] + 1;
}

void F77_ARR_max(int_f77 *dA, int_f77 *dAi, int_f77 *dAi2,
                              int_f77 *dAj, int_f77 *dAj2,
		 void *alpha, int_f77 *i, int_f77 *j)
{
  int index[2];
  DDI_ARR_max((int)*dA, (int)*dAi - 1, (int)*dAi2 - 1,
                        (int)*dAj - 1, (int)*dAj2 - 1, (double*)alpha, index);
  *i = (int_f77)index[0] + 1;
  *j = (int_f77)index[1] + 1;
}

void F77_ARR_dot(int_f77 *dA, int_f77 *dAi, int_f77 *dAi2,
                              int_f77 *dAj, int_f77 *dAj2,
		 int_f77 *dB, int_f77 *dBi, int_f77 *dBi2,
                              int_f77 *dBj, int_f77 *dBj2, void *x) 
{
  DDI_ARR_dot((int)*dA, (int)*dAi - 1, (int)*dAi2 - 1,
                        (int)*dAj - 1, (int)*dAj2 - 1,
	      (int)*dB, (int)*dBi - 1, (int)*dBi2 - 1,
                        (int)*dBj - 1, (int)*dBj2 - 1, x);
}

void F77_ARR_add(int_f77 *dA, int_f77 *dAi, int_f77 *dAi2,
                              int_f77 *dAj, int_f77 *dAj2, void *alpha,
		 int_f77 *dB, int_f77 *dBi, int_f77 *dBi2,
                              int_f77 *dBj, int_f77 *dBj2, void *beta,
		 int_f77 *dC, int_f77 *dCi, int_f77 *dCi2,
                              int_f77 *dCj, int_f77 *dCj2)
{
  DDI_ARR_add((int)*dA, (int)*dAi - 1, (int)*dAi2 - 1,
                        (int)*dAj - 1, (int)*dAj2 - 1, *(double*)alpha,
	      (int)*dB, (int)*dBi - 1, (int)*dBi2 - 1,
                        (int)*dBj - 1, (int)*dBj2 - 1, *(double*)beta,
	      (int)*dC, (int)*dCi - 1, (int)*dCi2 - 1,
                        (int)*dCj - 1, (int)*dCj2 - 1);
}

void F77_ARR_acc(int_f77 *dA, int_f77 *dAi, int_f77 *dAi2,
                              int_f77 *dAj, int_f77 *dAj2,
                              void *alpha, void *buf) 
{
  DDI_ARR_acc((int)*dA, (int)*dAi - 1, (int)*dAi2 - 1,
                        (int)*dAj - 1, (int)*dAj2 - 1, *(double*)alpha, buf);
}
