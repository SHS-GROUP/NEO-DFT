/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Collective operations used to create a distributed array.
 *
 * Author: Ryan M. Olson
 * 10 Jun 09 - RMO - allow for one data server/node on Cray
 * 29 Mar 10 - SS  - add missing 3rd argument to ddi-send-request
 * 18 Aug 10 - MWS - use 64 bit integer for printing of DM total size
 * 18 Feb 13 - MWS - trap for >2 Gword MEMDDI slice size
 *  2 Aug 13 - MWS - fix slice size check for TCP/IP runs
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"

/* -------------------------------------------------------------------- *\
   DDI_Create(idim,jdim,handle)
   ============================
   [IN]  idim   - Number of rows in the array to be created.
   [IN]  jdim   - Number of columns in the array to be created.
   [OUT] handle - Handle given to the newly created array.
   
   Creates a distributed array with the columns evenly divided amongst
   the processors.
\* -------------------------------------------------------------------- */
   void DDI_Create(int idim,int jdim,int *handle) {

   /* --------------- *\
      Local Variables
   \* --------------- */
      int i,np,me;
      int icol,mincol,lftcol;
      int jcols[MAX_PROCESSORS];
      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);

      DEBUG_ROOT(LVL1,(stdout," DDI: Entering DDI_Create.\n"))
      DEBUG_OUT(LVL3,(stdout,"%s: Entering DDI_Create.\n",DDI_Id()))
     
      np = comm->np;
      me = comm->me;
/*      
      if(jdim < np && me == 0) {
         fprintf(stdout," DDI Error: Trying to create an array with fewer columns than processors.\n");
         fprintf(stdout," DDI Error: Reduce the number of processors and try again.\n");
         Fatal_error(911);
      }
 */           
      mincol = jdim / np;
      lftcol = jdim % np;
      
      for(i=0,icol=0; i<np; i++) {
         jcols[i] = icol;
         icol += mincol;
         if(i<lftcol) icol++;
      }
      
      DDI_Create_custom(idim,jdim,jcols,handle);
      
      DEBUG_ROOT(LVL2,(stdout," DDI: Array[%i] successfully created.\n",*handle)) 
   }


/* -------------------------------------------------------------------- *\
   DDI_Create_custom(idim,jdim,jcols,handle)
   =========================================
   [IN]  idim   - Number of rows in the array to be created.
   [IN]  jdim   - Number of columns in the array to be created.
   [IN]  jcols  - Array holding the number of columns to be given to
                - each processor when creating the distributed array.
   [OUT] handle - Handle given to the newly created array.
   
   Creates a distributed array where the user can customize how the
   array is distributed across the processors.
\* -------------------------------------------------------------------- */
   void DDI_Create_custom(int idim,int jdim,int *jcols,int *handle) {
   
      int i,np,me,nn,my;
      int inode;
      DDI_INT64 totwrds;
      DDI_INT64 longrows,longcols,longslice,longnd,long2g;
    # ifndef USE_SYSV
      int remote_id;
    # endif
      DDI_Patch patch;
      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
      
      np = comm->np;
      me = comm->me;
      nn = comm->nn;
      my = comm->my;

      Comm_sync(3001,comm);

      *handle   = gv(ndda)++;
     
    # ifndef USE_SYSV
      remote_id = my;
    # endif

      DEBUG_ROOT(LVL2,(stdout," DDI: Entering DDI_Create_custom.\n"))
      DEBUG_ROOT(LVL2,(stdout," DDI: Creating Array [%i] - %ix%i=%i.\n",*handle,idim,jdim,idim*jdim))
      DEBUG_OUT(LVL3,(stdout,"%s: Entering DDI_Create_custom.\n",DDI_Id()))

    # ifdef DS_SIGNAL
      if(comm->me_local == 1) {
         signal(SIGALRM,DS_Thread_main);
      }
    # endif
      
      if(me == 0) {
         if(gv(dda_output)) {
            longrows = idim;
            longcols = jdim;
            totwrds = longrows*longcols;
            fprintf(stdout," DDI: Creating Array [%i] - %i x %i = %li words.\n",
                                    *handle,idim,jdim,totwrds);
            fflush(stdout);
         }
      }

   /*
       Make sure each slice of the distributed array will be under 2 GWords.

       Even on 64-bit hardware, most counting in this program is done
       with 32-bit data types, meaning we can't count higher than 2**31-1.

       If on 32-bit hardware, the 'long' data types here will be 32-bits,
       and so we'll see crazy products, including less than zero.
       In present form, nothing will be trapped here on a 32 bit machine!
   */
      longrows  = idim;
      longcols  = jdim;
      totwrds   = longrows*longcols;
   /*     Total distributed array over 2 Gwords is OK, but each  */
   /*     slice (MEMDDI per data server) must be under 2 GWords. */
   /*     TCP/IP has gv(nd)=-1 (uninitialized)                   */
   /*     Cray on one node has gv(nd)=0 since no d.s. exists.    */
      # if defined DDI_MPI
         longnd    = gv(nd);
         if (longnd <= 0) longnd=1;
      # endif
      # if defined DDI_SOC
         longnd = np;
      # endif
      longslice = totwrds/longnd;
   /*  next is largest signed 32 bit integer, stored as 64 bit quantity  */
      long2g   = 2147483643;
      if (longslice > long2g)
         {
            fprintf(stdout,"\n");
            fprintf(stdout," DDI: trouble creating distributed array!\n");
            fprintf(stdout," Current number of data servers is %li\n",longnd);
            fprintf(stdout," so each data server's slice of array");
            fprintf(stdout," [%i] is %li words\n",*handle,longslice);
            fprintf(stdout,"\n");
            fprintf(stdout," Add more processors so required total array");
            fprintf(stdout," size %li words\n",totwrds);
            fprintf(stdout," divided by no. of processors (data servers)");
            fprintf(stdout," is less than 2 Gwords= %li\n",long2g);
            fprintf(stdout," For example, %li or more data servers...\n",
                                1+totwrds/long2g);
            fprintf(stdout,"\n");
            fflush(stdout);
            Fatal_error(911);
         }

   /* ------------------------------------ *\
      Ensure 'jcols' is properly formatted
   \* ------------------------------------ */
      for(i=0; i<np; i++) {
         if(jcols[i] < 0 && me == 0) {
            fprintf(stdout," Error in argument 3 of DDI_Create_custom: Values must be >= 0.\n");
            Fatal_error(911);
         }
         
         if(i > 0)
         if(jcols[i] < jcols[i-1]) {
            fprintf(stdout," Error in argument 3 of DDI_Create_custom: Values must increase monotonically.\n");
            Fatal_error(911);
         }
      }
   
   /* ----------------------------------------------------------------- *\
      Check to ensure the maximum number of arrays hasn't been reached.
   \* ----------------------------------------------------------------- */
      if( gv(ndda) == MAX_DD_ARRAYS ) {
        if(me == 0) {
           fprintf(stderr," DDI Error:  The maximum number of distributed arrays [%i] has been reached.\n",MAX_DD_ARRAYS);
           fprintf(stderr," Information:  The maximum number of distributed arrays is a DDI compile-time option.\n");
        }
        Fatal_error(911);
      }

      gv(nrow)[*handle] = idim;
      gv(ncol)[*handle] = jdim;
      
 
   /* ---------------------------------------------------- *\
      Generate Column Mapping by Compute Process & by Node
   \* ---------------------------------------------------- */
      for(i=0,inode=-1; i<np; i++) {
        gv(pcmap)[*handle][i] = jcols[i];

     /* if(inode == gv(ddiprocs)[i].node) continue; */
        if(inode == comm->local_nid[i]) continue;
        gv(ncmap)[*handle][++inode] = gv(pcmap)[*handle][i];
      }

      gv(pcmap)[*handle][np] = jdim;
      gv(ncmap)[*handle][nn] = jdim;


   /* -------------------------- *\
      Get local patch dimensions
   \* -------------------------- */
      DDI_DistribP(*handle,me,&patch);
 
   /* ----------------------------- *\
      Create Distributed Data Array
   \* ----------------------------- */
      patch.handle = *handle;
# if defined WINTEL
      patch.oper   = DDI_CREATE_OP;
# else
      patch.oper   = DDI_CREATE;
# endif
      patch.size   = jdim;


# if defined USE_SYSV || defined DDI_ARMCI || defined DDI_MPI2
      DDI_Index_create(&patch);
# else
      DDI_Send_request(&patch,&remote_id,NULL);
# endif 


   /* ----------------------------- *\
      Synchronize Compute Processes
   \* ----------------------------- */
      Comm_sync(3002,comm);

      DEBUG_OUT(LVL3,(stdout,"%s: Leaving DDI_Create_custom.\n",DDI_Id()))
   }


   void DDI_Output(int flag) {
      gv(dda_output) = flag;
   }
