/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Subroutines associated with the scatter accumulate operation.
 *
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"

/* -------------------------------------------------------------- *\
   DDI_Scatter_AccS(handle,scattered,alpha,ibuff,buff)
   ============================
   [IN] handle     - Handle of the distributed array to be accessed.
   [IN] scattered  - structure containing information about scattered data.
   [IN] alpha      - scale factor.
   [IN] ibuff      - Coordinates of data segment.
   [IN] buff       - Data segment to be operated on.
\* -------------------------------------------------------------- */

void DDI_Scatter_AccS(int handle,DDI_Scattered *scattered,double alpha,long *ibuff,void *buff) {
   
   /* --------------- *\
      Local Variables
   \* --------------- */
      char ack=57;
      int i,np,me,nn,my,remote_id,nsubs;
      int ranks[MAX_NODES];
      DDI_Scattered subs[MAX_NODES];
      DDI_Patch trojen;
      char *working_buffer = (char *) buff;
      char *iworking_buffer = (char *) ibuff;
      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);

# if defined DDI_LAPI
      fprintf(stdout,"%s: DDI_Scatter_Acc for DDI_LAPI not implemented.\n",DDI_Id());
      Fatal_error(911);
# endif

# if defined CRAY_MPI
      fprintf(stdout,"%s: DDI_Scatter_Acc for CRAY_MPI not implemented.\n",DDI_Id());
      Fatal_error(911);
# endif

#if defined DDI_ARMCI
      fprintf(stdout,"%s: DDI_Scatter_Acc for DDI_ARMCI not implemented.\n",DDI_Id());
      Fatal_error(911);
# endif

#if defined DDI_MPI2
      fprintf(stdout,"%s: DDI_Scatter_Acc for DDI_MPI2 not implemented.\n",DDI_Id());
      Fatal_error(911);
#endif
      
# if defined WINTEL
      fprintf(stdout,"%s: DDI_Scatter_Acc for WINTEL not implemented.\n",DDI_Id());
      Fatal_error(911);
# endif

      STD_DEBUG((stdout,"%s: Entering DDI_Scatter_AccS.\n",DDI_Id()))

   /* -------------------- *\
      Process OR Node Rank
   \* -------------------- */
      DDI_NProc(&np,&me);
      DDI_NNode(&nn,&my);


   /* ------------------------------------- *\
      Ensure the scattered has the correct info
   \* ------------------------------------- */

# if defined WINTEL
      scattered->oper   = DDI_SCATTER_ACC_OP;
# else
      scattered->oper   = DDI_SCATTER_ACC;
# endif
      scattered->handle = handle;
# if defined WINTEL
      trojen.oper = DDI_SCATTER_ACC_OP;
# else
      trojen.oper = DDI_SCATTER_ACC;
# endif

   /* ---------------------------------- *\
      Check calling arguments for errors
   \* ---------------------------------- */
    # if defined DDI_CHECK_ARGS
      if(handle < 0 || handle >= gv(ndda)) {
         fprintf(stdout,"%s: Invalid handle [%i] in DDI_Scatter_Acc.\n",DDI_Id(),handle);
         Fatal_error(911);
      }
    # endif


   /* ------------------------------ *\
      Log some simple profiling info
   \* ------------------------------ */
    # if defined DDI_COUNTERS
      gv(scatter_acc_profile).ncalls++;
      gv(scatter_acc_profile).nbytes += DDI_Patch_sizeof(scattered);
    # endif


      /* --------------------------------------------- *\
	 Sort buffers from low->high values in ibuff
      \* --------------------------------------------- */

      if(np > 1 )
      DDI_Quicksort(ibuff,buff,0,(scattered->nelem)-1);


      /* ------------------------------------------------------- *\
         Determine where the pieces of scattered data exist
      \* ------------------------------------------------------- */

      DDI_Scattered_data(handle,scattered,ibuff,&nsubs,ranks,subs);

      MAX_DEBUG((stdout,"%s: %i subscatters.\n",DDI_Id(),nsubs));

      
   /* ------------------------------------------------------------------- *\
      Send data requests for all non-local pieces of the requested scattered.
      Operate immediately to Scatter Acc a local portion of the scattered.
   \* ------------------------------------------------------------------- */
      for(i=0; i<nsubs; i++) {
         ULTRA_DEBUG((stdout,"%s: Scatter Accumulating subscatter %i.\n",DDI_Id(),i))

      /* ------------------------------------------------------------- *\
         Using SysV, take advantage of shared-memory for a local scatter
      \* ------------------------------------------------------------- */
       # if defined USE_SYSV

      /* ------------------------------------------------ *\
         Determine if the ith scattered is local to 'my' node
      \* ------------------------------------------------ */
         if(ranks[i] == my) {
            MAX_DEBUG((stdout,"%s: Subscatter %i is local.\n",DDI_Id(),i))

         /* --------------------------------------------- *\
            Perform the local scatter acc immediately.
         \* --------------------------------------------- */

	  DDI_Scatter_Acc_local(&subs[i],alpha,working_buffer,iworking_buffer);

         /* ------------------------------------------------------- *\
            Move the working buffer to the next scatter and continue.
         \* ------------------------------------------------------- */
            working_buffer += subs[i].size;
            iworking_buffer += subs[i].size;
            /* jworking_buffer += subs[i].size; */
            continue;
         }
       # endif


      /* --------------------------------- *\
         If the current scattered is NOT local 
      \* --------------------------------- */
         remote_id = ranks[i];

      /* -------------------------------- *\
         Send data request for subscatter i
      \* -------------------------------- */
         MAX_DEBUG((stdout,"%s: Sending data request to node %i.\n",
                    DDI_Id(),remote_id))


         subs[i].alpha = alpha;

	 trojen.handle = subs[i].handle;

	 DDI_Send_request(&trojen,&remote_id,NULL);

         MAX_DEBUG((stdout,"%s: data request sent to global process %i.\n",
                    DDI_Id(),remote_id))

         DDI_Scatter_Acc_remote(alpha,working_buffer,iworking_buffer,&subs[i],remote_id);

      
      /* ------------ *\
         Shift buffer 
      \* ------------ */
         working_buffer += subs[i].size;
	 iworking_buffer +=subs[i].size;

      }

      MAX_DEBUG((stdout,"%s: Leaving DDI_Scatter_AccS.\n",DDI_Id()))
      
   }


/* -------------------------------------------------------------- *\
   DDI_Scatter_Acc_local(scattered,alpha,buff,ibuff)
   =========================
   [IN] scattered  - structure containing attributed of scattered data.
   [IN] alpha      - scale factor.
   [IN] buff       - Data segment to be operated on.
   [IN] ibuff      - coordinates of data segment ot be operated on.

   Accumulates the subscatter specified by scattered and stored in buff
   into the share-memory segment(s) of the local node.
\* -------------------------------------------------------------- */
void DDI_Scatter_Acc_local(const DDI_Scattered* scattered,double alpha,void *buff,void *ibuff) {

   /* --------------- *\
      Local Variables
   \* --------------- */
    DDA_Index *Index = gv(dda_index);
    double *dda,*dloc = (double *) buff;
    long *idloc = (long *) ibuff;
    size_t i,iloc_start,iloc_end;
    int handle = scattered->handle;
    long jlo ;
    long jhi ;
    long index;
# if FULL_SMP
    int icpu,smpme,smpnp;
    DDI_SMP_NProc(&smpnp,&smpme);
# endif
    
    MAX_DEBUG((stdout,"%s: Entering DDI_Scatter_Acc_local.\n",DDI_Id()))
	
	
	iloc_end = 0;  //initialization
   /* ------------------------------------------------------------------ *\
      For FULL SMP implementations, loop on the number of SMP processors 
   \* ------------------------------------------------------------------ */
# if FULL_SMP
    for(icpu=0; icpu<smpnp; icpu++) {
	
	Index = gv(smp_index)[icpu];
	
	jlo = Index[handle].jlo;
	jhi = Index[handle].jhi;
	
	iloc_start = iloc_end;
	for(iloc_end = iloc_start; iloc_end < scattered->nlocal; iloc_end++){
	    if((idloc[iloc_end]-1) > jhi) {break;} 
	}
	
# endif
	
   /* ---------------------------------------------------------- *\
      Now Scatter-Accumulate Data
   \* ---------------------------------------------------------- */

	DDI_Acquire(Index,handle,DDI_WRITE_ACCESS,(void **) &dda);
	
	if(alpha == (double)1.0){
	    for (i = iloc_start; i < iloc_end; i++){
		index = (idloc[i]-1-jlo);
		dda[index] += dloc[i];
	    }
	}else{ 
	    for (i = iloc_start; i < iloc_end; i++){
		index = (idloc[i]-1-jlo);
		dda[index] += alpha*dloc[i];
	    }}
	
	DDI_Release(Index,handle,DDI_WRITE_ACCESS);
	
# if FULL_SMP
    } /* end for-loop on local cpus */
# endif
    
    
   /* --------------------- *\
      Shared-memory counter 
   \* --------------------- */
# if defined DDI_COUNTERS
    gv(scatter_acc_profile).ncalls_shmem++;
    gv(scatter_acc_profile).nbytes_shmem += scattered->size;
# endif
      
    MAX_DEBUG((stdout,"%s: Leaving DDI_Scatter_Acc_local.\n",DDI_Id()))
	}



/* ----------------------------------------------------------------- *\
   DDI_Scatter_Acc_server(patch,from)
   ==========================
   [IN] patch - structure containing ilo, ihi, jlo, jhi, etc.
   [IN] from  - rank of DDI process sending data to be accumulated.
   
   Used by the data server to accept incoming data and perform a
   local accumulate.  Note, the fence is raised to protect the array
   from local get/put operations until the accumulate has finished.
\* ----------------------------------------------------------------- */
   void DDI_Scatter_Acc_server(const DDI_Patch *msg,int from) {
   
   /* --------------- *\
      Local Variables
   \* --------------- */
      char ack = 57;
      double alpha;
      char *buffer = NULL;
      char *ibuffer = NULL;
      DDI_Scattered scattered;
      size_t msg_size,msg_size_max,remaining = msg->size;
      size_t scattered_size = sizeof(DDI_Scattered);
      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
      long msg_left;
   /* -------------------------------------------------------------------- *\
      Raise protective fence.  This is necessary because a compute process
      can finish with the DDI_Scatter_Acc subroutine before the remote data server
      has finished accumulating the patch.
   \* -------------------------------------------------------------------- */
# if defined USE_SYSV
      DDI_Fence_acquire(msg->handle);
      Comm_send(&ack,1,from,comm);
# endif

      /* receive information about the scattered buffer */
      Comm_recv(&scattered,scattered_size,from,comm);

      alpha = scattered.alpha;
      scattered.handle = msg->handle;

      if(scattered.size > (MAX_DS_MSG_SIZE/3)) {

	 remaining = scattered.size;
	 msg_size_max = (MAX_DS_MSG_SIZE/3);
	 msg_size = msg_size_max - (msg_size_max % sizeof(msg_size_max));

	 DDI_Memory_push(msg_size_max,(void **)&buffer,NULL);
	 DDI_Memory_push(msg_size_max,(void **)&ibuffer,NULL);

	 while(remaining){

	     msg_left = remaining-msg_size;
	     if(msg_left < 1){msg_size = remaining;}
	     
	     Comm_recv(buffer,msg_size,from,comm);
	     Comm_recv(ibuffer,msg_size,from,comm);
	     /* Comm_recv(jbuffer,msg_size,from,comm); */

	     scattered.nlocal = msg_size/sizeof(msg_size_max);
	     DDI_Scatter_Acc_local(&scattered,alpha,buffer,ibuffer);
	     
	     remaining -= msg_size;
	     if(remaining) Comm_send(&ack,1,from,comm);

	 }

	 DDI_Memory_pop(msg_size_max);
	 DDI_Memory_pop(msg_size_max);

      } else {

      	  //	 DDI_Memory_push(msg->size,(void **)&ibuffer,NULL);
      	  DDI_Memory_push(scattered.size,(void **)&buffer,NULL);
      	  DDI_Memory_push(scattered.size,(void **)&ibuffer,NULL);
      	  /* DDI_Memory_push(scattered.size,(void **)&jbuffer,NULL); */
	 
      	  //	 Comm_recv(ibuffer,msg->size,from,comm);
      	  Comm_recv(buffer,scattered.size,from,comm);
      	  Comm_recv(ibuffer,scattered.size,from,comm);
      	  /* Comm_recv(jbuffer,scattered.size,from,comm); */
	  
      	  //  	 DDI_Scatter_Acc_local(msg,alpha,buffer,ibuffer);
      	  DDI_Scatter_Acc_local(&scattered,alpha,buffer,ibuffer);
	  
      	  //     DDI_Memory_pop(msg->size);
      	  DDI_Memory_pop(scattered.size);
      	  DDI_Memory_pop(scattered.size);
      	  /* DDI_Memory_pop(scattered.size); */

      }

   /* --------------- *\
      Take down fence
   \* --------------- */
    # if defined USE_SYSV
           DDI_Fence_release(msg->handle);
    # endif

   }


void DDI_Scatter_Acc_remote(double alpha, void *buff, void *ibuff,void *scattered,int remote_id) {

      char ack = 39;
      const DDI_Scattered *msg = (const DDI_Scattered *) scattered;
      char *buffer = (char *) buff;
      char *ibuffer = (char *) ibuff;
      //      char *jbuffer = (char *) jbuff;
      DDI_Scattered *imsg = (DDI_Scattered *) scattered;
      int from = remote_id;
      size_t msg_size,msg_size_max;
      size_t remaining = msg->size;
      size_t scattered_size = sizeof(DDI_Scattered);
      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
      long msg_left;

   /* ------------------------------------------------------------- *\
    * ensure the write fence is up before sending data to be acc'ed
   \* ------------------------------------------------------------- */
    # if defined USE_SYSV
      Comm_recv(&ack,1,from,comm);
    # endif

      /* send information about the scattered object */
      Comm_send(imsg,scattered_size,from,comm);

   /* ----------------------------------------------------------------------- *\
      large message can not be fully buffed on the data server due to memory
      limitations, therefore, large message from the data server must be 
      split into smaller pieces that can be easily buffered and sent.
   \* ----------------------------------------------------------------------- */

      if(msg->size > MAX_DS_MSG_SIZE/3) {

         DEBUG_OUT(LVL6,(stdout,"%s: nr=%i; nc=%i.\n",DDI_Id(),nr,nc));
      /* --------------------------------------------------- *\
         this message will be received in multiple segements
      \* --------------------------------------------------- */

	 msg_size_max = (MAX_DS_MSG_SIZE/3);
	 msg_size = msg_size_max - (msg_size_max % sizeof(msg_size_max));

	 while(remaining){

	     msg_left = remaining-msg_size;
	     if(msg_left < 1){msg_size = remaining;}
	     
	     Comm_send(buffer,msg_size,from,comm);
	     Comm_send(ibuffer,msg_size,from,comm);
	     /* Comm_send(jbuffer,msg_size,from,comm); */
	     
	     buffer    += msg_size;
	     ibuffer   += msg_size;
	     /* jbuffer   += msg_size; */
	     
	     remaining -= msg_size;
	     if(remaining) Comm_recv(&ack,1,from,comm);
	     
	 }

      } else {

         DEBUG_OUT(LVL6,(stdout,"%s: remote_scatter_acc - vanilla - from=%i; size=%li.\n",DDI_Id(),from,msg->size))
      	     Comm_send(buffer,msg->size,from,comm);
      	     Comm_send(ibuffer,imsg->size,from,comm);
      	     /* Comm_send(jbuffer,imsg->size,from,comm); */
      }

   }
