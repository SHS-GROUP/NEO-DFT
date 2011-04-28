/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Data Server Subroutine
 * 
 * Author: Ryan M. Olson
 * 10 Jun 09 - RMO - allow for one data server/node on Cray
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"

 # ifdef CRAY_MPI
   static int NRequests;
   static int *index_ds = NULL;
   static DDI_Patch *p = NULL;
   static MPI_Status *s = NULL;
   static MPI_Request *r = NULL;
 # endif

/* -------------------------------------------------------------------- *\
   DDI_Server()
   ============
   
   Called by DDI processes that specialize to become data servers.
\* -------------------------------------------------------------------- */
   void DDI_Server() {
   
   /* --------------- *\
      Local Variables
   \* --------------- */
      int from;
      char ack=57;
      char server=1;
      DDI_Patch *msg;
      DDI_Patch patch;
      size_t counter_value = 0;
      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_COMM_WORLD);
 
      DEBUG_OUT(LVL2,(stdout,"%s: Starting DDI data server.\n",DDI_Id()))

    # ifdef CRAY_MPI
      int i;
      int nfinished =  0;
      int last      = -1;
      int size      = sizeof(DDI_Patch);
      index_ds = (int *) Malloc(comm->np*sizeof(int));
      p = (DDI_Patch *) Malloc(comm->np*sizeof(DDI_Patch));
      s = (MPI_Status *) Malloc(comm->np*sizeof(MPI_Status));
      r = (MPI_Request *) Malloc(comm->np*sizeof(MPI_Request));

   /* ----------------------------------------------------------- *\
      Post IRecvs for remote data requests from all the processes
   \* ----------------------------------------------------------- */
      DEBUG_OUT(LVL2,(stdout,
          "%s: (DS) Posting MPI_IRecvs for data requests.\n",DDI_Id()))
      for(i=0; i<comm->np; i++) {
         MPI_Irecv(&p[i],size,MPI_BYTE,i,37,comm->world_comm,&r[i]);
      }
      NRequests = comm->np;
    # else

      int np,me,nn,my;
   /* ------------------------------------------------ *\
      Discovery Identity -- If it only were this easy!
   \* ------------------------------------------------ */
      DDI_NProc(&np,&me);
      DDI_NNode(&nn,&my);
/*
      if(me == np) DB_Init();
*/
    # endif

   /* -------------------- *\
      DDI Data Server Loop
   \* -------------------- */
      do {
  
       # ifdef CRAY_MPI
         MPI_Testsome(NRequests,r,&nfinished,index_ds,s);
         for(i=0; i<nfinished; i++) {
            msg = &p[index_ds[i]];
            from  = s[i].MPI_SOURCE;
       # else
         DDI_Recv_request(&patch,&from);
         msg = (DDI_Patch *) &patch;
       # endif
   
         switch(msg->oper) {

           case DDI_DEBUGFLAG:
              DebugOutput(msg->handle);
              break;

   
           case DDI_MEMORY:
              DEBUG_OUT(LVL2,(stdout,
                 "%s: (DS) Received DDI_MEMORY request.\n",DDI_Id()))
   
              DDI_Memory_server(msg->size);
              Comm_send(&ack,1,from,comm);
 
              DEBUG_OUT(LVL3,(stdout,
                 "%s: (DS) DDI_MEMORY requested completed.\n",DDI_Id()))
              break;
   
   
           case DDI_CREATE:
              DEBUG_OUT(LVL2,(stdout,
                 "%s: (DS) Received DDI_CREATE[%i] request.\n",
                 DDI_Id(),msg->handle))
   
              DDI_Index_create(msg);
   
              DEBUG_OUT(LVL3,(stdout,
                 "%s: (DS) Array[%i] created successfully.\n",
                 DDI_Id(),msg->handle))
              break;
   
   
           case DDI_DESTROY:
              DEBUG_OUT(LVL2,(stdout,
                 "%s: (DS) Received DDI_DESTROY[%i] request.\n",
                 DDI_Id(),msg->handle))

              DDI_Index_destroy(msg); 

              DEBUG_OUT(LVL3,(stdout,
                 "%s: (DS) Array[%i] destroyed successfully.\n",
                 DDI_Id(),msg->handle))
              break;
   
   
           case DDI_ZERO:
              DEBUG_OUT(LVL2,(stdout,
                "%s: (DS) Received DDI_ZERO request from %i.\n",DDI_Id(),from))
              DDI_Array_zero(msg->handle);
              DEBUG_OUT(LVL3,(stdout,
                "%s: (DS) Finished DDI_ZERO request from %i.\n",DDI_Id(),from))
              break;
   
   
           case DDI_GET:
              DEBUG_OUT(LVL2,(stdout,
                 "%s: (DS) Received DDI_GET request from %i.\n",DDI_Id(),from))
              DDI_Get_server(msg,from);
              DEBUG_OUT(LVL3,(stdout,
                 "%s: (DS) Finished DDI_GET request from %i.\n",DDI_Id(),from))
              break;
   
           
           case DDI_PUT:
              DEBUG_OUT(LVL2,(stdout,
                 "%s: (DS) Received DDI_PUT request from %i.\n",DDI_Id(),from))
              DDI_Put_server(msg,from);
              DEBUG_OUT(LVL3,(stdout,
                 "%s: (DS) Finished DDI_PUT request from %i.\n",DDI_Id(),from))
              break;
   
   
           case DDI_ACC:
              DEBUG_OUT(LVL2,(stdout,
                 "%s: (DS) Received DDI_ACC request from %i.\n",DDI_Id(),from))
              DDI_Acc_server(msg,from);
              DEBUG_OUT(LVL3,(stdout,
                 "%s: (DS) Finished DDI_ACC request from %i.\n",DDI_Id(),from))
              break;
              
              
           case DDI_GETACC:
              DEBUG_OUT(LVL2,(stdout,
                 "%s: (DS) Received DDI_GETACC request from %i.\n",
                 DDI_Id(),from))
              DDI_GetAcc_server(msg,from);
              DEBUG_OUT(LVL3,(stdout,
                 "%s: (DS) Finished DDI_GETACC request from %i.\n",
                 DDI_Id(),from))
              break;
              
              
           case DDI_DLBRESET:
              DEBUG_OUT(LVL2,(stdout,
                 "%s: (DS) Received DDI_DLBRESET request from %i.\n",
                 DDI_Id(),from))
              DDI_DLBReset_local();
              DEBUG_OUT(LVL3,(stdout,
                 "%s: (DS) Finished DDI_DLBRESET request from %i.\n",
                 DDI_Id(),from))
              break;
              break;
   
    
           case DDI_DLBNEXT:
              DEBUG_OUT(LVL2,(stdout,
                 "%s: (DS) Received DDI_DLBNEXT request from %i.\n",
                 DDI_Id(),from))
              DDI_DLBNext_local(&counter_value);
              Comm_send(&counter_value,sizeof(size_t),from,comm);
              DEBUG_OUT(LVL3,(stdout,
                 "%s: (DS) Finished DDI_DLBNEXT request from %i.\n",
                 DDI_Id(),from))
              break;


           case DDI_GDLBRESET:
              DDI_GDLBReset_local();
              break;
      
      
           case DDI_GDLBNEXT: 
              DDI_GDLBNext_local(&counter_value);
              DDI_Send(&counter_value,sizeof(size_t),from);
              break;

   
           case DDI_QUIT: /* Quit server loop, synchronize, and exit */
             DEBUG_OUT(LVL3,(stdout,
                 "%s: (DS) Received DDI_QUIT request\n",DDI_Id()))
          /* if(me == np) DB_Close(); */
             DDI_Memory_finalize(); 
             Comm_send(&ack,1,from,comm);
             server=0;
             break;

   
          /* --------------------------------------------- *\
             Clean-up distributed-memory & check for leaks
          \* --------------------------------------------- */
/*
           case DB_CREATE_ENTRY:
             DEBUG_OUT(LVL3,(stdout,
                   "%s: (DS) Received DB_CREATE_ENTRY request.\n",DDI_Id()))
             if(me != np) {
                fprintf(stdout,
                   "%s: received DB request but is not master data server.\n",
                   DDI_Id());
                Fatal_error(911);
             }
             DB_Create_server(msg,from);
             DEBUG_OUT(LVL3,(stdout,
                   "%s: (DS) Finished DB_CREATE_ENTRY request.\n",DDI_Id()))
             break;

           case DB_READ_ENTRY:
             DEBUG_OUT(LVL3,(stdout,
                   "%s: (DS) Received DB_READ_ENTRY request.\n",DDI_Id()))
             if(me != np) {
                fprintf(stdout,
                   "%s: received DB request but is not master data server.\n",
                   DDI_Id());
                Fatal_error(911);
             }
             DB_Read_server(msg,from);
             DEBUG_OUT(LVL3,(stdout,
                   "%s: (DS) Finished DB_READ_ENTRY request.\n",DDI_Id()))
             break;

           case DB_WRITE_ENTRY:
             DEBUG_OUT(LVL3,(stdout,
                   "%s: (DS) Received DB_WRITE_ENTRY request.\n",DDI_Id()))
             if(me != np) {
                fprintf(stdout,
                   "%s: received DB request but is not master data server.\n",
                   DDI_Id());
                Fatal_error(911);
             }
             DB_Write_server(msg,from);
             DEBUG_OUT(LVL3,(stdout,
                   "%s: (DS) Finished DB_WRITE_ENTRY request.\n",DDI_Id()))
             break;
*/

	 case DDI_ARR_ZERO:
	 case DDI_ARR_FILL:
	 case DDI_ARR_SCALE:
	   DDI_ARR_scalar_server(msg, from);
	   break;
	   
	 case DDI_ARR_MIN:
	 case DDI_ARR_MAX:
	   DDI_ARR_select_server(msg, from);
	   break;

	 case DDI_ARR_DOT:
	   DDI_ARR_dot_server(msg, from);
	   break;

	 case DDI_ARR_ADD:
	   DDI_ARR_add_server(msg, from);
	   break;

	 case DDI_ARR_ACC:
	   DDI_ARR_acc_server(msg, from);
	   break;
	   
         }

     # ifdef CRAY_MPI
       /* ----------------------------------------------------- *\
          Repost the asynchronus IRecv for remote data requests
       \* ----------------------------------------------------- */
          MPI_Irecv(&p[index_ds[i]],size,MPI_BYTE,from,37,
                     comm->world_comm,&r[index_ds[i]]);
       }
     # endif

     } while(server);


   /* -------------------------------------------------------------------------- *\
      If using MPI and not TCP socekts -- cancel/release the persistent receives
   \* -------------------------------------------------------------------------- */
    # if defined DDI_MPI && !defined DDI_SOC

      /* Working on this bit */

    # endif


   /* ------------------------------- *\
      Finalize communication and exit
   \* ------------------------------- */
      DDI_Finalize();
      exit(0);
   }
 
