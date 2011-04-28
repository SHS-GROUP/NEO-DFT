/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Sync/Barrier routines
 *
 * Author: Ryan M. Olson
 * 10 Jun 09 - RMO - allow for one data server/node on Cray
\* -------------------------------------------------------------------- */ 
 # include "ddi_base.h"

/* ----------------------------------------------------------------- *\
\* ----------------------------------------------------------------- */
   void DDI_Sync(int tag) {
      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);

      DEBUG_ROOT(LVL1,(stdout," DDI: Sync point (%i).\n",tag))
      DEBUG_OUT(LVL3,(stdout,"%s: entering DDI_Sync (tag=%i).\n",DDI_Id(),tag))

      Comm_sync(tag,comm);

      DEBUG_OUT(LVL3,(stdout,"%s: leaving DDI_Sync (tag=%i).\n",DDI_Id(),tag))
      DEBUG_ROOT(LVL2,(stdout," DDI: leaving DDI_Sync (tag=%i).\n",tag))
   }

   
   void DDI_Sync_smp(int tag) {
      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);

      DEBUG_ROOT(LVL1,(stdout," DDI: Sync point (%i).\n",tag))
      DEBUG_OUT(LVL3,(stdout,"%s: entering DDI_Sync_smp (tag=%i).\n",DDI_Id(),tag))

      Comm_sync_smp(comm);

      DEBUG_OUT(LVL3,(stdout,"%s: leaving DDI_Sync_smp (tag=%i).\n",DDI_Id(),tag))
      DEBUG_ROOT(LVL2,(stdout," DDI: leaving DDI_Sync_smp (tag=%i).\n",tag))
   }


   void DDI_Sync_node(int tag) {

      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);

      DEBUG_ROOT(LVL1,(stdout," DDI: Sync point (%i).\n",tag))
      DEBUG_OUT(LVL3,(stdout,"%s: entering DDI_Sync_node (tag=%i).\n",DDI_Id(),tag))

      Comm_sync_node(tag,comm);

      DEBUG_OUT(LVL3,(stdout,"%s: leaving DDI_Sync_node (tag=%i).\n",DDI_Id(),tag))
      DEBUG_ROOT(LVL2,(stdout," DDI: leaving DDI_Sync_node (tag=%i).\n",tag))
   }


   void DDI_Sync_comm(int tag,int comm_id) {

      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(comm_id);

      DEBUG_ROOT(LVL1,(stdout," DDI: Sync point (%i).\n",tag))
      DEBUG_OUT(LVL3,(stdout,"%s: entering DDI_Sync_comm (tag=%i).\n",DDI_Id(),tag))

      Comm_sync(tag,comm);

      DEBUG_OUT(LVL3,(stdout,"%s: leaving Sync_comm (tag=%i).\n",DDI_Id(),tag))
      DEBUG_ROOT(LVL2,(stdout," DDI: leaving Sync_comm (tag=%i).\n",tag))
   }

   
   void Comm_sync(int tag,const DDI_Comm *comm) {
      int sum_tag = tag;

      DEBUG_OUT(LVL4,(stdout,"%s: entering Comm_sync (tag=%i).\n",DDI_Id(),tag))

    # if defined DDI_ARMCI
      DDI_ARMCI_Barrier(comm->compute_comm);
    # elif defined DDI_MPI && !defined USE_DDI_COLLECTIVE_ROUTINES
      MPI_Barrier(comm->compute_comm);
    # else
   /* ---------------------------- *\
    * Single processes return now.
   \* ---------------------------- */
      if(comm->np == 1) return;

   /* ---------------------------- *\
    * Global Sum = Synchronization
   \* ---------------------------- */
      Comm_gsum(&sum_tag,1,sizeof(int),comm);

   /* -------------------- *\
    * Ensure matching tags
   \* -------------------- */
      sum_tag /= comm->np;
      if(sum_tag != tag) {
         fprintf(stdout," ERROR DDI_Sync: tag %i mismatched on rank %i (%i).\n",
                 tag,comm->me,sum_tag);
         Fatal_error(911);
      }
    # endif

      DEBUG_OUT(LVL4,(stdout,"%s: leaving Comm_sync (tag=%i).\n",DDI_Id(),tag))
   }


/* ----------------------------------------------------------------- *\
\* ----------------------------------------------------------------- */
   void Comm_sync_list(const DDI_List *list) {
      int buffer;
      int data = list->tag;

      DEBUG_OUT(LVL6,(stdout,"%s: entering Comm_sync_list.\n",DDI_Id()))

#if defined DDI_ARMCI
      DDI_ARMCI_Barrier(list->comm);
#elif defined DDI_MPI && !defined USE_DDI_COLLECTIVE_ROUTINES
      MPI_Barrier(list->comm);
#else
      Comm_gsum_list(&data,1,&buffer,1,DDI_INT,list);
      Comm_bcast_list(&data,sizeof(int),list);
      data /= list->np;
      if(data != list->tag) {
         fprintf(stdout," WARNING Comm_sync_list: mismatched message tag %i \
                          (on list rank %i)\n",list->tag,list->me);
         Fatal_error(911);
      }
#endif
      DEBUG_OUT(LVL6,(stdout,"%s: leaving Comm_sync_list.\n",DDI_Id()))
   }

   void Comm_sync_node(int tag,const DDI_Comm *comm) {
      DDI_List list;
      DEBUG_OUT(LVL6,(stdout,"%s: entering Comm_sync_node.\n",DDI_Id()))
    # if defined DDI_MPI && USE_MPI_BARRIER
      MPI_Barrier(comm->node_comm);
    # else
      list.np   = comm->nn;
      list.me   = comm->my;
      list.tag  = tag;
      list.root = 0;
      list.pids = comm->node_master;
    # if defined DDI_MPI
      list.comm = comm->node_comm;
    # endif
      Comm_sync_list(&list);
    # endif
      DEBUG_OUT(LVL6,(stdout,"%s: leaving Comm_sync_node.\n",DDI_Id()))
   }

   void Comm_sync_smp(const DDI_Comm *comm) {
      DDI_List list;
      DEBUG_OUT(LVL6,(stdout,"%s: entering Comm_sync_smp.\n",DDI_Id()))
    # if defined DDI_MPI && USE_MPI_BARRIER
      MPI_Barrier(comm->smp_comm);
    # else
      list.np   = comm->np_local;
      list.me   = comm->me_local;
      list.tag  = 345;
      list.root = 0;
      list.pids = comm->smp_pid;
    # if defined DDI_MPI
      list.comm = comm->smp_comm;
    # endif
      Comm_sync_list(&list);
    # endif
      DEBUG_OUT(LVL6,(stdout,"%s: leaving Comm_sync_smp.\n",DDI_Id()))
   }


/*
       a routine that is used only by code that is commented out
       but it seems reasonable to keep it here anyway...

   void Comm_sync_smp_fast(const DDI_Comm *comm,DDI_Data *data) {
      static int sync_num = 0;
      int i;
      CS_int *sync_array = (CS_int *) data->sync_array;

      sync_num = (sync_num+1) % 10;
      sync_array[comm->me_local].val = sync_num;
      for(i=0; i<comm->np_local; i++) {

         while(sync_array[i].val != sync_num && sync_array[i].val != (sync_num+1)%10) {
            sched_yield();
         }
      }
   }

 */
