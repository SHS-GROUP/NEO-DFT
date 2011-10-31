/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Specialized LAPI Subroutines
 * 
 * Author: Ryan M. Olson
 * CVS $Id: ddi_lapi.c,v 1.1.1.1 2007/05/26 01:42:29 andrey Exp $
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"
 # if defined DDI_LAPI
 
/* ---------------------------------- *\
   LAPI Active Message Header Handler
\* ---------------------------------- */
   void *DDI_AM_Header_hndlr(lapi_handle_t *hndl, void *uhdr, uint *uhdrlen,
         uint *msglen, compl_hndlr_t **cmpl_hndlr, void **saved_info) {
   
      size_t nbytes = 0;
      size_t free_mem = 0;
      void *buffer = NULL;
      char *data_buffer = NULL;
      DDI_Patch *patch = (DDI_Patch *) uhdr;

      MAX_DEBUG((stdout,"%s: Entering LAPI Active Message Header Handler.\n",DDI_Id()))

   /* ----------------------------------------- *\
      If handling an accumulate, raise a fence.
   \* ----------------------------------------- */
    # if defined WINTEL
      if(patch->oper == DDI_ACC_OP) DDI_Fence_acquire(patch->handle);
    # else
      if(patch->oper == DDI_ACC) DDI_Fence_acquire(patch->handle);
    # endif


   /* ------------------------------------------ *\
      Determine the size of the buffer to create
   \* ------------------------------------------ */
      nbytes  = sizeof(DDI_Patch);
      nbytes += (nbytes % sizeof(double));
      nbytes += patch->size;
      patch->ds_buffer_size = patch->size;


   /* --------------------------- *\
      Malloc the temporary buffer
   \* --------------------------- */
      DDI_Memory_heap_malloc(&buffer,nbytes);
      
      
   /* ------------------------------------------------- *\
      Save the DDI_Patch info at the head of the buffer
   \* ------------------------------------------------- */
      memcpy(buffer,patch,sizeof(DDI_Patch));
      

   /* -------------------------------------- *\
      Information for the completion handler
   \* -------------------------------------- */      
      *cmpl_hndlr = (compl_hndlr_t *) &DDI_AM_Compl_hndlr;
      *saved_info = (void *) buffer;


   /* -------------------------------------- *\
      Determine the start of the data buffer
   \* -------------------------------------- */
      nbytes      -= patch->size;
      data_buffer  = (char *) buffer;
      data_buffer += nbytes;
     
      MAX_DEBUG((stdout,"%s: Leaving AM Header Handler.\n",DDI_Id()))
      return (void *) data_buffer;
   }
   
   
 
/* -------------------------------------- *\
   LAPI Active Message Completion Handler
\* -------------------------------------- */
   void DDI_AM_Compl_hndlr(lapi_handle_t *handl,void *param) {
   
      size_t nbytes;
      char *stack = (char *) param;
      DDI_Patch *patch = (DDI_Patch *) param;
      double *buffer = NULL;

      MAX_DEBUG((stdout,"%s: Entering LAPI Active Message Completion Handler.\n",DDI_Id()))

      nbytes  = sizeof(DDI_Patch);
      nbytes += (nbytes % sizeof(double));
      
      stack += nbytes;
      buffer = (double *) stack;

      nbytes += patch->ds_buffer_size;
            
      switch(patch->oper) {
       # if defined WINTEL
         case DDI_GET_OP:  DDI_Get_lapi_server(patch,buffer);  break;
         case DDI_PUT_OP:  DDI_Put_lapi_server(patch,buffer);  break;
         case DDI_ACC_OP:  DDI_Acc_lapi_server(patch,buffer);  break;
       # else
         case DDI_GET:  DDI_Get_lapi_server(patch,buffer);  break;
         case DDI_PUT:  DDI_Put_lapi_server(patch,buffer);  break;
         case DDI_ACC:  DDI_Acc_lapi_server(patch,buffer);  break;
       # endif
         default:
            fprintf(stdout,"%s: unknown operation in data server handler.\n",DDI_Id());
            Fatal_error(911);
      }
      
      # if defined WINTEL
      if(patch->oper == DDI_ACC_OP) DDI_Fence_release(patch->handle);
      # else
      if(patch->oper == DDI_ACC) DDI_Fence_release(patch->handle);
      # endif
      DDI_Memory_heap_free(param,nbytes);

      MAX_DEBUG((stdout,"%s: Leaving LAPI AM Completion Handler.\n",DDI_Id()))
   }
   
 # else
   
   void DDI_LAPI_dummy(void) { return; }
   
 # endif
