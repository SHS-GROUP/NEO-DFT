/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Subroutines associated with signal handling.
 *
 * Author: Ryan M. Olson
 * 13 May 10 - SS  - accomodate Window's sleep differences
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"

/* -------------------------------------------------- *\
   Signal handling function for SIGURG.
   Recieve a message from the kickoff program, which:
   1) Resets the probe counter, or
   2) Kill the DDI process.
\* -------------------------------------------------- */
   void _sigurg_hndlr(int signo) {
   # if defined DDI_SOC
     char msg = 1;
     recv(gv(kickoffsock),&msg,1,MSG_OOB);
     if(msg != 1) Fatal_error(signo);
   # endif
   }

void Fatal_error(int signo) {
  DDI_Error(signo, NULL);
}

/** @see ddi_error.h */
void DDI_Error(int signo, char *message) {
  
  /* A fatal error has occurred and the program must quit.
     Trapping further signals might interfer with the cleanup */
  signal(SIGURG,SIG_IGN);
  signal(SIGPIPE,SIG_IGN);
  
  if (message == NULL) {
    switch(signo) {
    case SIGFPE:
      fprintf(DDI_STDERR, "%s: trapped a floating point error (SIGFPE).\n",DDI_Id());
      break;
    case SIGTERM:
      fprintf(DDI_STDERR, "%s: trapped a termination signal (SIGTERM).\n",DDI_Id());
      break;
    case SIGURG:
      fprintf(DDI_STDERR, "%s: terminated upon request.\n",DDI_Id());
      break;
    case SIGSEGV:
      fprintf(DDI_STDERR, "%s: trapped a segmentation fault (SIGSEGV).\n",DDI_Id());
      break;
    case SIGPIPE:
      fprintf(DDI_STDERR, "%s: SIGPIPE trapped.\n",DDI_Id());
      break;
    case SIGINT:
      fprintf(DDI_STDERR, "%s: SIGINT trapped.\n",DDI_Id());
      break;
    case SIGILL:
      fprintf(DDI_STDERR, "%s: SIGILL trapped.\n",DDI_Id());
      break;
    case SIGQUIT:
      fprintf(DDI_STDERR, "%s: SIGQUIT trapped.\n",DDI_Id());
      break;
    case SIGXCPU:
      fprintf(DDI_STDERR, "%s: process exceeded CPU time limit (SIGXCPU).\n",DDI_Id());
      break;
      
    /* DDI errors */
    case DDI_MAX_SMP_PROCS_ERROR:
      fprintf(DDI_STDERR, "%s: error code %i (MAX_SMP_PROCS exceeded).\n", DDI_Id(), signo);
      break;
    case DDI_MAX_NODES_ERROR:
      fprintf(DDI_STDERR, "%s: error code %i (MAX_NODES exceeded).\n", DDI_Id(), signo);
      break;
    case DDI_MAX_PROCESSORS_ERROR:
      fprintf(DDI_STDERR, "%s: error code %i (MAX_PROCESSORS exceeded).\n", DDI_Id(), signo);
      break;
      
    default:
      fprintf(DDI_STDERR, "%s: error code %i\n",DDI_Id(), signo);
      break;
    };
  }
  else {
    fprintf(DDI_STDERR, "%s: error code %i (%s)\n", DDI_Id(), signo, message);
  }
  
  DDI_Abort(signo);
}

/** @see ddi_error.h */
void DDI_Abort(int signo) {

# if defined DDI_MPI
  int i;
  char ack=0;
# endif
  
  int me,np;
  SMP_Data *smp_data = NULL;
  DDI_NProc(&np,&me);

  fflush(stdout);
  fflush(stderr);
  
   /* If running the mixed code, kill the other processes by SIGURG messages over
      TCP sockets, rather than dealing with the implementation specific MPI_ABORT */
#if defined DDI_MPI && defined DDI_SOC
  if(signo != SIGURG) { /* kill the compute processes */
    fprintf(stdout,"%s: Killing remaining DDI processes.\n",DDI_Id()); fflush(stdout);
    for(i=0; i<np; i++) {
      if(gv(sockets)[i] < 0 || i == me) continue;
      Send(gv(sockets)[i],&ack,1,MSG_OOB);
      STD_DEBUG((stdout,"%s: Sending SIGURG to %i\n",DDI_Id(),i));
    }
  }
  
  if(me < np && gv(sockets)[me+np] > 0) {
    Send(gv(sockets)[me+np],&ack,1,MSG_OOB);
    STD_DEBUG((stdout,"%s: Sending SIGURG to my data server %i.\n",DDI_Id(),me+np));
  }
  
 #if defined WINDOWS
   Sleep(1*1000);
 #else
   sleep(1);
 #endif

#endif
  
#if !(defined DDI_MPI2 || defined DDI_ARMCI)

  /* If used, delete the System V semaphore */
#if defined USE_SYSV
  DDI_Sem_remove(gv(dda_access));
  DDI_Sem_remove(gv(fence_access));
  if(gv(shmid) != 0) DDI_Shm_remove(gv(shmid));
#endif
  
#if FULL_SMP   
  DDI_Sem_remove(gv(dlb_access));
#endif

#endif /* !(defined DDI_MPI2 || defined DDI_ARMCI) */
  
  /* Clean up semaphores associated with shared-memory segments */
  while(smp_data = (SMP_Data *) SMP_find_end()) {
    fprintf(stdout,"%s: Cleaning up leftover SMP Array %i.\n",DDI_Id(),smp_data->handle);
    SMP_destroy(smp_data->handle);
  }
  
  /* If only using MPI, the cleanup is likely to be incomplete.  At the very least,
     the node on which the primary error occurred will be properly cleaned, however,
     some MPI implemenations send a hard kill (SIGKILL) to remote processes on
     MPI_Abort, which can not be trapped and those nodes will likely have semaphores
     remaining on them.  Note to administrators, if you have the ability to control
     your MPI implemenation(s), it is better to send a SIGTERM to the processes to
     allow proper cleanup, then a few seconds later, if they are not dead, send the
     very fatal SIGKILL (-9). */
  
#if defined DDI_ARMCI
  ARMCI_Cleanup();
#endif

# if defined DDI_MPI || !defined DDI_SOC
  MPI_Abort(MPI_COMM_WORLD,signo);
# endif
  
# if defined DDI_DUMP_CORE
  abort();
# else
  exit(signo);
# endif
}

