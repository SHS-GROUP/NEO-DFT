/* ------------------------------------------------------ *\
   Wrapper for common unix system calls w/ error messages
   ======================================================
   Author: Ryan M. Olson
   10 Jun 09 - RMO - add usleep routine
   13 May 10 - SS  - change memory allocation to work on Windows
\* ------------------------------------------------------ */
 # include "ddi_base.h"

 # if defined WINDOWS64
     # include <stdlib.h>
     # include <malloc.h>
 # endif

/* ------------------------- *\
   System Wrapper for malloc
\* ------------------------- */
   void *Malloc(size_t size) {
      void *addr;
      if((addr = malloc(size)) == NULL) {
         fprintf(stderr," Error: malloc(%lu) failed.\n",
                        (unsigned DDI_INT64) size);
         Fatal_error(911);
      }
      return addr;
   }

/* -------------------------------- *\
   System Wrapper for gethostbyname
\* -------------------------------- */
 # if defined DDI_SOC
   struct hostent *Gethostbyname(const char *name) {
      struct hostent *hp;
      if((hp = gethostbyname(name)) == NULL) {
         switch(h_errno) {
           case TRY_AGAIN:      fprintf(stdout," Warning: Gethostbyname returned TRY_AGAIN.\n");
                                fflush(stdout);  return Gethostbyname(name); break;
           case NO_DATA:        fprintf(stdout," Error: Gethostbyname(%s) returned NO_DATA.\n",name); break;
           case NO_RECOVERY:    fprintf(stdout," Error: Gethostbyname(%s) returned NO_RECOVERY.\n",name); break;
           case HOST_NOT_FOUND: fprintf(stdout," Error: Gethostbyname(%s) returned HOST_NOT_FOUND.\n",name); break;
           default:             fprintf(stdout," Error: Gethostbyname(%s) returned an unknown error (%s).\n",name,DDI_Id()); break;
         } 
         Fatal_error(911);
      }
      return hp;
   }
 # endif


/* ------------------------ *\
   System Wrapper for chdir
\* ------------------------ */
   int Chdir(const char *path) {
      int ret;
      if((ret=chdir(path)) < 0) {
         switch(errno) {
            case EFAULT:  fprintf(stdout,"Error: chdir(%s) failed (errno=EFAULT).\n",path); break;
            case EACCES:  fprintf(stdout,"Error: chdir(%s) failed (errno=EACCES).\n",path); break;
            case ENOTDIR:  fprintf(stdout,"Error: chdir(%s) failed (errno=ENOTDIR).\n",path); break;
            case ENAMETOOLONG:  fprintf(stdout,"Error: chdir(%s) failed (errno=ENAMETOOLONG).\n",path); break;
            default: fprintf(stdout,"Error: chdir(%s) failed (errno=unknown).\n",path); break;
         }
         Fatal_error(911);
      } 
      return ret;
   }


/* ------------------------- *\
   System wrapper for execvp
\* ------------------------- */
   int Execvp(const char *file,char *const argv[]) {
      int ret;
      if((ret=execvp(file,argv)) == -1) {
         switch(errno) {
            case EIO:    fprintf(stdout,"Error: execvp(%s,args) failed (EIO).\n",file); break;
            case ELOOP:  fprintf(stdout,"Error: execvp(%s,args) failed (ELOOP).\n",file); break;
            case E2BIG:  fprintf(stdout,"Error: execvp(%s,args) failed (E2BIG).\n",file); break;
            case ENOMEM: fprintf(stdout,"Error: execvp(%s,args) failed (ENOMEM).\n",file); break;
            case EACCES: fprintf(stdout,"Error: execvp(%s,args) failed (EACCES).\n",file); break;
            case ENOENT: fprintf(stdout,"Error: execvp(%s,args) failed (ENOENT).\n",file); break;
            case EFAULT: fprintf(stdout,"Error: execvp(%s,args) failed (EFAULT).\n",file); break;
            case ENOTDIR: fprintf(stdout,"Error: execvp(%s,args) failed (ENOTDIR).\n",file); break;
            case ENOEXEC: fprintf(stdout,"Error: execvp(%s,args) failed (ENOEXEC).\n",file); break;
            case ETXTBSY: fprintf(stdout,"Error: execvp(%s,args) failed (ETXTBSY).\n",file); break;
            case ENAMETOOLONG: fprintf(stdout,"Error: execvp(%s,args) failed (ENAMETOOLONG).\n",file); break;
            default:      fprintf(stdout,"Error: execvp(%s,args) failed (errno=unknown).\n",file); break;
         }
      }
      return ret;
}


/*
 * Portable implementation of usleep(3C), microsecond timer.
 * Kevin Thomas, kjt@cray.com, 12 Sep 2001.
 */
int usleep(unsigned int micro) {
  struct timeval timeout;
  fd_set fdset;

  /* Must be less than one second. */
  if (micro >= 1000000) { errno = EINVAL; return -1; }

  timeout.tv_sec = 0;
  timeout.tv_usec = micro;
  FD_ZERO (&fdset);
  return select(1,&fdset,NULL,NULL,&timeout);
}
