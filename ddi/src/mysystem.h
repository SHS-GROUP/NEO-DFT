/* ----------------------------------------------------------- *\
   Include file containing common include files and prototypes
   for wrappers to system calls.

   Author: Ryan M. Olson
   10 Jun 09 - RMO - microsleep timer added (needs param.h)
   13 May 10 - SS  - porting to Windows
\* ----------------------------------------------------------- */

/* --------------------------------------- *\
   System Dependent definitions: Sun SPARC
\* --------------------------------------- */
 # if (defined SUN32 || defined SUN64)
 # define _REENTRANT
 # endif

 # define NO_COLLECTIVE_SMP

/* --------------------------- *\
   Commonly used include files
\* --------------------------- */
 # include <stdio.h>
 # include <fcntl.h>
 # include <errno.h>
 # include <stdlib.h>
 # include <signal.h>
 # include <string.h>
 
 # ifndef WINDOWS
   # include <unistd.h>
   # include <strings.h>
 # endif

/* ---------------------------------------------- *\
   System includes: come before standard includes
\* ---------------------------------------------- */
   # include <sys/types.h>

 # ifndef WINDOWS
   # include <sys/param.h>
   # include <sys/time.h>
   # include <sys/resource.h>
 # endif 

/* ------------------------------------------ *\
   System V IPC Includes & Wrapper functions
   Subroutine definitions found in sysv_ipc.c
\* ------------------------------------------ */
 # if defined USE_SYSV

 # if defined __CYGWIN__ 
 #   define _KERNEL
 # endif

 # include <sys/ipc.h>
 # include <sys/shm.h>
 # include <sys/sem.h>

   int   Shmget(key_t,size_t,int);
   int   Shmctl(int,int,struct shmid_ds*);
   void *Shmat(int,void*,int);

   int   Semget(key_t,int,int);
   int   Semop(int,struct sembuf*,size_t);

 # if defined __CYGWIN__ 
 #   if !defined SHM_R 
 #     define SHM_R IPC_R
 #   endif
 #   if !defined SHM_W 
 #     define SHM_W IPC_W
 #   endif
 # endif

 # endif


/* -------- *\
   pThreads
\* -------- */
 # if (!defined CRAYXT3 && !defined IBMBG && !defined WINDOWS)
 # include <pthread.h>
 # endif


/* --------------------------------------------- *\
   TCP Sockets Includes & Wrapper functions
   Subroutine definitions found in tcp_sockets.c
\* --------------------------------------------- */
 # if defined DDI_SOC
 # include <netdb.h>
 # include <sys/socket.h>
 # include <netinet/in.h>
 # include <netinet/tcp.h>

   int     Accept(int,struct sockaddr*,socklen_t*);
   int     Connect(int,const struct sockaddr*,socklen_t);
   ssize_t Recv(int,void*,size_t,int);  /* Fully blocking receive -- MSG_WAITALL */
   ssize_t Send(int,const void*,size_t,int);
 # endif


/* ---------------- *\
   Microsleep timer
\* ---------------- */
/* int usleep(useconds_t); */


/* --- *\
   MPI
\* --- */
 # if defined DDI_MPI
 # include <mpi.h>
 # endif


/* ---------------- *\
   LAPI -- IBM SP's
\* ---------------- */
 # if defined DDI_LAPI
 # include <lapi.h>
 # endif


/* --------------------------------------------- *\
   Wrapper functions for standard system calls
   Subroutines definitions found in std_system.c
\* --------------------------------------------- */
   void *Malloc(size_t);
   int Chdir(const char *);
   int Execvp(const char *,char *const argv[]);

 # if defined DDI_SOC
   struct hostent *Gethostbyname(const char*);
 # endif

/* ----------------------- *\
   Define Max & Min Macros
\* ----------------------- */
 # define max(a,b) ((a) > (b) ? (a) : (b))
 # define min(a,b) ((a) < (b) ? (a) : (b))

/*        section creating Windows support      */
# if defined WINDOWS
 
 # include <signal.h>
 # include <limits.h> 
 
 # define SIGHUP          1
 /* 2 is used for SIGINT   on Windows */
 # define SIGQUIT         3
 /* 4 is used for SIGILL   on Windows */
 # define SIGTRAP         5
 # define SIGIOT          6
 # define SIGBUS          7
 /* 8 is used for SIGFPE   on Windows */
 # define SIGKILL         9
 # define SIGUSR1        10
 /* 11 is used for SIGSEGV on Windows */
 # define SIGUSR2        12
 # define SIGPIPE        13
 # define SIGALRM        14
 /* 15 is used for SIGTERM on Windows */
 # define SIGSTKFLT      16
 # define SIGCHLD        17 
 # define SIGCONT        18
 # define SIGSTOP        19
 # define SIGTSTP        20
 /* 21 is used for SIGBREAK on Windows */
 /* 22 is used for SIGABRT  on Windows */
 # define SIGTTIN        23
 # define SIGTTOU        24
 # define SIGURG         25
 # define SIGXCPU        26
 # define SIGXFSZ        27
 # define SIGVTALRM      28
 # define SIGPROF        29
 # define SIGWINCH       30
 # define SIGIO          31

 # define EIO             5
 # define ELOOP          31
 # define E2BIG           7
 # define ENOMEM         12
 # define EACCES         13
 # define ENOENT          2
 # define EFAULT         14
 # define ENOTDIR        20
 # define ENOEXEC         8
 # define ETXTBSY        16
 # define ENAMETOOLONG   38

 # if defined(_MSC_VER) || defined(_MSC_EXTENSIONS)
 # define DELTA_EPOCH_IN_MICROSECS  11644473600000000Ui64
 # else
 # define DELTA_EPOCH_IN_MICROSECS  11644473600000000ULL
 # endif

struct timezone 
{
  int  tz_minuteswest; /* minutes W of Greenwich */
  int  tz_dsttime;     /* type of dst correction */
};

 # define RUSAGE_SELF       0
 # define RUSAGE_CHILDREN (-1)
 
 # include <Winsock2.h>

struct rusage
{
  struct timeval ru_utime;        /* user time used */
  struct timeval ru_stime;        /* system time used */
};

extern int gettimeofday(struct timeval * tv, struct timezone * tz);
extern int getrusage(int who, struct rusage * rusage);

# endif
/*      end of a long Windows-specific passage   */
