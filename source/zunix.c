/* 12 June 2012

   Interface from FORTRAN/GAMESS to C/system calls
   ===============================================
 
   The amount of quiche is kept to a minimum, as there are only 4 routines!
          The first three routines amount to a FORTRAN wrapper around
          the C calls malloc and free, by Steve Elbert in January 1990.
   memory allocation: function invocation only: LOCMEM = MEMGET(NWORDS)
          NWORDS is the no. of 8-byte words of memory requested.
          The return value is the address of the allocated memory, if
          this is zero, the system doesn't have the memory available.
   memory release: subroutine call only: CALL MEMREL(LOCMEM)
          Deallocate all memory above the address LOCMEM.
          (deallocated memory is not returned to system)
   memory address: function invocation only: LOCX = LADDRS(X)
          Returns the absolute address of X, in bytes.  This is needed
          only if the FORTRAN library does not include the LOC function.
   process sleeping: subroutine call only: CALL NAPTIME(seconds)
          which is used when MPI is being used to pass messages, and
          allows DDI data servers to sleep when they are not in use.

   Please note: 
   There really shouldn't be anything here for a machine other than these 4,
   so one should refrain from adding anything else to this file at all.
   Other routines are present only because of inadequate system libraries.
*/

/*--------------- Digital/Compaq/HP "alpha" --------------*/

#ifdef AXP64
#include <stdlib.h>

long memget_(nwords) long *nwords;
   { size_t nbytes;
     nbytes = (size_t) (*nwords+2)*8;
     return (long) malloc(nbytes); }

void memrel_(locmem) long *locmem;
   { free((void*) *locmem); }

#include <unistd.h>
void naptime_(nptime) int *nptime;
   { unsigned int istat;
     istat=sleep(*nptime); }

#endif
 

/*-------------------- Cray --------------------------*/
/*      the Cray XT uses the 64 bit Linux clause,     */
/*  which makes everything here pretty much obsolete  */

/*  --  Cray PVP or T3E --   */
#ifdef CRAY

int MEMGET(nwords) int *nwords;
   { int nbytes;
     nbytes = (*nwords+2)*8;
     return malloc(nbytes); }

void MEMREL(locmem) int *locmem;
   { free(*locmem); }

#include <unistd.h>
void NAPTIME(nptime) int *nptime;
   { unsigned int istat;
     istat=sleep(*nptime); }

#endif

/*  --  Cray X1 --   */
#ifdef CRAYX1

#include <unistd.h>

long memget_(nwords) long *nwords;
   { size_t nbytes;
     nbytes = (size_t) (*nwords+2)*8;
     return (long) malloc(nbytes); }

void memrel_(locmem) long *locmem;
   { free(*locmem); }

void naptime_(nptime) int *nptime;
   { unsigned int istat;
     istat=sleep(*nptime); }

#endif

/*  --  Cray XD1 --   */
#ifdef CRAYXD1
   
#include <stdlib.h>
#include <malloc.h>
#define FORTINT long

FORTINT memget_(nwords) FORTINT *nwords;
   { size_t nbytes;
     nbytes = (*nwords+2)*8;
     return (FORTINT) malloc(nbytes); }

void memrel_(locmem) FORTINT *locmem;
   { free((void*)*locmem); }

#include <unistd.h>
void naptime_(nptime) FORTINT *nptime;
   { sleep((unsigned int) *nptime); }

#endif


/*  --  Cray XT3 --   */
#ifdef CRAYXT3

#include <stdlib.h>
#include <mpp/shmem.h>
#define FORTINT long

FORTINT memget_(nwords) FORTINT *nwords;
   { size_t nbytes;
     nbytes = (*nwords+2)*8;
     return (FORTINT) shmalloc(nbytes); }

void memrel_(locmem) FORTINT *locmem;
   { shfree((void*)*locmem); }

#include <unistd.h>
void naptime_(nptime) FORTINT *nptime;
   { sleep((unsigned int) *nptime); }

#endif

/*    we used this with MPI and data servers, above is for SHMEM-type DDI
#include <stdlib.h>
#include <malloc.h>
#define FORTINT long

FORTINT memget_(nwords) FORTINT *nwords;
   { size_t nbytes;
     nbytes = (*nwords+2)*8;
     return (FORTINT) malloc(nbytes); }

void memrel_(locmem) FORTINT *locmem;
   { free((void*)*locmem); }

#include <unistd.h>
void naptime_(nptime) FORTINT *nptime;
   { sleep((unsigned int) *nptime); }

#endif
*/


/*-------------- Fujitsu VPP and AP ---------------*/

#ifdef FUJITSU

int memget_(nwords) int *nwords;
   { int nbytes;
     nbytes = (*nwords+2)*8;
     return malloc(nbytes); }

void memrel_(locmem) int *locmem;
   { free(*locmem); }

int laddrs_(arg) int arg;
   { return(arg); }

/*   This is the only machine type missing a NAPTIME implementation   */
/*   If you are trying to link GAMESS, please send a copy of your     */
/*   version of NAPTIME (many examples here!) to Mike Schmidt         */

#endif


/*--------------------- Hewlett-Packard --------------------------*/
/*    this version uses FORTRAN callable library routines only    */
/*                 no need to compile this file                   */


/*----------------- International Business Machines ----------------*/
/*  IBM32    = AIX operating system, 32 bits                        */
/*  IBM64    = AIX operating system, 64 bits, includes SP machines  */
/*  IBMBG    = Linux operating system, 32 bits, used on Blue Gene.  */
/*  IBMPPC64 = Linux operating system, 64 bits, e.g. OpenPower      */
 
#ifdef IBM32

int memget(nwords) int *nwords;
   { int nbytes;
     nbytes = (*nwords+2)*8;
     return malloc(nbytes); }

void memrel(locmem) int *locmem;
   { free(*locmem); }

int laddrs(arg) int arg;
   { return(arg); }

#include <unistd.h>
void naptime(nptime) useconds_t *nptime;
   { int istat;
     istat=sleep(*nptime); }

#endif

#ifdef IBM64

#include <stdlib.h>
long memget(long *nwords)
  { size_t nbytes;
    nbytes = (*nwords + 2L) * 8L;
    return (long) malloc (nbytes); }

void memrel(long **locmem)
  { free (*locmem); }

long laddrs(void *arg)
  { return (long) arg; }

#include <unistd.h>
void naptime(long *nptime)
  { int istat;
    useconds_t delay;
    delay = (useconds_t) *nptime;
    istat=sleep(delay); }

#endif

/* it would be most excellent to get rid of this one, for LINUX32 */
#ifdef IBMBG

#include <unistd.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <sys/times.h>

int memget(int *nwords)
{
   int nbytes;
   nbytes = (*nwords+2)*8;
  return (int)malloc(nbytes);
}

void memrel(int *locmem)
{
   free((void *)*locmem);
}

int laddrs(arg) int arg;
   { return(arg); }

void naptime(nptime) int *nptime;
   { unsigned int istat;
     istat=sleep(*nptime); }

/*
     etime - return elapsed execution time - see comments elsewhere.
*/

double etime(float *a)
   { double        elapsed;
     struct tms    buf;
     elapsed= (float) times(&buf)/ (float) CLK_TCK;
     a[0]=(float)(buf.tms_utime + buf.tms_cutime)/CLK_TCK;
     a[1]=(float)(buf.tms_stime + buf.tms_cstime)/CLK_TCK;
     return(elapsed); }

/*
     fdate - return date and time in an ASCII string
     subroutine fdate(string)
     character*24 string

     returns the current date and time as a 24 character string in the
     format of ctime(3).  Neither `newline' nor NULL will be included.
     by Istvan Cserny, June 21, 1995
*/

void fdate(char *cht, int cht_len)
   { static time_t t;
     time(&t);
     strncpy(cht,ctime(&t),24L);
     return; }

#endif

/* IBMPPC64 differs from LINUX64 by including LADDRS,   */
/* and by having no underscores in the function names.  */
#ifdef IBMPPC64

#include <stdlib.h>
#include <malloc.h>

#define FORTINT long

FORTINT memget(nwords) FORTINT *nwords;
   { size_t nbytes;
     nbytes = (*nwords+2)*8;
     return (FORTINT) malloc(nbytes); }

void memrel(locmem) FORTINT *locmem;
   { free((void*)*locmem); }

long laddrs(void *arg)
  { return (long) arg; }

#include <unistd.h>
void naptime(nptime) FORTINT *nptime;
   { sleep((unsigned int) *nptime); }

/*  
     The 64 bit version of the etime clone was revised March 2007, after
     looking at "man 2 times".  This was motivated by Suse Enterprise 10,
     which is using a clock tick of 250, not 100!  Probably most Linux
     systems running at 64 bit were installed since 'times' was updated.

     The 32 bit routine was left untouched, for it gives us no trouble.
     Some of these systems may be quite a bit older, and the code has
     worked on everything we tried for many years, so let the dog sleep.

     See the 32 bit version of ETIME for more information about using
     it from FORTRAN.  Args here are identical to the 32 bit code.

     The structure 'tms' is set up in <sys/time.h>, as
        struct tms {clock_t tms_utime;   namely, user time
                    clock_t tms_stime;   namely, system time
                    clock_t tms_cutime;  namely, user time of children
                    clock_t tms_cstime;  namely, system time of children
     The return value of the function times, which also sets the structure,
     is the number of clock ticks since an arbitrary time in the past:
        # include <sys/times.h>
        clock_t times(struct tms *buf);
*/
#include <sys/times.h>
#include <unistd.h>

double etime(float *a)
   { double        elapsed;
     clock_t       elapticks;
     struct tms    buf;
     long          sysconf(int name);
     long          POSIX_CLK_TCK;

     elapticks = times(&buf);
     POSIX_CLK_TCK = sysconf(_SC_CLK_TCK);

     elapsed= (double) elapticks/ (double) POSIX_CLK_TCK;
     a[0] = (float) (buf.tms_utime + buf.tms_cutime) / (float) POSIX_CLK_TCK;
     a[1] = (float) (buf.tms_stime + buf.tms_cstime) / (float) POSIX_CLK_TCK;
     return(elapsed);
   }
#endif


/*----- 32 bit PC (e.g. Pentium, Athlon, ...) running Linux -------*/
/*--------- this also includes the Apple running MAC OS X ---------*/
/*   ETIME is included below, used to replace the one that might   */
/*   be found in Linux libraries, with a different return value.   */

#ifdef LINUX32

#include <stdlib.h>
int memget_(nwords) int *nwords;
   { int nbytes;
     nbytes = (*nwords+2)*8;
     return (int) malloc(nbytes); }

void memrel_(locmem) int *locmem;
   { free((void*)*locmem); }

int laddrs_(arg) int arg;
   { return(arg); }

#include <unistd.h>
void naptime_(nptime) int *nptime;
   { unsigned int istat;
     istat=sleep(*nptime); }

/*
     etime - return elapsed execution time - usage is:

        double precision function etime (tarray)
        real*4 tarray(2)

     the REAL*4 argument returns user and system CPU times, as is traditional
     for ETIME in vendor Unix.  The REAL*8 return value, however, is not the 
     sum of these, but rather it gives the elapsed wall clock time (thus it
     replaces the traditional TIME() function in vendor Unix libraries).

     written by Klaus-Peter Gulden, January 3, 1994

     note: CLK_TCK counts hundredths of seconds, usually, and is
           normally set by the system, through time.h's includes.
           gcc 4.x moved the definition, or maybe changed the name
           of the value?  If we don't find it, set it to 100, and
           hope that's really the correct value!
*/

#include <time.h>
#include <sys/times.h>

# ifndef CLK_TCK
# define CLK_TCK 100
# endif

double etime_(float *a)
   { double        elapsed;
     struct tms    buf;
     elapsed= (float) times(&buf)/CLK_TCK;
     a[0]= (float) (buf.tms_utime + buf.tms_cutime)/CLK_TCK;
     a[1]= (float) (buf.tms_stime + buf.tms_cstime)/CLK_TCK;
     return(elapsed); }

/*
     fdate - return date and time in an ASCII string
     subroutine fdate(string)
     character*24 string

     returns the current date and time as a 24 character string in the
     format of ctime(3).  Neither `newline' nor NULL will be included.
     by Istvan Cserny, June 21, 1995
*/

#include <strings.h>
void fdate_(char *cht, int cht_len)
   { static time_t t;
     time(&t);
     strncpy(cht,ctime(&t),24L);
     return; }

#endif

/*----------------- 64 bit systems running Linux -----------------*/
/*    This is tested for Itanium2, and Opteron based machines.    */
/*    This also used by Cray XT systems' Compute Node Linux.      */

/*    For more information, see 32 bit equivalent just above.     */
/*    Note that so far, -laddrs- has not been copied here.        */
/*    Note that so far, -fdate- has not been copied here.         */

#ifdef LINUX64

#include <stdlib.h>
#include <malloc.h>

#define FORTINT long

FORTINT memget_(nwords) FORTINT *nwords;
   { size_t nbytes;
     nbytes = (*nwords+2)*8;
     return (FORTINT) malloc(nbytes); }

void memrel_(locmem) FORTINT *locmem;
   { free((void*)*locmem); }

#include <unistd.h>
void naptime_(nptime) FORTINT *nptime;
   { sleep((unsigned int) *nptime); }

/*  
     The 64 bit version of the etime clone was revised March 2007, after
     looking at "man 2 times".  This was motivated by Suse Enterprise 10,
     which is using a clock tick of 250, not 100!  Probably most Linux
     systems running at 64 bit were installed since 'times' was updated.

     The 32 bit routine was left untouched, for it gives us no trouble.
     Some of these systems may be quite a bit older, and the code has
     worked on everything we tried for many years, so let the dog sleep.

     See the 32 bit version of ETIME for more information about using
     it from FORTRAN.  Args here are identical to the 32 bit code.

     The structure 'tms' is set up in <sys/time.h>, as
        struct tms {clock_t tms_utime;   namely, user time
                    clock_t tms_stime;   namely, system time
                    clock_t tms_cutime;  namely, user time of children
                    clock_t tms_cstime;  namely, system time of children
     The return value of the function times, which also sets the structure,
     is the number of clock ticks since an arbitrary time in the past:
        # include <sys/times.h>
        clock_t times(struct tms *buf);
*/
#include <sys/times.h>
#include <unistd.h>

double etime_(float *a)
   { double        elapsed;
     clock_t       elapticks;
     struct tms    buf;
     long          sysconf(int name);
     long          POSIX_CLK_TCK;

     elapticks = times(&buf);
     POSIX_CLK_TCK = sysconf(_SC_CLK_TCK);

     elapsed= (double) elapticks/ (double) POSIX_CLK_TCK;
     a[0] = (float) (buf.tms_utime + buf.tms_cutime) / (float) POSIX_CLK_TCK;
     a[1] = (float) (buf.tms_stime + buf.tms_cstime) / (float) POSIX_CLK_TCK;
     return(elapsed);
   }

#endif
 

/*----------------- NEC SX series -----------------*/

#ifdef NECSX
#include <unistd.h>

void memget_(nwords,sxaddr)
 long long *nwords, **sxaddr ;
   {char *malloc(),*card;
    unsigned nbytes;
    nbytes = (unsigned) ((*nwords+2)*8);
    card =  malloc (nbytes);
    *sxaddr = (long long *) card; }

void memrel_(locmem)
 long long *locmem;
 {void free();
  char* ptrmem;
  ptrmem = (char *)locmem;
  free(ptrmem); }

void laddrs_(arg,addrarg)
long long *arg, **addrarg;
   { *addrarg = arg ; }

void naptime_(long *nptime)
  { long istat;
    istat=(long) sleep((unsigned)*nptime); }

#endif


/*-------------- Silicon Graphics, Incorporated --------------*/
/*                 ancient MIPS systems, only                 */

#ifdef SGI32
#define FORTINT int
#endif
#ifdef SGI64
#define FORTINT long
#endif

#if (defined SGI32) || (defined SGI64)

int memget_(nwords) FORTINT *nwords;
   { FORTINT nbytes;
     nbytes = (*nwords+2)*8;
     return malloc(nbytes); }

void memrel_(locmem) FORTINT *locmem;
   { free(*locmem); }

#include <unistd.h>
void naptime_(FORTINT *nptime)
   { sginap((long) *nptime); }

#endif


/*----------------- SUN -----------------*/

#ifdef SUN32
#define FORTINT int
#endif
#ifdef SUN64
#define FORTINT long
#endif

#if (defined SUN32) || (defined SUN64)

FORTINT memget_(nwords) FORTINT *nwords;
   { FORTINT nbytes; FORTINT malloc();
     nbytes = (*nwords+2)*8;
     return malloc(nbytes); }

void memrel_(locmem) FORTINT *locmem;
   { void free();
     free(*locmem); }

FORTINT laddrs_(arg) FORTINT arg;
   { return(arg); }

#include <unistd.h>
void naptime_(nptime) FORTINT *nptime;
   { unsigned int istat,delay;
     delay = (unsigned int) *nptime;
     istat = sleep(delay); }

#endif


/*-------------- Windows --------------*/

#ifdef WINDOWS32

#include <stdlib.h>

int memget_(nwords) int *nwords;
   { int nbytes;
     nbytes = (*nwords+2)*8;
     return (int) malloc(nbytes); }

void memrel_(locmem) int *locmem;
   { free((void*) *locmem); }

int laddrs_(arg) int arg;
   { return(arg); }

#include <Windows.h>
void naptime_(nptime) int *nptime;
   { Sleep((DWORD) *nptime*1000); }

#include <mpi.h>
double windowstimer_(double *a) 
{
  /* 
     Visit this webpage for more information on what is being done
     http://cplus.about.com/od/howtodothingsin1/a/timing.htm
  */
  
  LARGE_INTEGER clockticks;
  LARGE_INTEGER frequency;
  QueryPerformanceFrequency( &frequency );
  QueryPerformanceCounter( &clockticks );
  *a = MPI_Wtime();
  return ((double)clockticks.QuadPart/(double)frequency.QuadPart);
}

#endif

#ifdef WINDOWS64

#include <stdlib.h>
 #ifdef WINTEL
  #include <malloc.h>
 #endif
#define FORTINT __int64

FORTINT memget_(nwords) FORTINT *nwords;
   { size_t nbytes;
     nbytes = (*nwords+2)*8;
     return (FORTINT) malloc(nbytes); }

void memrel_(locmem) FORTINT *locmem;
   { free((void*) *locmem); }

#include <Windows.h>
void naptime_(nptime) FORTINT *nptime;
   { Sleep((DWORD) *nptime*1000); }

#ifndef WINDOWSAZURE

#include <mpi.h>
double windowstimer_(double *a) 
{
  /* 
     Visit this webpage for more information on what is being done
     http://cplus.about.com/od/howtodothingsin1/a/timing.htm
  */
  
  LARGE_INTEGER clockticks;
  LARGE_INTEGER frequency;
  QueryPerformanceFrequency( &frequency );
  QueryPerformanceCounter( &clockticks );
  *a = MPI_Wtime();
  return ((double)clockticks.QuadPart/(double)frequency.QuadPart);
}
double walltimer_()
{
  return (MPI_Wtime());
}
#endif

#endif
