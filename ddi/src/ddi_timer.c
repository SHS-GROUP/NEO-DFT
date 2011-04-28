/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Subroutines associated with synchronizing the compute processes.
 *
 * Author: Ryan M. Olson
 * 13 May 10 - SS  - add timing passage for Windows
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"


/* -------------------------------------------------------------------- *\
   DDI_Timer_reset()
   =================
   Reset the cpu timer within the current operating scope.
\* -------------------------------------------------------------------- */
   void DDI_Timer_reset() {
      getrusage(RUSAGE_SELF,&gv(cpu_timer));
      gettimeofday(&gv(wall_timer),NULL);
   }


/* -------------------------------------------------------------------- *\
   DDI_Timer_output()
   ==================
   Synchronous barrier on compute processes, but also collects total
   cpu time from each compute process and prints the totals stdout.
\* -------------------------------------------------------------------- */
   void DDI_Timer_output() {

      int i,me,np;
      struct rusage mycputime;
      struct rusage *timings = NULL;
      struct timeval cpu_total;
      struct timeval wall_total;

      DDI_NProc(&np,&me);
      DDI_Sync(3081);

      if(me == 0) {
         timings = (struct rusage *) Malloc(np*sizeof(struct rusage));
         getrusage(RUSAGE_SELF,timings);
         gettimeofday(&wall_total,NULL);
      } else {
         getrusage(RUSAGE_SELF,&mycputime);
         timings = &mycputime;
      }

      timings->ru_utime.tv_sec  -= gv(cpu_timer).ru_utime.tv_sec;
      timings->ru_utime.tv_usec -= gv(cpu_timer).ru_utime.tv_usec;
      if(timings->ru_utime.tv_usec < 0) {
         timings->ru_utime.tv_sec--;
         timings->ru_utime.tv_usec += 1000000;
      }

      timings->ru_stime.tv_sec  -= gv(cpu_timer).ru_stime.tv_sec;
      timings->ru_stime.tv_usec -= gv(cpu_timer).ru_stime.tv_usec;
      if(timings->ru_stime.tv_usec < 0) {
         timings->ru_stime.tv_sec--;
         timings->ru_stime.tv_usec += 1000000;
      } 

      wall_total.tv_sec  -= gv(wall_timer).tv_sec;
      wall_total.tv_usec -= gv(wall_timer).tv_usec;
      if(wall_total.tv_usec < 0) {
         wall_total.tv_sec--;
         wall_total.tv_usec += 1000000;
      }

      if(me == 0) {
         for(i=1; i<np; i++) DDI_Recv(&timings[i],sizeof(struct rusage),i);

         fprintf(stdout,"\n ------------------------------------------------");
         fprintf(stdout,"\n CPU timing information for all compute processes");
         fprintf(stdout,"\n ================================================");

         for(i=0; i<np; i++) {


            cpu_total.tv_sec  = timings[i].ru_utime.tv_sec  + timings[i].ru_stime.tv_sec;
            cpu_total.tv_usec = timings[i].ru_utime.tv_usec + timings[i].ru_stime.tv_usec;
            if(cpu_total.tv_usec > 1000000) {
               cpu_total.tv_sec++;
               cpu_total.tv_usec -= 1000000;
            }

            fprintf(stdout,"\n %4i: %d.%.6d + %d.%.6d = %d.%.6d",i,
               (int)timings[i].ru_utime.tv_sec,(int)timings[i].ru_utime.tv_usec,
               (int)timings[i].ru_stime.tv_sec,(int)timings[i].ru_stime.tv_usec,
               (int)cpu_total.tv_sec,(int)cpu_total.tv_usec);
         }

         fprintf(stdout,"\n Wall: %d.%.6d", 
                (int) wall_total.tv_sec, (int) wall_total.tv_usec);
         fprintf(stdout,"\n ================================================\n\n");
                  
         fflush(stdout);
         free(timings);

      } else {

         DDI_Send(&mycputime,sizeof(struct rusage),0);

      }

      DDI_Sync(3082);
   }

# if defined WINDOWS

/*  Windows doesn't provide the usual Unix function 'gettimeofday'  */

# include <time.h> 
# include <Windows.h>
int gettimeofday(struct timeval *tv, struct timezone *tz)
{
  FILETIME ft;
  unsigned __int64 tmpres = 0;
  static int tzflag;
 
  if (NULL != tv)
  {
    GetSystemTimeAsFileTime(&ft);
 
    tmpres |= ft.dwHighDateTime;
    tmpres <<= 32;
    tmpres |= ft.dwLowDateTime;
 
    /*converting file time to unix epoch*/
    tmpres -= DELTA_EPOCH_IN_MICROSECS; 
    tmpres /= 10;  /*convert into microseconds*/
    tv->tv_sec  = (long)(tmpres / 1000000UL);
    tv->tv_usec = (long)(tmpres % 1000000UL);
  }
 
  if (NULL != tz)
  {
    if (!tzflag)
    {
      _tzset();
      tzflag++;
    }
    tz->tz_minuteswest = _timezone / 60;
    tz->tz_dsttime = _daylight;
  }
 
  return 0;
}
# endif
