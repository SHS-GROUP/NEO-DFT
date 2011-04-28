/*      this file used only on Microsoft Windows
 *
 * 13 May 10 - SS  - new file to support Windows/pgcc
 */

 # ifdef WINDOWS
 # include "mysystem.h" 
int getrusage(int who, struct rusage * rusage)
{
  FILETIME        starttime;
  FILETIME        exittime;
  FILETIME        kerneltime;
  FILETIME        usertime;
  ULARGE_INTEGER li;

  if (rusage == (struct rusage *) NULL)
  {
    errno = EFAULT;
    return -1;
  }
  memset(rusage, 0, sizeof(struct rusage));
  if (GetProcessTimes(GetCurrentProcess(),
     &starttime, &exittime, &kerneltime, &usertime) == 0)
  {
    _dosmaperr(GetLastError());
    return -1;
  }

  /* Convert FILETIMEs (0.1 us) to struct timeval */
  memcpy(&li, &kerneltime, sizeof(FILETIME));
  li.QuadPart /= 10L;                     /* Convert to microseconds */
  rusage->ru_stime.tv_sec = li.QuadPart / 1000000L;
  rusage->ru_stime.tv_usec = li.QuadPart % 1000000L;

  memcpy(&li, &usertime, sizeof(FILETIME));
  li.QuadPart /= 10L;                     /* Convert to microseconds */
  rusage->ru_utime.tv_sec = li.QuadPart / 1000000L;
  rusage->ru_utime.tv_usec = li.QuadPart % 1000000L;
  return 0;
}
 # endif
