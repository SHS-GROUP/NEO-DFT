/* ------------------------------------------------------------------ *\
 * Distributed Data Interface
 * ==========================
 *
 * print strings about version numbers and machine characteristics
 *
 * Author: Andrey Asadchev
 * 13 May 10 - SS  - Windows can't do the printing
\* ------------------------------------------------------------------ */

#ifndef WINDOWS
  # include <sys/utsname.h>
#endif
  # include "ddi_base.h"
  # include "ddi_runtime.h"
  # include "ddi_util.h"

void DDI_Get_version(int *version, int *subversion) {
    *version = DDI_VERSION;
    *subversion = DDI_SUBVERSION;
}

/** DDI build string */
void DDI_Get_build(char *build) {

#if defined USE_SYSV
    join(build,"SysV","/");
#endif

#if FULL_SMP
    join(build,"SMP","/");
#endif

#if defined DDI_SOC
    join(build,"Sockets","/");
#endif

#if defined DDI_MPI
#if defined DDI_MPI2
    join(build,"MPI2","/");
#else
    join(build,"MPI","/");
#endif
#endif

#if defined DDI_LAPI
    join(build,"LAPI","/");
#endif

#if defined DDI_ARMCI
    join(build,"ARMCI","/");
#endif

#ifdef DDI_MPI_FILE
    join(build,"MPI-IO","/");
#endif

}

# if !defined WINDOWS

/** @see ddi_runtime.h */
void DDI_Runtime_print(FILE *stream) {
    DDI_RUNTIME_PRINT_(stream);
}

/** @see ddi_runtime.h */
void DDI_POSIX_Runtime_print(FILE *stream) {
    int version, subversion;
    char build[DDI_MAX_BUILD];
    struct utsname uts;

    build[0] = '\0';

    DDI_Get_version(&version, &subversion);
    DDI_Get_build(build);
    uname(&uts);

    fprintf(stream,"DDI %i.%i %s %s %s %s\n",
	    version, subversion, build, uts.sysname, uts.release, uts.machine);
}

#endif
