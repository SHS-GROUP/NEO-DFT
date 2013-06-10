@REM
@REM May   05, 2011 :: Sarom Sok :: File assignments for GAMESS version 1 Oct. 2010 R3
@REM                                ERICFMT, MCPPATH, BASPATH, QUANPOL
@REM
@REM Oct   19, 2010 :: Sarom Sok :: File assignments for GAMESS version 1 Oct. 2010 R1
@REM                                COSCAV, COSDATA, COSPOT
@REM                                RIVMAT, RIT2A, RIT3A, RIT2B, RIT3B
@REM
@REM May    8, 2010 :: Sarom Sok :: Minimized the amount of polling by adding the following
@REM                                flag to the mpiexec call:
@REM                                -env MPICH_DISABLE_SHM 1
@REM
@REM March 31, 2010 :: Sarom Sok :: Windows Server 2008 HPC R2 Edition (beta release)
@REM                                batch file version of the rungms script.
@REM
@REM We make sure that all environmental variables set in this batch file are
@REM done locally and not globally.  We ensure this with SETLOCAL.
@REM
@SETLOCAL
@REM
@REM If you are interested in knowing what is happening in this batch file
@REM then set DEBUG=TRUE.  If SUPRESS=TRUE this flag will override the
@REM printing of the job's environmental variable.
@REM [DEFAULT: FALSE] = NO DEBUGGING INFORMATION FROM THIS SCRIPT
@REM
@SET DEBUG=FALSE
@REM
@REM Control if we print out the GAMESS ASCII text art below.
@REM [DEFAULT: TRUE] = PRINT OUT GAMESS BANNER
@REM
@SET HEADER=TRUE
@REM
@REM Control if we supress the GAMESS environmental variable print out.
@REM The environmental variable will always be made available in the
@REM file location specified by the %ENVFIL% variable in this script.
@REM [DEFAULT: FALSE] = OUTPUT THE ENVIRONMENTAL VARIABLE SETTINGS TO THE SCREEN
@REM
@SET SUPRESS=FALSE
@REM
@REM Control if we want to prevent the script from warning about restart
@REM files from a previous run.  This is done by deleting them before
@REM the script checks for them.
@REM [DEFAULT: FALSE] = DO NOT DELETE PREVIOUS RESTART FILES IF FOUND.
@REM
@SET ERASEOLDFILES=FALSE
@REM
@REM Control if we want to remove %ENVFIL% and %PROCFILE% used by this job
@REM at the end of the run.
@REM [DEFAULT: TRUE] = DELETE %ENVFIL% AND %PROCFILE% AT THE END OF THE RUN
@REM
@SET CLEANALLLOCAL=TRUE
@REM
@REM Controls if we will be redirecting the output of the PURE GAMESS output
@REM a variables defined by %LOGFILE%  This value by default is FALSE.  The
@REM variable is set by the prescence of the 5th argument passed into this
@REM script via the command line.  So there is no need in changing the value
@REM of this variable.
@REM [DEFAULT:FALSE] = NO REDIRECTION. OUTPUT WILL BE TO THE SCREEN FOR PPN=0 OR
@REM                   TO %JOB%.log IF PPN>0.
@REM
@SET REDIRECTION=FALSE
@REM
@REM We need this parameter list!  If it don't exist then teach the user how
@REM to make one.  I hate error messages that don't provide a solution.
@REM
@REM Future development considerations in this script include:
@REM [1] The ability to allow the user to use the parameters.hpc.gms file to
@REM     over-write the values set in this script ^(and eventually
@REM     exported to a the job's environmental variable list^).
@REM
@IF NOT EXIST parameters.hpc.gms (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: parameters.hpc.gms file not found. Please make one.
  @ECHO =========================================================================
  @ECHO.
  @ECHO Create a new txt file called parameters.hpc.gms in this directory.
  @ECHO This file will contain information about your GAMESS setup.
  @ECHO.
  @ECHO The format of this file is:
  @ECHO VARIABLE=VALUE
  @ECHO.
  @ECHO The file must contain:
  @ECHO GAMESSDIR=C:\Path\to\this\directory
  @ECHO UNC_GAMESSDIR=\\server\Path\to\this\directory
  @ECHO UNC_AUXDATADIR=\\server\Path\to\auxdata\directory
  @ECHO.
  @ECHO GAMESSDIR is the path to the GAMESS directory.  You should know it
  @ECHO since you just executed the 'rungms.bat' file that is contained in
  @ECHO this directory.
  @ECHO.
  @ECHO UNC_GAMESSDIR is the path to the GAMESS directory in UNC format.
  @ECHO.
  @ECHO Additional parameters:
  @ECHO         DEBUG=[TRUE,FALSE] Option for extra printing from this script.
  @ECHO        HEADER=[TRUE,FALSE] Option to 'show'/'not show' the GAMESS banner.
  @ECHO       SUPRESS=[TRUE,FALSE] Option to 'not show'/'show' the job's
  @ECHO                            environmental variable list.
  @ECHO ERASEOLDFILES=[TRUE,FALSE] Option to allow the script to remove old restart
  @ECHO                            files if they are found.
  @ECHO CLEANALLLOCAL=[TRUE,FALSE] Option to remove the job's environmental variable
  @ECHO                            and MPI config file.
  @ECHO =========================================================================
  @ECHO.
  @ECHO Now exiting.
  @EXIT /B
)
@REM
@REM We now set these environmental variables.
@REM
@FOR /F "tokens=1,2 delims==" %%A IN (parameters.hpc.gms) DO @SET %%A=%%B
@REM
@REM Print out what the user entered to invoke this batch file.
@REM
@IF %DEBUG%==TRUE (
  @ECHO GAMESS USER^> %0 %1 %2 %3 %4
  @ECHO.
  @ECHO         DEBUG=%DEBUG%
  @ECHO        HEADER=%HEADER%
  @ECHO       SUPRESS=%SUPRESS%
  @ECHO ERASEOLDFILES=%ERASEOLDFILES%
  @ECHO CLEANALLLOCAL=%CLEANALLLOCAL%
)
@REM
@REM Lets keep it "OLD SCHOOL" with some ASCII text art.  GAUSSIAN ain't
@REM got nothing on this!
@REM
@IF %HEADER%==TRUE (
  @ECHO -------------------------------------------------------------------------
  @ECHO.
  @ECHO   .g8"""bgd       db      `7MMM.     ,MMF'`7MM"""YMM   .M"""bgd  .M"""bgd
  @ECHO .dP'     `M      ;MM:       MMMb    dPMM    MM    `7  ,MI    "Y ,MI    "Y
  @ECHO dM'       `     ,V^^MML      M YM   ,M MM    MM   d    `MMb.     `MMb.
  @ECHO MM             ,M  `MM      M  Mb  M' MM    MMmmMM      `YMMNq.   `YMMNq.
  @ECHO MM.    `7MMF'  A     MA     M  YM.P'  MM    MM   Y  , .     `MM .     `MM
  @ECHO `Mb.     MM   A'     VML    M  `YM'   MM    MM     ,M Mb     dM Mb     dM
  @ECHO   `"bmmm!GO .YAN.   .KEES!.JML. `'  .JMML..JMMmmmmMMM P"Ybmmd"  P"Ybmmd"
  @ECHO.
)
@ECHO -------------------------------------------------------------------------
@ECHO.
@ECHO       GAMESS for Microsoft Windows Server 2008 HPC Edition Clusters
@ECHO                     http://www.msg.ameslab.gov/gamess
@ECHO.
@ECHO =========================================================================
@REM
@REM We print the contents of parameters.hpc.gms that may affect this script.
@REM
@IF NOT %SUPRESS%==TRUE (
  @ECHO.
  @ECHO Contents of parameters.hpc.gms file:
  @ECHO.
  @FOR /F "tokens=1,2 delims==" %%A IN (parameters.hpc.gms) DO @ECHO %%A=%%B
  @ECHO.
  @ECHO =========================================================================
)
@REM
@REM Let's see if we have acquired GAMESSDIR, UNC_GAMESSDIR
@REM from that last read.
@REM
@IF NOT DEFINED GAMESSDIR (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: You need to define GAMESSDIR in parameters.hpc.gms
  @ECHO.
  @ECHO          GAMESSDIR is the path to the GAMESS directory.  You should know
  @ECHO          it since you just executed the 'rungms.bat' file that is
  @ECHO          contained in this directory.
  @ECHO =========================================================================
  @ECHO Now exiting.
  @EXIT /B
)
@IF NOT DEFINED UNC_GAMESSDIR (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: You need to define UNC_GAMESSDIR in parameters.hpc.gms
  @ECHO.
  @ECHO          UNC_GAMESSDIR is the path to the GAMESS directory in UNC format.
  @ECHO =========================================================================
  @ECHO Now exiting.
  @EXIT /B
)
@IF NOT DEFINED UNC_AUXDATADIR (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: You need to define UNC_AUXDATADIR in parameters.hpc.gms
  @ECHO.
  @ECHO          UNC_AUXDATADIR is the path to the AUXDATA directory in UNC format.
  @ECHO =========================================================================
  @ECHO Now exiting.
  @EXIT /B
)
@IF NOT DEFINED LOCALSCRATCH (
  @ECHO -------------------------------------------------------------------------
  @ECHO You did not define LOCALSCRATCH in parameters.hpc.gms
  @ECHO.
  @ECHO LOCALSCRATCH will be set to the default value: 'C:\Temp'
  @ECHO =========================================================================
  @SET  LOCALSCRATCH=C:\Temp
)
@REM
@REM We set some of the PATHS now
@REM
@SET RESTART_=scr
@SET SCRATCH_=tmp
@SET   BATCH_=windows
@SET UNC_BATCHDIR=%UNC_GAMESSDIR%\%BATCH_%
@SET     RESTARTDIR=%GAMESSDIR%\%RESTART_%
@SET UNC_RESTARTDIR=%UNC_GAMESSDIR%\%RESTART_%
@SET     SCRATCHDIR=%GAMESSDIR%\%SCRATCH_%
@SET UNC_SCRATCHDIR=%UNC_GAMESSDIR%\%SCRATCH_%
@REM
@REM Debugging will show what these values are.
@REM
@IF %DEBUG%==TRUE (
  @ECHO [Setting up directory path^(s^)]
  @ECHO.
  @ECHO GAMESSDIR=%GAMESSDIR%
  @ECHO UNC_GAMESSDIR=%UNC_GAMESSDIR%
  @ECHO.
  @ECHO UNC_BATCHDIR=%UNC_BATCHDIR%
  @ECHO.
  @ECHO RESTARTDIR=%RESTARTDIR%
  @ECHO UNC_RESTARTDIR=%UNC_RESTARTDIR%
  @ECHO.
  @ECHO SCRATCHDIR=%SCRATCHDIR%
  @ECHO UNC_SCRATCHDIR=%UNC_SCRATCHDIR%
  @ECHO.
  @ECHO LOCALSCRATCH=%LOCALSCRATCH%
  @ECHO.
  @ECHO LOCALSCRATCHDIR will be set later in the script to: C:\Temp\[JOB].[JOBID]
  @ECHO =========================================================================
)
@REM
@REM At minimum we need the input filename in order to proceed further.
@REM
@IF [%1]==[] (
  @ECHO -------------------------------------------------------------------------
  @ECHO.
  @ECHO Usage: rungms.bat [input] [version] [ncpus] [ppn]
  @ECHO  [input]   = The filename of the input ^(with or without a file extension^)
  @ECHO  [version] = The GAMESS version number                              ^(default: 00^)
  @ECHO  [ncpus]   = The number of compute processes requested for this job ^(default:  1^)
  @ECHO  [ppn]     = The number of compute processes per node               ^(default:  1^)
  @ECHO  [logfile] = If a 5th argument is passed then the output of the GAMESS run
  @ECHO              is redirected from STDOUT to [logfile].  The presence of
  @ECHO              this argument will set the variable REDIRECTION to TRUE.
  @ECHO              This will only work if [ppn]=0
  @ECHO =========================================================================
  @ECHO Now exiting.
  @EXIT /B
) ELSE (
  @REM
  @REM Input filename minus the file extension.  The ~n inside the %1 variable
  @REM allows us to do this for arguments passed in at the command-line.
  @REM
  @SET JOB=%~n1
  @REM
)
@REM
@SET LOGFILE=%JOB%.log
@SET ERRFILE=%JOB%.err
@REM
@REM Now process the remaining arguments if they are passed.
@REM
@REM GAMESS version number.
@REM
@IF [%2]==[] (
  @SET VERSION=00
) ELSE (
  @SET VERSION=%2
)
@REM
@REM Number of requested CPUS.
@REM
@IF [%3]==[] (
  @SET NCPUS=1
) ELSE (
  @SET NCPUS=%3
)
@REM
@REM Number of processors per node.
@REM
@IF [%4]==[] (
  @SET PPN=1
) ELSE (
  @SET PPN=%4
)
@REM
@REM Option to redirect the output from the run to %JOB%.log
@REM
@IF [%5]==[] (
  @SET REDIRECTION=FALSE
) ELSE (
  @SET REDIRECTION=TRUE
  @SET LOGFILE=%5
)
@IF %REDIRECTION%==TRUE (
  @IF NOT %PPN%==0 (
    @ECHO REDIRECTION of the output to %LOGFILE% only works if PPN=0.
    @ECHO You specified PPN=%PPN%.  The output of this run will be placed in
    @ECHO the file %JOB%.log
    @ECHO.
    @SET LOGFILE=%JOB%.log
    @ECHO =========================================================================
  )
)
@REM
@REM We set these values early on so that we can perform the check below.  We
@REM perform these checks early on to avoid bombing-out and leaving a job in the
@REM 'configuring' state.  So we do all possible checks before 'job new' is
@REM invoked.
@REM
@SET  MAKEFP=%UNC_RESTARTDIR%\%JOB%.efp
@SET TRAJECT=%UNC_RESTARTDIR%\%JOB%.trj
@SET RESTART=%UNC_RESTARTDIR%\%JOB%.rst
@SET   PUNCH=%UNC_RESTARTDIR%\%JOB%.dat
@REM
@REM Data left over from a previous run might be precious, stop if found.
@REM
@ECHO.
@ECHO Searching for restart files from a previous %JOB% run.
@ECHO.
@SET ERROR=FALSE
@IF EXIST %PUNCH% (
  @IF %ERASEOLDFILES%==TRUE (
    @DEL /Q %PUNCH%
    @SET ERROR=FALSE
  ) ELSE (
    @ECHO Found old PUNCH   : %PUNCH%
    @SET ERROR=TRUE
  )
)
@IF EXIST %MAKEFP% (
  @IF %ERASEOLDFILES%==TRUE (
    @DEL /Q %MAKEFP%
    @SET ERROR=FALSE
  ) ELSE (
    @ECHO Found old MAKEEFP : %MAKEFP%
    @SET ERROR=TRUE
  )
)
@IF EXIST %TRAJECT% (
  @IF %ERASEOLDFILES%==TRUE (
    @DEL /Q %TRAJECT%
    @SET ERROR=FALSE
  ) ELSE (
    @ECHO Found old TRAJECT : %TRAJECT%
    @SET ERROR=TRUE
  )
)
@IF EXIST %RESTART% (
  @IF %ERASEOLDFILES%==TRUE (
    @DEL /Q %RESTART%
    @SET ERROR=FALSE
  ) ELSE (
    @ECHO Found old RESTART : %RESTART%
    @SET ERROR=TRUE
  )
)
@IF %ERROR%==TRUE (
  @ECHO.
  @ECHO Please save, rename, or erase the file^(s^) ^(shown above^) and re-run the job submission.
  @ECHO Now exiting.
  @EXIT /B
) ELSE (
  @ECHO None found.
  @ECHO.
  @ECHO =========================================================================
)
@REM
@REM Lets make sure the file exists
@REM
@IF NOT EXIST %JOB%.inp (
  @REM
  @REM In the case of EXAMnn jobs, this file might be in the "tests" subdirectory.
  @REM
  @IF NOT EXIST %UNC_GAMESSDIR%\tests\%JOB%.inp (
    @ECHO -------------------------------------------------------------------------
    @ECHO "Oh no you didn't!"
    @ECHO ERROR :: Input file %JOB%.inp does not exist.
    @ECHO =========================================================================
    @ECHO.
    @ECHO This job expected the input file to be in the directory:
    @ECHO %~dp0
    @ECHO Please fix your file name problem, and resubmit.
    @ECHO =========================================================================
    @ECHO Now exiting.
    @EXIT /B
  )
)
@REM
@REM Lets make sure the binary exists
@REM
@IF NOT EXIST %UNC_GAMESSDIR%\gamess.%VERSION%.exe (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: The binary gamess.%VERSION%.exe does not exist
  @ECHO =========================================================================
  @ECHO.
  @ECHO Now exiting.
  @EXIT /B
)
@REM
@REM Set the number of nodes required for the job.
@REM
@IF %PPN%==0 (
  @SET /A NODESNEEDED=1
  @REM
  @REM Give it a fake JOBID for local runs
  @REM
  @SET JOBID=00
) ELSE (
  @SET /A NODESNEEDED=NCPUS/PPN
)
@REM
@REM Lets obtain a JOBID.
@REM
@IF NOT %PPN%==0 (
   job new /jobname:"%JOB%"                          ^
           /numnodes:%NODESNEEDED%                   ^
           %NODERESTRICTION%                         ^
           /jobenv:JOBNAME=%JOB%                     ^
           /jobenv:LOGFILE=%UNC_GAMESSDIR%\%LOGFILE% ^
           > .batch.temp
  @FOR /F "tokens=2 delims=:" %%A IN (.batch.temp) DO @ECHO %%A > .batch.temp
  @FOR /F "tokens=1 delims= " %%A IN (.batch.temp) DO @SET JOBID=%%A
)
@IF NOT %SUPRESS%==TRUE (
  @ECHO.
  @ECHO Job Properties
  @ECHO.
  @ECHO Input Filename               : %JOB%.inp
  @ECHO GAMESS Version               : gamess.%VERSION%.exe
  @ECHO Compute Processes Requested  : %NCPUS%
  @ECHO Compute Processes Per Node   : %PPN% ^(0 = running locally^)
  @ECHO Nodes Required               : %NODESNEEDED%
  @ECHO Node Restriction             : %NODERESTRICTION%
  @ECHO.
  @ECHO =========================================================================
)
@REM
@REM Set up the local scratch directory.
@REM
@IF %PPN%==0 (
  @REM
  @REM For local jobs we use the user's scratch space.
  @REM
  @SET LOCALSCRATCHDIR=%UNC_SCRATCHDIR%
) ELSE (
  @REM
  @REM We now have a JOBID lets define the unique local scratch directory.
  @REM
  @SET LOCALSCRATCHDIR=%LOCALSCRATCH%\%JOB%.%JOBID%
)
@REM
@REM Debugging
@REM
@IF %DEBUG%==TRUE (
  @ECHO [Setting up directory path^(s^)]
  @ECHO.
  @ECHO LOCALSCRATCHDIR=%LOCALSCRATCHDIR%
  @ECHO =========================================================================
)
@ECHO.
@ECHO ------------------------ GAMESS EXECUTION SCRIPT ------------------------
@REM
@REM Get host name.
@REM
@FOR /F "usebackq" %%A IN (`HOSTNAME`) DO @SET HOST=%%A
@REM
@REM Get OS version.
@REM
@FOR /F "usebackq tokens=1 delims=+" %%A IN (`VER`) DO @SET MSVERSION=%%A
@REM
@REM Get free disk space - the result will be DISK="   [number] bytes free".
@REM
@DIR %SCRATCHDIR% | FIND "bytes free" > .batch.temp
@FOR /F "tokens=2 delims=)" %%A IN (.batch.temp) DO @SET DISK=%%A
@REM
@REM Print out some useful information.
@REM
@ECHO This job is running on host %HOST% under operating system
@ECHO %MSVERSION% at %DATE% - %TIME%
@ECHO Available user scratch disk at beginning of the job is %DISK:~2,-5%
@ECHO Local scratch disk location on each compute node is: %LOCALSCRATCHDIR%
@REM
@REM
@REM What did I do in the above line to %DISK%
@REM
@REM ~2 means to skip the first two characters of the string stored in %DISK%
@REM this removes extra spaces from the parsing done to get this variable.
@REM -5 means to skip the last 5 characters of the string stored in %DISK%
@REM this removes the word ' free' from the parsing done to get this variable.
@REM
@REM
@REM We will be passing this job's environmental variables to a file.
@REM GAMESS is setup to use the FGEhack so that these environmental
@REM variables are available locally on each compute node.
@REM
@REM (FGE = File Get Environment [variables])
@REM
@SET  ENVFIL=%UNC_SCRATCHDIR%\%JOB%.%JOBID%.GMS.ENV
@REM
@REM
@REM Copy the input.  If it does not exist then complain and exit.
@REM But that shouldn't happen since we just checked for it's existence already.
@REM
@IF EXIST %JOB%.inp (
  COPY /Y %JOB%.inp %UNC_SCRATCHDIR%\%JOB%.%JOBID%.F05 > NUL
) ELSE (
  @IF EXIST tests\%JOB%.inp (
    COPY /Y %UNC_GAMESSDIR%\tests\%JOB%.inp %UNC_SCRATCHDIR%\%JOB%.%JOBID%.F05 > NUL
  ) ELSE (
    @ECHO -------------------------------------------------------------------------
    @ECHO "Oh no you didn't!"
    @ECHO ERROR :: Input file %JOB%.inp does not exist.
    @ECHO =========================================================================
    @ECHO.
    @ECHO This job expected the input file to be in the directory:
    @ECHO %~dp0
    @ECHO Please fix your file name problem, and resubmit.
    @ECHO =========================================================================
    @ECHO Now exiting.
    @EXIT /B
  )
)
@REM File Assignments.
@REM
@REM All binary files should be put on a node's local disk:
@REM
@REM (%LOCALSCRATCHDIR%),
@REM
@REM for the highest speed access possible.  These .Fxx files are typically
@REM not saved for the next run, but they may be big and/or I/O intensive.
@REM
@REM It is convenient to write ASCII output files (PUNCH, RESTART, TRAJECT,
@REM and MAKEFP) to the user's permanent disk, on your file server:
@REM
@REM (%UNC_RESTARTDIR% - cluster jobs)
@REM (%RESTARTDIR% - local jobs)
@REM
@REM They are small, written only by the master process, and are useful outputs
@REM for further runs.
@REM
@REM Some data files may be read by a run, each is read only once, so
@REM that storage of one (1) copy on your file server is appropriate.
@REM a. The ERICFMT file is the Fm(t) data for ERI computations and is provided
@REM    with GAMESS.  It is essential for correct functioning of 2e- integrals.
@REM b. The MCPPATH is just a path name, the program appends the specific
@REM    file names inside that path before opening.  The model core potentials
@REM    and basis set files come as part of GAMESS's normal distribution tree.
@REM c. The EXTBAS file is user-supplied basis sets.
@REM d. see NEO plug-in code's documentation regarding the NUCBAS file.
@REM        Note that you must edit a+b, but will probably skip c+d.
@REM
@REM               ASCII input files (see explanation above)
@REM
@SET ERICFMT=%UNC_AUXDATADIR%\ericfmt.dat
@SET MCPPATH=%UNC_AUXDATADIR%\MCP
@SET BASPATH=%UNC_AUXDATADIR%\BASES
@SET QUANPOL=%UNC_AUXDATADIR%\QUANPOL
@SET  EXTBAS=/dev/null
@SET  NUCBAS=/dev/null
@REM
@REM  MAKEFP - already set above
@SET   GAMMA=%UNC_RESTARTDIR%\%JOB%.gamma
@REM TRAJECT - already set above
@REM RESTART - already set above
@SET   INPUT=%UNC_SCRATCHDIR%\%JOB%.%JOBID%.F05
@REM   PUNCH - already set above
@SET  AOINTS=%LOCALSCRATCHDIR%\%JOB%.F08
@SET  MOINTS=%LOCALSCRATCHDIR%\%JOB%.F09
@SET DICTNRY=%LOCALSCRATCHDIR%\%JOB%.F10
@SET DRTFILE=%LOCALSCRATCHDIR%\%JOB%.F11
@SET CIVECTR=%LOCALSCRATCHDIR%\%JOB%.F12
@SET CASINTS=%LOCALSCRATCHDIR%\%JOB%.F13
@SET  CIINTS=%LOCALSCRATCHDIR%\%JOB%.F14
@SET  WORK15=%LOCALSCRATCHDIR%\%JOB%.F15
@SET  WORK16=%LOCALSCRATCHDIR%\%JOB%.F16
@SET CSFSAVE=%LOCALSCRATCHDIR%\%JOB%.F17
@SET FOCKDER=%LOCALSCRATCHDIR%\%JOB%.F18
@SET  WORK19=%LOCALSCRATCHDIR%\%JOB%.F19
@SET  DASORT=%LOCALSCRATCHDIR%\%JOB%.F20
@SET DIABDAT=%LOCALSCRATCHDIR%\%JOB%.F21
@SET DFTINTS=%LOCALSCRATCHDIR%\%JOB%.F21
@SET DFTGRID=%LOCALSCRATCHDIR%\%JOB%.F22
@SET  JKFILE=%LOCALSCRATCHDIR%\%JOB%.F23
@SET  ORDINT=%LOCALSCRATCHDIR%\%JOB%.F24
@SET  EFPIND=%LOCALSCRATCHDIR%\%JOB%.F25
@SET PCMDATA=%LOCALSCRATCHDIR%\%JOB%.F26
@SET PCMINTS=%LOCALSCRATCHDIR%\%JOB%.F27
@SET SVPWRK1=%LOCALSCRATCHDIR%\%JOB%.F26
@SET SVPWRK2=%LOCALSCRATCHDIR%\%JOB%.F27
@REM
@SET  COSCAV=%LOCALSCRATCHDIR%\%JOB%.F26
@SET COSDATA=%UNC_RESTARTDIR%\%JOB%.cosmo
@SET  COSPOT=%UNC_RESTARTDIR%\%JOB%.pot
@REM
@SET   MLTPL=%LOCALSCRATCHDIR%\%JOB%.F28
@SET  MLTPLT=%LOCALSCRATCHDIR%\%JOB%.F29
@SET  DAFL30=%LOCALSCRATCHDIR%\%JOB%.F30
@SET  SOINTX=%LOCALSCRATCHDIR%\%JOB%.F31
@SET  SOINTY=%LOCALSCRATCHDIR%\%JOB%.F32
@SET  SOINTZ=%LOCALSCRATCHDIR%\%JOB%.F33
@SET  SORESC=%LOCALSCRATCHDIR%\%JOB%.F34
@REM
@REM 35 is used by RESTART, see above
@REM
@SET GCILIST=%LOCALSCRATCHDIR%\%JOB%.F37
@SET HESSIAN=%LOCALSCRATCHDIR%\%JOB%.F38
@SET QMMMTEI=%LOCALSCRATCHDIR%\%JOB%.F39
@SET SOCCDAT=%LOCALSCRATCHDIR%\%JOB%.F40
@SET  AABB41=%LOCALSCRATCHDIR%\%JOB%.F41
@SET  BBAA42=%LOCALSCRATCHDIR%\%JOB%.F42
@SET  BBBB43=%LOCALSCRATCHDIR%\%JOB%.F43
@SET  REMD  =%LOCALSCRATCHDIR%\%JOB%.F44
@SET  MCQD50=%LOCALSCRATCHDIR%\%JOB%.F50
@SET  MCQD51=%LOCALSCRATCHDIR%\%JOB%.F51
@SET  MCQD52=%LOCALSCRATCHDIR%\%JOB%.F52
@SET  MCQD53=%LOCALSCRATCHDIR%\%JOB%.F53
@SET  MCQD54=%LOCALSCRATCHDIR%\%JOB%.F54
@SET  MCQD55=%LOCALSCRATCHDIR%\%JOB%.F55
@SET  MCQD56=%LOCALSCRATCHDIR%\%JOB%.F56
@SET  MCQD57=%LOCALSCRATCHDIR%\%JOB%.F57
@SET  MCQD58=%LOCALSCRATCHDIR%\%JOB%.F58
@SET  MCQD59=%LOCALSCRATCHDIR%\%JOB%.F59
@SET  MCQD60=%LOCALSCRATCHDIR%\%JOB%.F60
@SET  MCQD61=%LOCALSCRATCHDIR%\%JOB%.F61
@SET  MCQD62=%LOCALSCRATCHDIR%\%JOB%.F62
@SET  MCQD63=%LOCALSCRATCHDIR%\%JOB%.F63
@SET  MCQD64=%LOCALSCRATCHDIR%\%JOB%.F64
@SET NMRINT1=%LOCALSCRATCHDIR%\%JOB%.F61
@SET NMRINT2=%LOCALSCRATCHDIR%\%JOB%.F62
@SET NMRINT3=%LOCALSCRATCHDIR%\%JOB%.F63
@SET NMRINT4=%LOCALSCRATCHDIR%\%JOB%.F64
@SET NMRINT5=%LOCALSCRATCHDIR%\%JOB%.F65
@SET NMRINT6=%LOCALSCRATCHDIR%\%JOB%.F66
@SET DCPHFH2=%LOCALSCRATCHDIR%\%JOB%.F67
@SET DCPHF21=%LOCALSCRATCHDIR%\%JOB%.F68
@SET ELNUINT=%LOCALSCRATCHDIR%\%JOB%.F67
@SET NUNUINT=%LOCALSCRATCHDIR%\%JOB%.F68
@SET   GVVPT=%LOCALSCRATCHDIR%\%JOB%.F69
@SET NUMOIN =%LOCALSCRATCHDIR%\%JOB%.F69
@SET NUMOCAS=%LOCALSCRATCHDIR%\%JOB%.F70
@SET NUELMO =%LOCALSCRATCHDIR%\%JOB%.F71
@SET NUELCAS=%LOCALSCRATCHDIR%\%JOB%.F72
@REM
@REM Next files are for RI-MP2
@REM
@SET  RIVMAT=%LOCALSCRATCHDIR%\%JOB%.F51
@SET   RIT2A=%LOCALSCRATCHDIR%\%JOB%.F52
@SET   RIT3A=%LOCALSCRATCHDIR%\%JOB%.F53
@SET   RIT2B=%LOCALSCRATCHDIR%\%JOB%.F54
@SET   RIT3B=%LOCALSCRATCHDIR%\%JOB%.F55
@REM
@REM Next files are for GMCQDPT
@REM
@SET  GMCREF=%LOCALSCRATCHDIR%\%JOB%.F70
@SET  GMCO2R=%LOCALSCRATCHDIR%\%JOB%.F71
@SET  GMCROC=%LOCALSCRATCHDIR%\%JOB%.F72
@SET  GMCOOC=%LOCALSCRATCHDIR%\%JOB%.F73
@SET  GMCCC0=%LOCALSCRATCHDIR%\%JOB%.F74
@SET  GMCHMA=%LOCALSCRATCHDIR%\%JOB%.F75
@SET  GMCEI1=%LOCALSCRATCHDIR%\%JOB%.F76
@SET  GMCEI2=%LOCALSCRATCHDIR%\%JOB%.F77
@SET  GMCEOB=%LOCALSCRATCHDIR%\%JOB%.F78
@SET  GMCEDT=%LOCALSCRATCHDIR%\%JOB%.F79
@SET  GMCERF=%LOCALSCRATCHDIR%\%JOB%.F80
@SET  GMCHCR=%LOCALSCRATCHDIR%\%JOB%.F81
@SET  GMCGJK=%LOCALSCRATCHDIR%\%JOB%.F82
@SET  GMCGAI=%LOCALSCRATCHDIR%\%JOB%.F83
@SET  GMCGEO=%LOCALSCRATCHDIR%\%JOB%.F84
@SET  GMCTE1=%LOCALSCRATCHDIR%\%JOB%.F85
@SET  GMCTE2=%LOCALSCRATCHDIR%\%JOB%.F86
@SET  GMCHEF=%LOCALSCRATCHDIR%\%JOB%.F87
@SET  GMCMOL=%LOCALSCRATCHDIR%\%JOB%.F88
@SET  GMCMOS=%LOCALSCRATCHDIR%\%JOB%.F89
@SET  GMCWGT=%LOCALSCRATCHDIR%\%JOB%.F90
@SET  GMCRM2=%LOCALSCRATCHDIR%\%JOB%.F91
@SET  GMCRM1=%LOCALSCRATCHDIR%\%JOB%.F92
@SET  GMCR00=%LOCALSCRATCHDIR%\%JOB%.F93
@SET  GMCRP1=%LOCALSCRATCHDIR%\%JOB%.F94
@SET  GMCRP2=%LOCALSCRATCHDIR%\%JOB%.F95
@SET  GMCVEF=%LOCALSCRATCHDIR%\%JOB%.F96
@SET  GMCDIN=%LOCALSCRATCHDIR%\%JOB%.F97
@SET  GMC2SZ=%LOCALSCRATCHDIR%\%JOB%.F98
@SET  GMCCCS=%LOCALSCRATCHDIR%\%JOB%.F99
@REM
@REM Next files are used only during closed shell coupled cluster runs.
@REM
@SET  CCREST=%LOCALSCRATCHDIR%\%JOB%.F70
@SET  CCDIIS=%LOCALSCRATCHDIR%\%JOB%.F71
@SET  CCINTS=%LOCALSCRATCHDIR%\%JOB%.F72
@SET CCT1AMP=%LOCALSCRATCHDIR%\%JOB%.F73
@SET CCT2AMP=%LOCALSCRATCHDIR%\%JOB%.F74
@SET CCT3AMP=%LOCALSCRATCHDIR%\%JOB%.F75
@SET    CCVM=%LOCALSCRATCHDIR%\%JOB%.F76
@SET    CCVE=%LOCALSCRATCHDIR%\%JOB%.F77
@SET CCQUADS=%LOCALSCRATCHDIR%\%JOB%.F78
@SET QUADSVO=%LOCALSCRATCHDIR%\%JOB%.F79
@SET EOMSTAR=%LOCALSCRATCHDIR%\%JOB%.F80
@SET EOMVEC1=%LOCALSCRATCHDIR%\%JOB%.F81
@SET EOMVEC2=%LOCALSCRATCHDIR%\%JOB%.F82
@SET  EOMHC1=%LOCALSCRATCHDIR%\%JOB%.F83
@SET  EOMHC2=%LOCALSCRATCHDIR%\%JOB%.F84
@SET EOMHHHH=%LOCALSCRATCHDIR%\%JOB%.F85
@SET EOMPPPP=%LOCALSCRATCHDIR%\%JOB%.F86
@SET EOMRAMP=%LOCALSCRATCHDIR%\%JOB%.F87
@SET EOMRTMP=%LOCALSCRATCHDIR%\%JOB%.F88
@SET EOMDG12=%LOCALSCRATCHDIR%\%JOB%.F89
@SET    MMPP=%LOCALSCRATCHDIR%\%JOB%.F90
@SET   MMHPP=%LOCALSCRATCHDIR%\%JOB%.F91
@SET MMCIVEC=%LOCALSCRATCHDIR%\%JOB%.F92
@SET MMCIVC1=%LOCALSCRATCHDIR%\%JOB%.F93
@SET MMCIITR=%LOCALSCRATCHDIR%\%JOB%.F94
@SET  EOMVL1=%LOCALSCRATCHDIR%\%JOB%.F95
@SET  EOMVL2=%LOCALSCRATCHDIR%\%JOB%.F96
@SET EOMLVEC=%LOCALSCRATCHDIR%\%JOB%.F97
@SET  EOMHL1=%LOCALSCRATCHDIR%\%JOB%.F98
@SET  EOMHL2=%LOCALSCRATCHDIR%\%JOB%.F99
@SET  CCVVVV=%LOCALSCRATCHDIR%\%JOB%.F80
@REM
@REM Next files are used only during open shell coupled cluster runs.
@REM
@SET AMPROCC=%LOCALSCRATCHDIR%\%JOB%.F70
@SET ITOPNCC=%LOCALSCRATCHDIR%\%JOB%.F71
@SET FOCKMTX=%LOCALSCRATCHDIR%\%JOB%.F72
@SET  LAMB23=%LOCALSCRATCHDIR%\%JOB%.F73
@SET   VHHAA=%LOCALSCRATCHDIR%\%JOB%.F74
@SET   VHHBB=%LOCALSCRATCHDIR%\%JOB%.F75
@SET   VHHAB=%LOCALSCRATCHDIR%\%JOB%.F76
@SET    VMAA=%LOCALSCRATCHDIR%\%JOB%.F77
@SET    VMBB=%LOCALSCRATCHDIR%\%JOB%.F78
@SET    VMAB=%LOCALSCRATCHDIR%\%JOB%.F79
@SET    VMBA=%LOCALSCRATCHDIR%\%JOB%.F80
@SET  VHPRAA=%LOCALSCRATCHDIR%\%JOB%.F81
@SET  VHPRBB=%LOCALSCRATCHDIR%\%JOB%.F82
@SET  VHPRAB=%LOCALSCRATCHDIR%\%JOB%.F83
@SET  VHPLAA=%LOCALSCRATCHDIR%\%JOB%.F84
@SET  VHPLBB=%LOCALSCRATCHDIR%\%JOB%.F85
@SET  VHPLAB=%LOCALSCRATCHDIR%\%JOB%.F86
@SET  VHPLBA=%LOCALSCRATCHDIR%\%JOB%.F87
@SET    VEAA=%LOCALSCRATCHDIR%\%JOB%.F88
@SET    VEBB=%LOCALSCRATCHDIR%\%JOB%.F89
@SET    VEAB=%LOCALSCRATCHDIR%\%JOB%.F90
@SET    VEBA=%LOCALSCRATCHDIR%\%JOB%.F91
@SET   VPPPP=%LOCALSCRATCHDIR%\%JOB%.F92
@SET INTERM1=%LOCALSCRATCHDIR%\%JOB%.F93
@SET INTERM2=%LOCALSCRATCHDIR%\%JOB%.F94
@SET INTERM3=%LOCALSCRATCHDIR%\%JOB%.F95
@SET ITSPACE=%LOCALSCRATCHDIR%\%JOB%.F96
@SET INSTART=%LOCALSCRATCHDIR%\%JOB%.F97
@SET  ITSPC3=%LOCALSCRATCHDIR%\%JOB%.F98
@REM
@REM Next files are used only during elongation method runs.
@REM
@FINDSTR /I /R /C:"NELONG=" %UNC_SCRATCHDIR%\%JOB%.%JOBID%.F05
@IF NOT ERRORLEVEL 1 (
    @IF [%4]==[] (
      @SET ELGNAME=ELGFILE
    ) ELSE (
      @SET ELGNAME=%4
    )
)
@IF NOT ERRORLEVEL 1 (
  @SET AOINTS=%LOCALSCRATCHDIR%\%ELGNAME%.F08
  @SET ELGDOS=%UNC_RESTARTDIR%\%JOB%.ldos
  @SET ELGDAT=%LOCALSCRATCHDIR%\%ELGNAME%.F71
  @SET ELGPAR=%LOCALSCRATCHDIR%\%ELGNAME%.F72
  @SET ELGCUT=%LOCALSCRATCHDIR%\%ELGNAME%.F74
  @SET ELGVEC=%LOCALSCRATCHDIR%\%ELGNAME%.F75
  @SET EGINTA=%LOCALSCRATCHDIR%\%ELGNAME%.F77
  @SET EGINTB=%LOCALSCRATCHDIR%\%ELGNAME%.F78
  @SET EGTDHF=%LOCALSCRATCHDIR%\%ELGNAME%.F79
  @SET EGTEST=%LOCALSCRATCHDIR%\%ELGNAME%.F80
)
@REM
@REM Next files are used only during extended TDHF package runs.
@REM
@SET  OLI201=%LOCALSCRATCHDIR%\%JOB%.F201
@SET  OLI202=%LOCALSCRATCHDIR%\%JOB%.F202
@SET  OLI203=%LOCALSCRATCHDIR%\%JOB%.F203
@SET  OLI204=%LOCALSCRATCHDIR%\%JOB%.F204
@SET  OLI205=%LOCALSCRATCHDIR%\%JOB%.F205
@SET  OLI206=%LOCALSCRATCHDIR%\%JOB%.F206
@SET  OLI207=%LOCALSCRATCHDIR%\%JOB%.F207
@SET  OLI208=%LOCALSCRATCHDIR%\%JOB%.F208
@SET  OLI209=%LOCALSCRATCHDIR%\%JOB%.F209
@SET  OLI210=%LOCALSCRATCHDIR%\%JOB%.F210
@SET  OLI211=%LOCALSCRATCHDIR%\%JOB%.F211
@SET  OLI212=%LOCALSCRATCHDIR%\%JOB%.F212
@SET  OLI213=%LOCALSCRATCHDIR%\%JOB%.F213
@SET  OLI214=%LOCALSCRATCHDIR%\%JOB%.F214
@SET  OLI215=%LOCALSCRATCHDIR%\%JOB%.F215
@SET  OLI216=%LOCALSCRATCHDIR%\%JOB%.F216
@SET  OLI217=%LOCALSCRATCHDIR%\%JOB%.F217
@SET  OLI218=%LOCALSCRATCHDIR%\%JOB%.F218
@SET  OLI219=%LOCALSCRATCHDIR%\%JOB%.F219
@SET  OLI220=%LOCALSCRATCHDIR%\%JOB%.F220
@SET  OLI221=%LOCALSCRATCHDIR%\%JOB%.F221
@SET  OLI222=%LOCALSCRATCHDIR%\%JOB%.F222
@SET  OLI223=%LOCALSCRATCHDIR%\%JOB%.F223
@SET  OLI224=%LOCALSCRATCHDIR%\%JOB%.F224
@SET  OLI225=%LOCALSCRATCHDIR%\%JOB%.F225
@SET  OLI226=%LOCALSCRATCHDIR%\%JOB%.F226
@SET  OLI227=%LOCALSCRATCHDIR%\%JOB%.F227
@SET  OLI228=%LOCALSCRATCHDIR%\%JOB%.F228
@SET  OLI229=%LOCALSCRATCHDIR%\%JOB%.F229
@SET  OLI230=%LOCALSCRATCHDIR%\%JOB%.F230
@SET  OLI231=%LOCALSCRATCHDIR%\%JOB%.F231
@SET  OLI232=%LOCALSCRATCHDIR%\%JOB%.F232
@SET  OLI233=%LOCALSCRATCHDIR%\%JOB%.F233
@SET  OLI234=%LOCALSCRATCHDIR%\%JOB%.F234
@SET  OLI235=%LOCALSCRATCHDIR%\%JOB%.F235
@SET  OLI236=%LOCALSCRATCHDIR%\%JOB%.F236
@SET  OLI237=%LOCALSCRATCHDIR%\%JOB%.F237
@SET  OLI238=%LOCALSCRATCHDIR%\%JOB%.F238
@SET  OLI239=%LOCALSCRATCHDIR%\%JOB%.F239
@REM
@REM Next files are used only during divide-and-conquer runs
@REM
@SET   DCSUB=%LOCALSCRATCHDIR%\%JOB%.F250
@SET   DCVEC=%LOCALSCRATCHDIR%\%JOB%.F251
@SET   DCEIG=%LOCALSCRATCHDIR%\%JOB%.F252
@SET    DCDM=%LOCALSCRATCHDIR%\%JOB%.F253
@SET   DCDMO=%LOCALSCRATCHDIR%\%JOB%.F254
@SET     DCQ=%LOCALSCRATCHDIR%\%JOB%.F255
@SET     DCW=%LOCALSCRATCHDIR%\%JOB%.F256
@SET   DCEDM=%LOCALSCRATCHDIR%\%JOB%.F257
@REM
@REM Next files are used only during LMO hyperpolarizability analysis
@REM
@SET LHYPWRK=%LOCALSCRATCHDIR%\%JOB%.F297
@SET LHYPWK2=%LOCALSCRATCHDIR%\%JOB%.F298
@SET BONDDPF=%LOCALSCRATCHDIR%\%JOB%.F299
@REM
@REM Next value is used only within the VB2000 add-on code
@REM
@SET VB2000PATH=%UNC_GAMESSDIR%\vb2000
@SET GMSJOBNAME=%JOB%
@REM
@REM Next files are used during EFMO runs
@REM
@FINDSTR /I /R /C:"^ IEFMO=" %UNC_SCRATCHDIR%\%JOB%.%JOBID%.F05 > NUL
@IF NOT ERRORLEVEL 1 (
  @SET EFMOI=%LOCALSCRATCHDIR%\$JOB.F102
  @SET EFMOF=%LOCALSCRATCHDIR%\$JOB.F103
)
@REM
@REM Next files are used only during CIM runs
@REM
@SET CIMFILE=%UNC_RESTARTDIR%\%JOB%.cim
@SET CIMDMN=%UNC_RESTARTDIR%\%JOB%.dmn
@SET CIMDCT=%LOCALSCRATCHDIR%\%JOB%.Fcdt
@SET CIMAOI=%LOCALSCRATCHDIR%\%JOB%.Fcao
@SET CIMMOI=%LOCALSCRATCHDIR%\%JOB%.Fcmo
@REM
@REM "Group DDI" will be set to true, below, if=$GDDI input is found
@REM
@SET GDDIJOB=FALSE
@REM
@FINDSTR /I /R /C:"^ \$GDDI" %UNC_SCRATCHDIR%\%JOB%.%JOBID%.F05 > NUL
@IF NOT ERRORLEVEL 1 (
  @SET GDDIJOB=TRUE
  @REM
  @ECHO.
  @ECHO This is a GDDI run, keeping output files on local disks
  @ECHO until the very end of the run, when they'll be saved from
  @ECHO the master process in the first group, only.
  @REM
  @SET  MAKEFP=%LOCALSCRATCHDIR%\%JOB%.F01
  @SET TRAJECT=%LOCALSCRATCHDIR%\%JOB%.F04
  @SET  OUTPUT=%LOCALSCRATCHDIR%\%JOB%.F06
  @SET   PUNCH=%LOCALSCRATCHDIR%\%JOB%.F07
  @SET RESTART=%LOCALSCRATCHDIR%\%JOB%.F35
  @REM
  @IF NOT %PPN%==0 (
    @REM
    @REM If we are running across multiple nodes then we will need to update
    @REM the environmental variable for the location of the job's INPUT file
    @REM which is located on the compute node's local scratch directory.
    @REM
    @SET  INPUT=%LOCALSCRATCHDIR%\%JOB%.F05
  )
  @REM
) ELSE (
  @SET GDDIJOB=FALSE
)
@REM
@REM Store the job's environmental variable to a file.
@REM
@SET > %ENVFIL%
@REM
@REM Echo out the job's environmental variable if we are not SUPPRESSing
@REM the output.
@REM
@IF NOT %SUPRESS%==TRUE (
  @ECHO.
  @CALL %UNC_BATCHDIR%\GMSENV.HPC.BAT
) ELSE (
  @REM
  @REM Maybe the DEBUG flag is on.
  @REM
  @IF %DEBUG%==TRUE (
    @ECHO.
    @CALL %UNC_BATCHDIR%\GMSENV.HPC.BAT
  )
)
@ECHO.
@ECHO =========================================================================
@REM
@REM We don't need to have a HOSTFILE because the Windows HPC Cluster Manager takes
@REM care of those particulars.
@REM
@SET PROCFILE=%UNC_SCRATCHDIR%\%JOB%.%JOBID%.processes.mpd
@IF EXIST %PROCFILE% DEL /Q %PROCFILE%
@REM
@REM Allow for compute process and data server (one pair per core). (We use this for PPN=0)
@REM
@IF %PPN%==0 (
  @SET /A NPROCS=NCPUS+NCPUS
)
@IF %PPN%==0 (
  @ECHO -env ENVFIL %ENVFIL% -n %NPROCS% %UNC_GAMESSDIR%\gamess.%VERSION%.exe >> %PROCFILE%
)
@REM
@REM Running on the cluster
@REM
@IF NOT %PPN%==0 (
  @REM
  @REM Get maxnodes from parsing the list of 'Online' nodes from the command `node list`.
  @REM
  @SET /A MAXNODES=0
  @FOR /F "usebackq tokens=2 delims= " %%A IN (`node list`) DO @IF %%A==Online @SET /A MAXNODES+=1
  @REM
  @REM For each compute process we need a data server
  @REM
  @SET /A PPN2=PPN+PPN
  @REM
) ELSE (
  @SET /A MAXNODES=1
)
@IF NOT %PPN%==0 (
  @IF %NODESNEEDED% GTR %MAXNODES% (
    @ECHO.
    @ECHO   This run requires %NODESNEEDED% nodes.
    @ECHO   Only %MAXNODES% nodes are ONLINE.
    @ECHO.
    @ECHO   Your job might not run immediately
    @ECHO   and just sit in the queue forever!
    @ECHO.
    @ECHO =========================================================================
  )
  @FOR /L %%A IN (1,1,%NODESNEEDED%) DO @ECHO -env MPICH_DISABLE_SHM 1 -n %PPN2% %UNC_GAMESSDIR%\gamess.%VERSION%.exe >> %PROCFILE%
)
@IF %DEBUG%==TRUE (
  @ECHO [Parameters]
  @IF %PPN%==0 (
    @ECHO NCPUS            : %NCPUS%
    @ECHO NPROCS           : %NPROCS% ^(Should be 2x the value of NCPUS^)
    @ECHO MAXNODES         : %MAXNODES%
    @ECHO.
  ) ELSE (
    @ECHO NCPUS            : %NCPUS%
    @ECHO PPN              : %PPN%
    @ECHO PPN2             : %PPN2% ^(Should be 2x the value of PPN^)
    @ECHO NODESNEEDED      : %NODESNEEDED%
    @ECHO MAXNODES         : %MAXNODES% ^(Number of nodes with status 'Online'^)
    @ECHO Node Restriction : %NODERESTRICTION%
    @ECHO.
  )
  @ECHO [Contents of procfile]
  @ECHO %PROCFILE%
  @ECHO.
  @TYPE %PROCFILE%
  @ECHO =========================================================================
)
@REM
@REM  At this point we are ready to run!
@REM
@ECHO Microsoft MPI ^(MS-MPI^) will be running GAMESS on %NODESNEEDED% node^(s^).
@ECHO The binary will be kicked off by 'mpiexec' is gamess.%VERSION%.exe
@ECHO MS-MPI will run %NCPUS% compute process^(es^) and %NCPUS% data server^(s^).
@ECHO.
@IF %PPN%==0 (
  @REM
  @REM We are running locally.
  @REM
  @IF %GDDIJOB%==TRUE (
    @ECHO.
    @ECHO This is a GDDI run the output will be dumped to this screen when
    @ECHO the run has completed which may be a while. But since this is
    @ECHO running on the current computer you can tail the progress of
    @ECHO the following file:
    @ECHO.
    @ECHO %SCRATCHDIR%\%JOB%.F06
    @ECHO.
  )
  @IF %REDIRECTION%==TRUE (
    @mpiexec -configfile %PROCFILE% >%LOGFILE%
  ) ELSE (
    @mpiexec -configfile %PROCFILE%
  )
) ELSE (
  @REM
  @REM We are running on the cluster.
  @REM
  @REM  We setup the job with tasks.
  @REM
  @REM  Task :: Get the compute nodes ready.
  @REM
  IF %GDDIJOB%==TRUE (
    @REM
    @REM  This NodePrep task will create a %JOB%.%JOBID% directory in each
    @REM  compute node's local scratch directory.
    @REM
    @REM  And because this is a GDDI job this task will copy the
    @REM  %JOB.%JOBID%.F05 (job's input file) from the user's local scratch
    @REM  directory to %JOB%.F05 in the local scratch directry on every
    @REM  compute node.
    @REM
    @ECHO =========================================================================
    @ECHO Setting up NodePrep task for GAMESS run using GDDI.
    @REM
    @job add  %JOBID% /name:"GAMESSGDDIJOBNodePreparation"   ^
                      -Type:NodePrep                         ^
                      /env:JOBID=%JOBID%                     ^
                      /env:LOCALDIR=%LOCALSCRATCHDIR%        ^
                      /env:WORKDIR=%UNC_SCRATCHDIR%          ^
     %UNC_BATCHDIR%\gamess.start.gddi.bat
    @REM
    @ECHO =========================================================================
  ) ELSE (
    @REM
    @REM  This isn't a GDDI job so just create a %JOB%.%JOBID% directory in
    @REM  each compute nodes local scratch directory.
    @REM
    @ECHO =========================================================================
    @ECHO Setting up NodePrep task for GAMESS run.
    @REM
     job add  %JOBID% /name:"GAMESSJobNodePreparation"       ^
                      -Type:NodePrep                         ^
                      /env:JOBID=%JOBID%                     ^
                      /env:LOCALDIR=%LOCALSCRATCHDIR%        ^
                      /env:WORKDIR=%UNC_SCRATCHDIR%          ^
     %UNC_BATCHDIR%\gamess.start.bat
    @REM
    @ECHO =========================================================================
  )
  @REM
  @REM  Task :: Run the job.
  @REM
  @ECHO.
  @ECHO =========================================================================
  @ECHO Setting up RunGAMESS task for GAMESS run.
  @REM
   job add  %JOBID% /name:"RunGAMESS"                        ^
                    /numnodes:%NODESNEEDED%                  ^
                    /workdir:%UNC_GAMESSDIR%                 ^
                    /stdout:%LOGFILE%                        ^
                    /stderr:%ERRFILE%                        ^
                    /env:ENVFIL=%ENVFIL%                     ^
   mpiexec -configfile %PROCFILE%
  @REM
  @ECHO =========================================================================
  @REM
  @REM  Task :: Each compute node will create a disk usage information
  @REM          log file in the user's local scratch directory.
  @REM
  @REM          If disk usage information is not important then you may
  @REM          want to comment out this task.  Especially, if you are
  @REM          performing a parametric sweep job with GAMESS. But you
  @REM          must change the /depend: flag on each task below or else
  @REM          your job won't be submitted.  You will also want to
  @REM          comment out the GAMESSAppendDiskUsage task as well.
  @REM
  @ECHO.
  @ECHO =========================================================================
  @ECHO Setting up GAMESSDiskUsage task for GAMESS run.
  @REM
   job add  %JOBID% /name:"GAMESSDiskUsage"                  ^
                    /depend:"RunGAMESS"                      ^
                    /numnodes:%NODESNEEDED%                  ^
                    /env:JOBID=%JOBID%                       ^
                    /env:LOCALDIR=%LOCALSCRATCHDIR%          ^
                    /env:WORKDIR=%UNC_SCRATCHDIR%            ^
   mpiexec -n %NODESNEEDED% %UNC_BATCHDIR%\gamess.diskusage.bat
  @REM
  @ECHO =========================================================================
  @REM
  @REM  Task :: Append disk usage information to the log file.
  @REM
  IF %GDDIJOB%==TRUE (
    @REM
    @REM  We need to save some files from each compute nodes local scratch
    @REM  to the user's restart directory.
    @REM
    @ECHO.
    @ECHO =========================================================================
    @ECHO Setting up GAMESSGDDIOutput task for GAMESS run using GDDI.
    @REM
     job add  %JOBID% /name:"GAMESSGDDIOutput"               ^
                      /depend:"GAMESSDiskUsage"              ^
                      /numnodes:%NODESNEEDED%                ^
                      /env:JOBID=%JOBID%                     ^
                      /env:LOCALDIR=%LOCALSCRATCHDIR%        ^
                      /env:LOGFILE=%UNC_GAMESSDIR%\%LOGFILE% ^
                      /env:WORKDIR=%UNC_SCRATCHDIR%          ^
     mpiexec -n %NODESNEEDED% %UNC_BATCHDIR%\gamess.gddi.output.bat
    @REM
    @ECHO =========================================================================
    @REM
    @REM  This task is no different from the one directly below it after the
    @REM  'else'.  The only difference is in the task's /depend: flag.
    @REM  This is to unsure that
    @REM
    @ECHO.
    @ECHO =========================================================================
    @ECHO Setting up GAMESSAppendDiskUsage task for GAMESS run using GDDI.
    @REM
     job add  %JOBID% /name:"GAMESSAppendDiskUsage"          ^
                      /depend:"GAMESSGDDIOutput"             ^
                      /numnodes:1                            ^
                      /env:JOBID=%JOBID%                     ^
                      /env:LOGFILE=%UNC_GAMESSDIR%\%LOGFILE% ^
                      /env:LOCALDIR=%LOCALSCRATCHDIR%        ^
                      /env:WORKDIR=%UNC_SCRATCHDIR%          ^
     %UNC_BATCHDIR%\gamess.append.diskusage.bat
    @REM
    @ECHO =========================================================================
  ) ELSE (
    @REM
    @ECHO.
    @ECHO =========================================================================
    @ECHO Setting up GAMESSAppendDiskUsage task for GAMESS run.
    @REM
    @job add  %JOBID% /name:"GAMESSAppendDiskUsage"          ^
                      /depend:"GAMESSDiskUsage"              ^
                      /numnodes:1                            ^
                      /env:JOBID=%JOBID%                     ^
                      /env:LOGFILE=%UNC_GAMESSDIR%\%LOGFILE% ^
                      /env:LOCALDIR=%LOCALSCRATCHDIR%        ^
                      /env:WORKDIR=%UNC_SCRATCHDIR%          ^
     %UNC_BATCHDIR%\gamess.append.diskusage.bat
    @REM
    @ECHO =========================================================================
  )
  @REM
  @REM  Task :: Clean up the users's local scratch directory.  This is
  @REM          done in a folder (default: tmp) within the WORKDIR
  @REM          specified below in the task.  This directory contains the
  @REM          following files:
  @REM
  @REM          %JOB%.%JOBID%.ENV
  @REM          %JOB%.%JOBID%.F05
  @REM        *.%JOB%.%JOBID%.diskusage.log
  @REM          %JOB%.%JOBID%.nodes.mpd
  @REM          %JOB%.%JOBID%.processes.mpd
  @REM
  @REM          This task also merges two files produced during a GDDI run
  @REM          into a single %LOGFILE%.  This is done because the output
  @REM          the in the user's GAMESS directory that is produced during
  @REM          the run only contains the GDDI timing information.  The
  @REM          real output file (%JOB%.F06) is located on the master
  @REM          compute node.  So the merging is done here.  The GDDI
  @REM          conditional is done in the batch file.
  @REM
  @ECHO.
  @ECHO =========================================================================
  @ECHO Setting up CleanGAMESS task for GAMESS run.
  @REM
   job add  %JOBID% /name:"CleanGAMESS"                      ^
                    /depend:"GAMESSAppendDiskUsage"          ^
                    /numnodes:1                              ^
                    /workdir:%UNC_GAMESSDIR%                 ^
                    /env:JOBID=%JOBID%                       ^
                    /env:WORKDIR=%UNC_SCRATCHDIR%            ^
   %UNC_BATCHDIR%\gamess.clean.bat
  @REM
  @ECHO =========================================================================
  @REM
  @REM  Task :: Clean up the compute nodes.
  @REM
  IF %GDDIJOB%==TRUE (
    @REM
    @REM  We need to save some files from each compute nodes local scratch
    @REM  to the user's restart directory.
    @REM
    @ECHO.
    @ECHO =========================================================================
    @ECHO Setting up NodeRelease task for GAMESS run using GDDI.
    @REM
     job add  %JOBID% /name:"GAMESSJobNodeCleanup"           ^
                      -Type:NodeRelease                      ^
                      /env:LOCALDIR=%LOCALSCRATCHDIR%        ^
                      /env:WORKDIR=%UNC_RESTARTDIR%          ^
     %UNC_BATCHDIR%\gamess.finish.gddi.bat
    @REM
    @ECHO =========================================================================
  ) ELSE (
    @REM
    @ECHO.
    @ECHO =========================================================================
    @ECHO Setting up NodeRelease task for GAMESS run.
    @REM
     job add  %JOBID% /name:"GAMESSJobNodeCleanup"           ^
                      -Type:NodeRelease                      ^
                      /env:LOCALDIR=%LOCALSCRATCHDIR%        ^
     %UNC_BATCHDIR%\gamess.finish.bat
    @REM
    @ECHO =========================================================================
  )
  @REM
  @REM  Now, submit the job.
  @REM
   job submit /id:%JOBID%
)
@REM
@REM Tidy up some GDDI stuff if we ran locally (PPN=0)
@REM
@IF %PPN%==0 (
  @IF %GDDIJOB%==TRUE (
    @REM
    @REM If we are using REDIRECTION then do some file juggling for GDDI runs.
    @REM
    @IF %REDIRECTION%==TRUE (
      @TYPE %LOGFILE% > .batch.temp
      @TYPE %OUTPUT% > %LOGFILE%
      @TYPE .batch.temp >> %LOGFILE%
    ) ELSE (
      @TYPE %OUTPUT%
    )
        @COPY %LOCALSCRATCHDIR%\%JOB%.F06* %RESTARTDIR% > NUL
        @COPY %LOCALSCRATCHDIR%\%JOB%.F07  %RESTARTDIR%\%JOB%.dat > NUL
    @IF EXIST %LOCALSCRATCHDIR%\%JOB%.F04 @COPY %LOCALSCRATCHDIR%\%JOB%.F04 %RESTARTDIR%\%JOB%.trj > NUL
    @IF EXIST %LOCALSCRATCHDIR%\%JOB%.F35 @COPY %LOCALSCRATCHDIR%\%JOB%.F35 %RESTARTDIR%\%JOB%.rst > NUL
  )
  @REM
  @IF %REDIRECTION%==TRUE (
    @ECHO            Disk usage ^(in bytes^)       >>%LOGFILE%
    @ECHO  ========================================>>%LOGFILE%
    @DIR %LOCALSCRATCHDIR%\%JOB%.* | FIND "/"      >>%LOGFILE%
    @ECHO  ---------------------------------------->>%LOGFILE%
    @ECHO             Local date and time          >>%LOGFILE%
    @ECHO  ========================================>>%LOGFILE%
    @ECHO %DATE% - %TIME%                          >>%LOGFILE%
  ) ELSE (
    @REM
    @ECHO            Disk usage ^(in bytes^)
    @ECHO  ========================================
    @REM
    @REM Print disk usage without all that extra stuff on top that comes
    @REM with the DIR command.
    @REM
    @DIR %LOCALSCRATCHDIR%\%JOB%.* | FIND "/"
    @REM
    @ECHO  ----------------------------------------
    @ECHO             Local date and time
    @ECHO  ========================================
    @ECHO %DATE% - %TIME%
  )
  @REM
  @REM Time to clean up.
  @REM
  @REM Remove scratch files
  @REM
  @DEL /Q %LOCALSCRATCHDIR%\%JOB%.F*
  @DEL /Q %LOCALSCRATCHDIR%\%JOB%.%JOBID%.F*
  @REM
  @REM Remove temporary file used by this batch to store variables.
  @REM
  @IF %CLEANALLLOCAL%==TRUE (
    @DEL /Q %ENVFIL%
    @DEL /Q %PROCFILE%
  )
)
@IF EXIST .batch.temp @DEL /Q .batch.temp
@REM
@REM FIN
@REM
@EXIT /B
