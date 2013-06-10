@REM
@REM
@REM May   28, 2013 :: Sarom Leang :: File assignments for VB2000, CIM, EFMO
@REM
@REM May   05, 2011 :: Sarom Sok :: File assignments for GAMESS version 1 Oct. 2010 R3
@REM                                ERICFMT, MCPPATH, BASPATH, QUANPOL
@REM
@REM Oct   19, 2010 :: Sarom Sok :: File assignments for GAMESS version 1 Oct. 2010 R1
@REM                                COSCAV, COSDATA, COSPOT
@REM                                RIVMAT, RIT2A, RIT3A, RIT2B, RIT3B
@REM
@REM March 31, 2010 :: Sarom Sok :: Windows batch file version
@REM                                of the rungms script.
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
@REM [1] The ability to allow the user to use the parameters.gms file to
@REM     over-write the values set in this script ^(and eventually
@REM     exported to a the job's environmental variable list^).
@REM
@REM [2] Permit MPI jobs across multiple nodes.
@REM
@IF NOT EXIST parameters.gms (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: parameters.gms file not found. Please make one.
  @ECHO =========================================================================
  @ECHO.
  @ECHO Create a new txt file called parameters.gms in this directory.
  @ECHO This file will contain information about your GAMESS setup.
  @ECHO.
  @ECHO The format of this file is:
  @ECHO VARIABLE=VALUE
  @ECHO.
  @ECHO The file must contain:
  @ECHO GAMESSDIR=C:\Path\to\this\directory
  @ECHO AUXDATADIR=C:\Path\to\your\local\auxdata\directory
  @ECHO RESTARTDIR=C:\Path\to\your\local\restart\directory
  @ECHO SCRATCHDIR=C:\Path\to\your\local\scratch\directory
  @ECHO.
  @ECHO GAMESSDIR is the path to the GAMESS directory.  You should know it
  @ECHO since you just executed the 'rungms.bat' file that is contained in
  @ECHO this directory.
  @ECHO.
  @ECHO RESTARTDIR is the path to the directory where you wish to save
  @ECHO the GAMESS restart files ^(.dat, .rst, .trj, .efp, .gamma, .ldos^)
  @ECHO.
  @ECHO SCRATCHDIR is the path to the directory GAMESS will use
  @ECHO to store scratch files ^(which gets deleted at the end of a run^).
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
@FOR /F "tokens=1,2 delims==" %%A IN (parameters.gms) DO @SET %%A=%%B
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
@ECHO                        GAMESS for Microsoft Windows
@ECHO                     http://www.msg.ameslab.gov/gamess
@ECHO.
@ECHO =========================================================================
@REM
@REM We print the contents of parameters.gms that may affect this script.
@REM
@IF NOT %SUPRESS%==TRUE (
  @ECHO.
  @ECHO Contents of parameters.gms file:
  @ECHO.
  @FOR /F "tokens=1,2 delims==" %%A IN (parameters.gms) DO @ECHO %%A=%%B
  @ECHO.
  @ECHO =========================================================================
)
@REM
@REM Let's see if we have acquired GAMESSDIR, AUXDATADIR, RESTARTDIR, SCRATCHDIR
@REM from that last read.
@REM
@IF NOT DEFINED GAMESSDIR (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: You need to define GAMESSDIR in parameters.gms
  @ECHO =========================================================================
  @ECHO Now exiting.
  @EXIT /B
)
@IF NOT DEFINED AUXDATADIR (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: You need to define AUXDATADIR in parameters.gms
  @ECHO =========================================================================
  @ECHO Now exiting.
  @EXIT /B
)
@IF NOT DEFINED RESTARTDIR (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: You need to define RESTARTDIR in parameters.gms
  @ECHO =========================================================================
  @ECHO Now exiting.
  @EXIT /B
)
@IF NOT DEFINED SCRATCHDIR (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: You need to define SCRATCHDIR in parameters.gms
  @ECHO =========================================================================
  @ECHO Now exiting.
  @EXIT /B
)
@REM
@REM Debugging will show what these values are.
@REM
@IF %DEBUG%==TRUE (
  @ECHO [Setting up directory path^(s^)]
  @ECHO.
  @ECHO  GAMESSDIR=%GAMESSDIR%
  @ECHO AUXDATADIR=%AUXDATADIR%
  @ECHO RESTARTDIR=%RESTARTDIR%
  @ECHO SCRATCHDIR=%SCRATCHDIR%
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
  @ECHO  [version] = The GAMESS version number                      ^(default: 00^)
  @ECHO  [ncpus]   = The number of compute processes requested for this job ^(default:  1^)
  @ECHO  [ppn]     = The number of compute processes per node               ^(default:  0^)
  @ECHO              The default means that we are running a locally.
  @ECHO              This script is not set to run across nodes. The value of
  @ECHO              [ppn] must be 0.  So either pass 0 as your 4th argument or
  @ECHO              leave it blank.  If you wish to run across multiple nodes then
  @ECHO              please use 'rungms.hpc.bat' which is configured specifically
  @ECHO              for Windows Server 2008 HPC Edition clusters.
  @ECHO  [logfile] = If a 5th argument is passed then the output of the GAMESS run is
  @ECHO              redirected from STDOUT to the name of this file.  The presence of
  @ECHO              this argument will set the variable REDIRECTION to TRUE.
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
@REM We set this variable just incase a users decides to play with the %REDIRECTION%
@REM variable.
@REM
@SET LOGFILE=%JOB%.log
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
  @SET PPN=0
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
@REM
@IF NOT %PPN%==0 (
    @ECHO.
    @ECHO This script is not set to run across nodes. The value of
    @ECHO [ppn] must be 0.  So either pass 0 as your 4th argument or
    @ECHO leave it blank.  If you wish to run across multiple nodes then
    @ECHO please consider the Windows Server 2008 HPC distribution of GAMESS.
    @ECHO.
    @ECHO =========================================================================
    @ECHO Now exiting.
    @EXIT /B
)
@REM
@REM Lets make sure the file exists
@REM
@IF NOT EXIST %JOB%.inp (
  @REM
  @REM In the case of EXAMnn jobs, this file might be in the "tests" subdirectory.
  @REM
  @IF NOT EXIST %GAMESSDIR%\tests\%JOB%.inp (
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
  ) ELSE (
    @REM
    @REM The file exists in the test directory. Copy it over to the local scratch.
    @REM
    @COPY /Y %GAMESSDIR%\tests\%JOB%.inp %SCRATCHDIR%\%JOB%.F05 > NUL
  )
) ELSE (
    @REM
    @REM The file exists in the current directory. Copy it over to the local scratch.
    @REM
    @COPY /Y %JOB%.inp %SCRATCHDIR%\%JOB%.F05 > NUL
)
@IF NOT %SUPRESS%==TRUE (
  @ECHO.
  @ECHO Job Properties
  @ECHO.
  @ECHO Input Filename               : %JOB%.inp
  @ECHO GAMESS Version               : gamess.%VERSION%.exe
  @ECHO Compute Processes Requested  : %NCPUS%
  @ECHO Compute Processes Per Node   : %PPN% ^(0 = running locally^)
  @ECHO Nodes Required               : 1
  @ECHO.
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
@ECHO Available scratch disk at beginning of the job is %DISK:~2,-5%
@REM
@REM What did I do in the above line to %DISK%
@REM
@REM ~2 means to skip the first two characters of the string stored in %DISK%
@REM this removes extra spaces from the parsing done to get this variable.
@REM -5 means to skip the last 5 characters of the string stored in %DISK%
@REM this removes the word ' free' from the parsing done to get this variable.
@REM
@REM File Assignments.
@REM
@REM All binary files should be put on a node's local disk (%SCRATCHDIR%),
@REM for the highest speed access possible.  These .Fxx files are typically
@REM not saved for the next run, but they may be big and/or I/O intensive.
@REM
@REM It is convenient to write ASCII output files (PUNCH, RESTART, TRAJECT,
@REM and MAKEFP) to the user's permanent disk, on your file server (%RESTARTDIR%).
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
@SET  ENVFIL=%SCRATCHDIR%\%JOB%.GMS.ENV
@SET ERICFMT=%AUXDATADIR%\ericfmt.dat
@SET MCPPATH=%AUXDATADIR%\MCP
@SET BASPATH=%AUXDATADIR%\BASES
@SET QUANPOL=%AUXDATADIR%\QUANPOL
@SET  EXTBAS=/dev/null
@SET  NUCBAS=/dev/null
@REM
@SET  MAKEFP=%RESTARTDIR%\%JOB%.efp
@SET   GAMMA=%RESTARTDIR%\%JOB%.gamma
@SET TRAJECT=%RESTARTDIR%\%JOB%.trj
@SET RESTART=%RESTARTDIR%\%JOB%.rst
@SET   INPUT=%SCRATCHDIR%\%JOB%.F05
@SET   PUNCH=%RESTARTDIR%\%JOB%.dat
@SET  AOINTS=%SCRATCHDIR%\%JOB%.F08
@SET  MOINTS=%SCRATCHDIR%\%JOB%.F09
@SET DICTNRY=%SCRATCHDIR%\%JOB%.F10
@SET DRTFILE=%SCRATCHDIR%\%JOB%.F11
@SET CIVECTR=%SCRATCHDIR%\%JOB%.F12
@SET CASINTS=%SCRATCHDIR%\%JOB%.F13
@SET  CIINTS=%SCRATCHDIR%\%JOB%.F14
@SET  WORK15=%SCRATCHDIR%\%JOB%.F15
@SET  WORK16=%SCRATCHDIR%\%JOB%.F16
@SET CSFSAVE=%SCRATCHDIR%\%JOB%.F17
@SET FOCKDER=%SCRATCHDIR%\%JOB%.F18
@SET  WORK19=%SCRATCHDIR%\%JOB%.F19
@SET  DASORT=%SCRATCHDIR%\%JOB%.F20
@SET DIABDAT=%SCRATCHDIR%\%JOB%.F21
@SET DFTINTS=%SCRATCHDIR%\%JOB%.F21
@SET DFTGRID=%SCRATCHDIR%\%JOB%.F22
@SET  JKFILE=%SCRATCHDIR%\%JOB%.F23
@SET  ORDINT=%SCRATCHDIR%\%JOB%.F24
@SET  EFPIND=%SCRATCHDIR%\%JOB%.F25
@SET PCMDATA=%SCRATCHDIR%\%JOB%.F26
@SET PCMINTS=%SCRATCHDIR%\%JOB%.F27
@SET SVPWRK1=%SCRATCHDIR%\%JOB%.F26
@SET SVPWRK2=%SCRATCHDIR%\%JOB%.F27
@REM
@SET  COSCAV=%SCRATCHDIR%\%JOB%.F26
@SET COSDATA=%RESTARTDIR%\%JOB%.cosmo
@SET  COSPOT=%RESTARTDIR%\%JOB%.pot
@REM
@SET   MLTPL=%SCRATCHDIR%\%JOB%.F28
@SET  MLTPLT=%SCRATCHDIR%\%JOB%.F29
@SET  DAFL30=%SCRATCHDIR%\%JOB%.F30
@SET  SOINTX=%SCRATCHDIR%\%JOB%.F31
@SET  SOINTY=%SCRATCHDIR%\%JOB%.F32
@SET  SOINTZ=%SCRATCHDIR%\%JOB%.F33
@SET  SORESC=%SCRATCHDIR%\%JOB%.F34
@REM
@REM 35 is used by RESTART, see above
@REM
@SET GCILIST=%SCRATCHDIR%\%JOB%.F37
@SET HESSIAN=%SCRATCHDIR%\%JOB%.F38
@SET QMMMTEI=%SCRATCHDIR%\%JOB%.F39
@SET SOCCDAT=%SCRATCHDIR%\%JOB%.F40
@SET  AABB41=%SCRATCHDIR%\%JOB%.F41
@SET  BBAA42=%SCRATCHDIR%\%JOB%.F42
@SET  BBBB43=%SCRATCHDIR%\%JOB%.F43
@SET  REMD  =%SCRATCHDIR%\%JOB%.F44
@SET  MCQD50=%SCRATCHDIR%\%JOB%.F50
@SET  MCQD51=%SCRATCHDIR%\%JOB%.F51
@SET  MCQD52=%SCRATCHDIR%\%JOB%.F52
@SET  MCQD53=%SCRATCHDIR%\%JOB%.F53
@SET  MCQD54=%SCRATCHDIR%\%JOB%.F54
@SET  MCQD55=%SCRATCHDIR%\%JOB%.F55
@SET  MCQD56=%SCRATCHDIR%\%JOB%.F56
@SET  MCQD57=%SCRATCHDIR%\%JOB%.F57
@SET  MCQD58=%SCRATCHDIR%\%JOB%.F58
@SET  MCQD59=%SCRATCHDIR%\%JOB%.F59
@SET  MCQD60=%SCRATCHDIR%\%JOB%.F60
@SET  MCQD61=%SCRATCHDIR%\%JOB%.F61
@SET  MCQD62=%SCRATCHDIR%\%JOB%.F62
@SET  MCQD63=%SCRATCHDIR%\%JOB%.F63
@SET  MCQD64=%SCRATCHDIR%\%JOB%.F64
@SET NMRINT1=%SCRATCHDIR%\%JOB%.F61
@SET NMRINT2=%SCRATCHDIR%\%JOB%.F62
@SET NMRINT3=%SCRATCHDIR%\%JOB%.F63
@SET NMRINT4=%SCRATCHDIR%\%JOB%.F64
@SET NMRINT5=%SCRATCHDIR%\%JOB%.F65
@SET NMRINT6=%SCRATCHDIR%\%JOB%.F66
@SET DCPHFH2=%SCRATCHDIR%\%JOB%.F67
@SET DCPHF21=%SCRATCHDIR%\%JOB%.F68
@SET ELNUINT=%SCRATCHDIR%\%JOB%.F67
@SET NUNUINT=%SCRATCHDIR%\%JOB%.F68
@SET   GVVPT=%SCRATCHDIR%\%JOB%.F69
@SET NUMOIN =%SCRATCHDIR%\%JOB%.F69
@SET NUMOCAS=%SCRATCHDIR%\%JOB%.F70
@SET NUELMO =%SCRATCHDIR%\%JOB%.F71
@SET NUELCAS=%SCRATCHDIR%\%JOB%.F72
@REM
@REM Next files are for RI-MP2
@REM
@SET  RIVMAT=%SCRATCHDIR%\%JOB%.F51
@SET   RIT2A=%SCRATCHDIR%\%JOB%.F52
@SET   RIT3A=%SCRATCHDIR%\%JOB%.F53
@SET   RIT2B=%SCRATCHDIR%\%JOB%.F54
@SET   RIT3B=%SCRATCHDIR%\%JOB%.F55
@REM
@REM Next files are for GMCQDPT
@REM
@SET  GMCREF=%SCRATCHDIR%\%JOB%.F70
@SET  GMCO2R=%SCRATCHDIR%\%JOB%.F71
@SET  GMCROC=%SCRATCHDIR%\%JOB%.F72
@SET  GMCOOC=%SCRATCHDIR%\%JOB%.F73
@SET  GMCCC0=%SCRATCHDIR%\%JOB%.F74
@SET  GMCHMA=%SCRATCHDIR%\%JOB%.F75
@SET  GMCEI1=%SCRATCHDIR%\%JOB%.F76
@SET  GMCEI2=%SCRATCHDIR%\%JOB%.F77
@SET  GMCEOB=%SCRATCHDIR%\%JOB%.F78
@SET  GMCEDT=%SCRATCHDIR%\%JOB%.F79
@SET  GMCERF=%SCRATCHDIR%\%JOB%.F80
@SET  GMCHCR=%SCRATCHDIR%\%JOB%.F81
@SET  GMCGJK=%SCRATCHDIR%\%JOB%.F82
@SET  GMCGAI=%SCRATCHDIR%\%JOB%.F83
@SET  GMCGEO=%SCRATCHDIR%\%JOB%.F84
@SET  GMCTE1=%SCRATCHDIR%\%JOB%.F85
@SET  GMCTE2=%SCRATCHDIR%\%JOB%.F86
@SET  GMCHEF=%SCRATCHDIR%\%JOB%.F87
@SET  GMCMOL=%SCRATCHDIR%\%JOB%.F88
@SET  GMCMOS=%SCRATCHDIR%\%JOB%.F89
@SET  GMCWGT=%SCRATCHDIR%\%JOB%.F90
@SET  GMCRM2=%SCRATCHDIR%\%JOB%.F91
@SET  GMCRM1=%SCRATCHDIR%\%JOB%.F92
@SET  GMCR00=%SCRATCHDIR%\%JOB%.F93
@SET  GMCRP1=%SCRATCHDIR%\%JOB%.F94
@SET  GMCRP2=%SCRATCHDIR%\%JOB%.F95
@SET  GMCVEF=%SCRATCHDIR%\%JOB%.F96
@SET  GMCDIN=%SCRATCHDIR%\%JOB%.F97
@SET  GMC2SZ=%SCRATCHDIR%\%JOB%.F98
@SET  GMCCCS=%SCRATCHDIR%\%JOB%.F99
@REM
@REM Next files are used only during closed shell coupled cluster runs.
@REM
@SET  CCREST=%SCRATCHDIR%\%JOB%.F70
@SET  CCDIIS=%SCRATCHDIR%\%JOB%.F71
@SET  CCINTS=%SCRATCHDIR%\%JOB%.F72
@SET CCT1AMP=%SCRATCHDIR%\%JOB%.F73
@SET CCT2AMP=%SCRATCHDIR%\%JOB%.F74
@SET CCT3AMP=%SCRATCHDIR%\%JOB%.F75
@SET    CCVM=%SCRATCHDIR%\%JOB%.F76
@SET    CCVE=%SCRATCHDIR%\%JOB%.F77
@SET CCQUADS=%SCRATCHDIR%\%JOB%.F78
@SET QUADSVO=%SCRATCHDIR%\%JOB%.F79
@SET EOMSTAR=%SCRATCHDIR%\%JOB%.F80
@SET EOMVEC1=%SCRATCHDIR%\%JOB%.F81
@SET EOMVEC2=%SCRATCHDIR%\%JOB%.F82
@SET  EOMHC1=%SCRATCHDIR%\%JOB%.F83
@SET  EOMHC2=%SCRATCHDIR%\%JOB%.F84
@SET EOMHHHH=%SCRATCHDIR%\%JOB%.F85
@SET EOMPPPP=%SCRATCHDIR%\%JOB%.F86
@SET EOMRAMP=%SCRATCHDIR%\%JOB%.F87
@SET EOMRTMP=%SCRATCHDIR%\%JOB%.F88
@SET EOMDG12=%SCRATCHDIR%\%JOB%.F89
@SET    MMPP=%SCRATCHDIR%\%JOB%.F90
@SET   MMHPP=%SCRATCHDIR%\%JOB%.F91
@SET MMCIVEC=%SCRATCHDIR%\%JOB%.F92
@SET MMCIVC1=%SCRATCHDIR%\%JOB%.F93
@SET MMCIITR=%SCRATCHDIR%\%JOB%.F94
@SET  EOMVL1=%SCRATCHDIR%\%JOB%.F95
@SET  EOMVL2=%SCRATCHDIR%\%JOB%.F96
@SET EOMLVEC=%SCRATCHDIR%\%JOB%.F97
@SET  EOMHL1=%SCRATCHDIR%\%JOB%.F98
@SET  EOMHL2=%SCRATCHDIR%\%JOB%.F99
@SET  CCVVVV=%SCRATCHDIR%\%JOB%.F80
@REM
@REM Next files are used only during open shell coupled cluster runs.
@REM
@SET AMPROCC=%SCRATCHDIR%\%JOB%.F70
@SET ITOPNCC=%SCRATCHDIR%\%JOB%.F71
@SET FOCKMTX=%SCRATCHDIR%\%JOB%.F72
@SET  LAMB23=%SCRATCHDIR%\%JOB%.F73
@SET   VHHAA=%SCRATCHDIR%\%JOB%.F74
@SET   VHHBB=%SCRATCHDIR%\%JOB%.F75
@SET   VHHAB=%SCRATCHDIR%\%JOB%.F76
@SET    VMAA=%SCRATCHDIR%\%JOB%.F77
@SET    VMBB=%SCRATCHDIR%\%JOB%.F78
@SET    VMAB=%SCRATCHDIR%\%JOB%.F79
@SET    VMBA=%SCRATCHDIR%\%JOB%.F80
@SET  VHPRAA=%SCRATCHDIR%\%JOB%.F81
@SET  VHPRBB=%SCRATCHDIR%\%JOB%.F82
@SET  VHPRAB=%SCRATCHDIR%\%JOB%.F83
@SET  VHPLAA=%SCRATCHDIR%\%JOB%.F84
@SET  VHPLBB=%SCRATCHDIR%\%JOB%.F85
@SET  VHPLAB=%SCRATCHDIR%\%JOB%.F86
@SET  VHPLBA=%SCRATCHDIR%\%JOB%.F87
@SET    VEAA=%SCRATCHDIR%\%JOB%.F88
@SET    VEBB=%SCRATCHDIR%\%JOB%.F89
@SET    VEAB=%SCRATCHDIR%\%JOB%.F90
@SET    VEBA=%SCRATCHDIR%\%JOB%.F91
@SET   VPPPP=%SCRATCHDIR%\%JOB%.F92
@SET INTERM1=%SCRATCHDIR%\%JOB%.F93
@SET INTERM2=%SCRATCHDIR%\%JOB%.F94
@SET INTERM3=%SCRATCHDIR%\%JOB%.F95
@SET ITSPACE=%SCRATCHDIR%\%JOB%.F96
@SET INSTART=%SCRATCHDIR%\%JOB%.F97
@SET  ITSPC3=%SCRATCHDIR%\%JOB%.F98
@REM
@REM Next files are used only during elongation method runs.
@REM
@FINDSTR /I /R /C:"NELONG=" %SCRATCHDIR%\%JOB%.F05
@IF NOT ERRORLEVEL 1 (
    @IF [%4]==[] (
      @SET ELGNAME=ELGFILE
    ) ELSE (
      @SET ELGNAME=%4
    )
)
@IF NOT ERRORLEVEL 1 (
  @SET AOINTS=%SCRATCHDIR%\%ELGNAME%.F08
  @SET ELGDOS=%RESTARTDIR%\%JOB%.ldos
  @SET ELGDAT=%SCRATCHDIR%\%ELGNAME%.F71
  @SET ELGPAR=%SCRATCHDIR%\%ELGNAME%.F72
  @SET ELGCUT=%SCRATCHDIR%\%ELGNAME%.F74
  @SET ELGVEC=%SCRATCHDIR%\%ELGNAME%.F75
  @SET EGINTA=%SCRATCHDIR%\%ELGNAME%.F77
  @SET EGINTB=%SCRATCHDIR%\%ELGNAME%.F78
  @SET EGTDHF=%SCRATCHDIR%\%ELGNAME%.F79
  @SET EGTEST=%SCRATCHDIR%\%ELGNAME%.F80
)
@REM
@REM Next files are used only during extended TDHF package runs.
@REM
@SET  OLI201=%SCRATCHDIR%\%JOB%.F201
@SET  OLI202=%SCRATCHDIR%\%JOB%.F202
@SET  OLI203=%SCRATCHDIR%\%JOB%.F203
@SET  OLI204=%SCRATCHDIR%\%JOB%.F204
@SET  OLI205=%SCRATCHDIR%\%JOB%.F205
@SET  OLI206=%SCRATCHDIR%\%JOB%.F206
@SET  OLI207=%SCRATCHDIR%\%JOB%.F207
@SET  OLI208=%SCRATCHDIR%\%JOB%.F208
@SET  OLI209=%SCRATCHDIR%\%JOB%.F209
@SET  OLI210=%SCRATCHDIR%\%JOB%.F210
@SET  OLI211=%SCRATCHDIR%\%JOB%.F211
@SET  OLI212=%SCRATCHDIR%\%JOB%.F212
@SET  OLI213=%SCRATCHDIR%\%JOB%.F213
@SET  OLI214=%SCRATCHDIR%\%JOB%.F214
@SET  OLI215=%SCRATCHDIR%\%JOB%.F215
@SET  OLI216=%SCRATCHDIR%\%JOB%.F216
@SET  OLI217=%SCRATCHDIR%\%JOB%.F217
@SET  OLI218=%SCRATCHDIR%\%JOB%.F218
@SET  OLI219=%SCRATCHDIR%\%JOB%.F219
@SET  OLI220=%SCRATCHDIR%\%JOB%.F220
@SET  OLI221=%SCRATCHDIR%\%JOB%.F221
@SET  OLI222=%SCRATCHDIR%\%JOB%.F222
@SET  OLI223=%SCRATCHDIR%\%JOB%.F223
@SET  OLI224=%SCRATCHDIR%\%JOB%.F224
@SET  OLI225=%SCRATCHDIR%\%JOB%.F225
@SET  OLI226=%SCRATCHDIR%\%JOB%.F226
@SET  OLI227=%SCRATCHDIR%\%JOB%.F227
@SET  OLI228=%SCRATCHDIR%\%JOB%.F228
@SET  OLI229=%SCRATCHDIR%\%JOB%.F229
@SET  OLI230=%SCRATCHDIR%\%JOB%.F230
@SET  OLI231=%SCRATCHDIR%\%JOB%.F231
@SET  OLI232=%SCRATCHDIR%\%JOB%.F232
@SET  OLI233=%SCRATCHDIR%\%JOB%.F233
@SET  OLI234=%SCRATCHDIR%\%JOB%.F234
@SET  OLI235=%SCRATCHDIR%\%JOB%.F235
@SET  OLI236=%SCRATCHDIR%\%JOB%.F236
@SET  OLI237=%SCRATCHDIR%\%JOB%.F237
@SET  OLI238=%SCRATCHDIR%\%JOB%.F238
@SET  OLI239=%SCRATCHDIR%\%JOB%.F239
@REM
@REM Next files are used only during divide-and-conquer runs
@REM
@SET   DCSUB=%SCRATCHDIR%\%JOB%.F250
@SET   DCVEC=%SCRATCHDIR%\%JOB%.F251
@SET   DCEIG=%SCRATCHDIR%\%JOB%.F252
@SET    DCDM=%SCRATCHDIR%\%JOB%.F253
@SET   DCDMO=%SCRATCHDIR%\%JOB%.F254
@SET     DCQ=%SCRATCHDIR%\%JOB%.F255
@SET     DCW=%SCRATCHDIR%\%JOB%.F256
@SET   DCEDM=%SCRATCHDIR%\%JOB%.F257
@REM
@REM Next files are used only during LMO hyperpolarizability analysis
@REM
@SET LHYPWRK=%SCRATCHDIR%\%JOB%.F297
@SET LHYPWK2=%SCRATCHDIR%\%JOB%.F298
@SET BONDDPF=%SCRATCHDIR%\%JOB%.F299
@REM
@REM Next value is used only within the VB2000 add-on code
@REM
@SET VB2000PATH=%GAMESSDIR%\vb2000
@SET GMSJOBNAME=%JOB%
@REM
@REM Next files are used during EFMO runs
@REM
@FINDSTR /I /R /C:"^ IEFMO=" %SCRATCHDIR%\%JOB%.F05 > NUL
@IF NOT ERRORLEVEL 1 (
  @SET EFMOI=%SCRATCHDIR%\$JOB.F102
  @SET EFMOF=%SCRATCHDIR%\$JOB.F103
)
@REM
@REM Next files are used only during CIM runs
@REM
@SET CIMFILE=%RESTARTDIR%\%JOB%.cim
@SET CIMDMN=%RESTARTDIR%\%JOB%.dmn
@SET CIMDCT=%SCRATCHDIR%\%JOB%.Fcdt
@SET CIMAOI=%SCRATCHDIR%\%JOB%.Fcao
@SET CIMMOI=%SCRATCHDIR%\%JOB%.Fcmo
@REM
@REM "Group DDI" will be set to true, below, if=$GDDI input is found
@REM
@SET GDDIJOB=FALSE
@REM
@FINDSTR /I /R /C:"^ \$GDDI" %SCRATCHDIR%\%JOB%.F05 > NUL
@IF NOT ERRORLEVEL 1 (
  @SET GDDIJOB=TRUE
  @REM
  @ECHO.
  @ECHO This is a GDDI run, keeping output files on local disks
  @ECHO until the very end of the run, when they'll be saved from
  @ECHO the master process in the first group, only.
  @REM
  @SET  MAKEFP=%SCRATCHDIR%\%JOB%.F01
  @SET TRAJECT=%SCRATCHDIR%\%JOB%.F04
  @SET  OUTPUT=%SCRATCHDIR%\%JOB%.F06
  @SET   PUNCH=%SCRATCHDIR%\%JOB%.F07
  @SET RESTART=%SCRATCHDIR%\%JOB%.F35
  @REM
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
  @CALL %GAMESSDIR%\WINDOWS\GMSENV.BAT
) ELSE (
  @REM
  @REM Maybe the DEBUG flag is on.
  @REM
  @IF %DEBUG%==TRUE (
    @ECHO.
    @CALL %GAMESSDIR%\WINDOWS\GMSENV.BAT
  )
)
@REM
@REM Data left over from a previous run might be precious, stop if found.
@REM
@ECHO.
@ECHO =========================================================================
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
@REM We don't need to have a HOSTFILE because running locally.
@REM
@SET PROCFILE=%SCRATCHDIR%\%JOB%.processes.mpd
@IF EXIST %PROCFILE% DEL /Q %PROCFILE%
@REM
@REM Allow for compute process and data server (one pair per core).
@REM
@SET /A NPROCS=NCPUS+NCPUS
@REM
@ECHO -env ENVFIL %ENVFIL% -n %NPROCS% %GAMESSDIR%\gamess.%VERSION%.exe >> %PROCFILE%
@REM
@IF %DEBUG%==TRUE (
  @ECHO [Contents of procfile]
  @ECHO %PROCFILE%
  @ECHO.
  @TYPE %PROCFILE%
  @ECHO =========================================================================
)
@REM
@REM  At this point we are ready to run!
@REM
@ECHO Microsoft MPI ^(MS-MPI^) will be running GAMESS on 1 node^(s^).
@ECHO The binary will be kicked off by 'mpiexec' is gamess.%VERSION%.exe
@ECHO MS-MPI will run %NCPUS% compute process^(es^) and %NCPUS% data server^(s^).
@ECHO.
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
@REM
@REM We stall for 5 seconds. Hopefully this is enough time for
@REM mpiexec to release hold on %LOGFILE%
@REM
@PING -n 5 127.0.0.1 > NUL
@REM
@REM Tidy up some GDDI stuff
@REM
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
      @COPY %SCRATCHDIR%\%JOB%.F06.* %RESTARTDIR% > NUL
      @COPY %SCRATCHDIR%\%JOB%.F07 %RESTARTDIR%\%JOB%.dat > NUL
  @IF EXIST %SCRATCHDIR%\%JOB%.F04 @COPY %SCRATCHDIR%\%JOB%.F04 %RESTARTDIR%\%JOB%.trj > NUL
  @IF EXIST %SCRATCHDIR%\%JOB%.F35 @COPY %SCRATCHDIR%\%JOB%.F35 %RESTARTDIR%\%JOB%.rst > NUL
)
@REM
@IF %REDIRECTION%==TRUE (
  @ECHO            Disk usage ^(in bytes^)       >>%LOGFILE%
  @ECHO  ========================================>>%LOGFILE%
  @DIR %SCRATCHDIR%\%JOB%.* | FIND "/"      >>%LOGFILE%
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
  @DIR %SCRATCHDIR%\%JOB%.* | FIND "/"
  @REM
  @ECHO  ----------------------------------------
  @ECHO             Local date and time
  @ECHO  ========================================
  @ECHO %DATE% - %TIME%
)
@REM
@REM Time to clean up.
@REM
@REM
@REM Remove scratch files
@REM
@DEL /Q %SCRATCHDIR%\%JOB%.F*
@REM
@REM Remove temporary file used by this batch to store variables.
@REM
@IF %CLEANALLLOCAL%==TRUE (
  @DEL /Q %ENVFIL%
  @DEL /Q %PROCFILE%
)
@IF EXIST .batch.temp @DEL /Q .batch.temp
@REM
@REM FIN
@REM
@EXIT /B
