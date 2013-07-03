@REM
@REM
@REM May   28, 2013 :: Sarom Leang :: New exam files
@REM
@REM
@REM March 25, 2010 :: Sarom Sok :: Windows batch file version of the runall script.
@REM
@SETLOCAL
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
  @ECHO       SUPRESS=[TRUE,FALSE] Option to 'not show'/'show' the job's environmental
  @ECHO                            variable list.
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
@REM
@ECHO -------------------------------------------------------------------------
@ECHO.
@ECHO Usage: runall.bat [version] [ncpus]
@ECHO  [version] = The GAMESS version number                      ^(default: 00^)
@ECHO  [ncpus]   = The number of CPUs requested for this job      ^(default:  1^)
@ECHO              [ncpus] = 1 will run all the exam files in serial
@ECHO              [ncpus] ^> 1 will run all the exam files in parallel except
@ECHO              for 05, 23, 25, 27, 32, 39, 42, and 45-47 which only run in serial
@ECHO =========================================================================
@REM
@IF [%1]==[] (
  @SET VERSION=00
) ELSE (
  @SET VERSION=%1
)
@REM
@IF [%2]==[] (
  @SET NCPUS=1
) ELSE (
  @SET NCPUS=%2
)
@REM
@SET JOB=exam
@REM
@IF NOT EXIST %GAMESSDIR%\gamess.%VERSION%.exe (
  @ECHO -------------------------------------------------------------------------
  @ECHO "Oh no you didn't!"
  @ECHO ERROR :: The binary gamess.%VERSION%.exe does not exist
  @ECHO =========================================================================
  @ECHO.
  @ECHO Now exiting.
  @EXIT /B
)
@REM
@IF %ncpus%==1 (
  @CALL :runallserial
) ELSE (
  @CALL :runallparallel
)
@CALL :finish
@EXIT /B

:runallserial
@REM
@ECHO ----------------------------------------------
@ECHO   Running GAMESS exam files 1 - 47 in serial
@ECHO   Using gamess.%VERSION%.exe binary
@ECHO   NCPUS = %NCPUS%
@ECHO ----------------------------------------------
@SET  CURRENTJOB=%JOB%01
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%02
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%03
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%04
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%05
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%06
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%07
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%08
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%09
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%10
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%11
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%12
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%13
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%14
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%15
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%16
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%17
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%18
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%19
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%20
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%21
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%22
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%23
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%24
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%25
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%26
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%27
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%28
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%29
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%30
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%31
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%32
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%33
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%34
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%35
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%36
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%37
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%38
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%39
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%40
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%41
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%42
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%43
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%44
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%45
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%46
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%47
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 5 127.0.0.1 > NUL
@REM
@GOTO:eof

:runallparallel
@REM
@ECHO ----------------------------------------------
@ECHO  Running GAMESS exam files 1 - 47 in parallel
@ECHO.
@ECHO   The following serial-only jobs are omitted
@ECHO   05, 23, 25, 27, 32, 39, 42, 45, 46, and 47
@ECHO.
@ECHO   Using gamess.%VERSION%.exe binary
@ECHO   NCPUS = %NCPUS%
@ECHO ----------------------------------------------
@SET  CURRENTJOB=%JOB%01
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%02
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%03
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%04
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%06
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%07
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%08
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%09
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%10
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%11
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%12
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%13
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%14
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%15
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%16
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%17
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%18
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%19
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%20
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%21
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%22
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%24
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%26
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%28
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%29
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%30
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%31
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%33
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%34
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%35
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%36
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%37
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%38
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%40
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%41
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%43
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@SET  CURRENTJOB=%JOB%44
@ECHO   Running %CURRENTJOB%
@CALL rungms %CURRENTJOB% %VERSION% %NCPUS% 0 %CURRENTJOB%.log > NUL
@PING -n 10 127.0.0.1 > NUL
@REM
@GOTO:eof

:finish
@ECHO ----------------------------------------------
@ECHO -                    FIN                     -
@ECHO ----------------------------------------------
@REM
@GOTO:eof