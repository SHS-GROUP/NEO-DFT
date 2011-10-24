setlocal
for /F "usebackq" %%i in (`hostname`) do set HOSTNAME=%%i
echo %HOSTNAME%
@REM So we know where we are when looking at the HPC Cluster Manager
@REM Take care of possible GDDI runs
if exist %JOBNAME%.log.master echo  ----------------------------------------------------- >> %JOBNAME%.log.master
if exist %JOBNAME%.log.master echo             Merging GDDI Timing Information            >> %JOBNAME%.log.master
if exist %JOBNAME%.log.master echo        Please ignore any group information below       >> %JOBNAME%.log.master
if exist %JOBNAME%.log.master echo  ===================================================== >> %JOBNAME%.log.master
if exist %JOBNAME%.log.master type  %JOBNAME%.log >> %JOBNAME%.log.master
if exist %JOBNAME%.log.master copy /Y %JOBNAME%.log.master  %JOBNAME%.log
if exist %JOBNAME%.log.master erase %JOBNAME%.log.master
@REM Clean up user's local scratch directory
erase %WORKDIR%\%JOBNAME%.%JOBID%.*
