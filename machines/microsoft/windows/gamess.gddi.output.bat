setlocal
for /F "usebackq" %%i in (`hostname`) do set HOSTNAME=%%i
echo %HOSTNAME%
@REM So we know where we are when looking at the HPC Cluster Manager
cd %LOCALDIR%
if exist %JOBNAME%.F06  type %JOBNAME%.F06 > %LOGFILE%.master
