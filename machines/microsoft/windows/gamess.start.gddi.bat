setlocal
for /F "usebackq" %%i in (`hostname`) do set HOSTNAME=%%i
echo %HOSTNAME%
@REM So we know where we are when looking at the HPC Cluster Manager
mkdir %LOCALDIR%
cd %LOCALDIR%
copy %WORKDIR%\%JOBNAME%.%JOBID%.F05 %LOCALDIR%\%JOBNAME%.F05
