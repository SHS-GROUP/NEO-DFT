setlocal
for /F "usebackq" %%i in (`hostname`) do set HOSTNAME=%%i
echo %HOSTNAME%
@REM So we know where we are when looking at the HPC Cluster Manager
cd %LOCALDIR%
echo  Files used on %HOSTNAME%                 >      %JOBNAME%.%JOBID%.%HOSTNAME%.diskusage
echo  ---------------------------------------->>      %JOBNAME%.%JOBID%.%HOSTNAME%.diskusage
dir %JOBNAME%.F* | FIND "/"                   >>      %JOBNAME%.%JOBID%.%HOSTNAME%.diskusage
echo  ---------------------------------------->>      %JOBNAME%.%JOBID%.%HOSTNAME%.diskusage
copy %JOBNAME%.%JOBID%.%HOSTNAME%.diskusage %WORKDIR%\%JOBNAME%.%JOBID%.%HOSTNAME%.diskusage
