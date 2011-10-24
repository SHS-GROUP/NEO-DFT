setlocal
for /F "usebackq" %%i in (`hostname`) do set HOSTNAME=%%i
echo %HOSTNAME%
@REM So we know where we are when looking at the HPC Cluster Manager
cd %LOCALDIR%
if exist %JOBNAME%.F06.* copy %JOBNAME%.F06.* %WORKDIR%
if exist %JOBNAME%.F07   copy %JOBNAME%.F07 %WORKDIR%\%JOBNAME%.dat
if exist %JOBNAME%.F04   copy %JOBNAME%.F04 %WORKDIR%\%JOBNAME%.trj
if exist %JOBNAME%.F35   copy %JOBNAME%.F35 %WORKDIR%\%JOBNAME%.rst
cd ..
rmdir /S /Q %LOCALDIR%
