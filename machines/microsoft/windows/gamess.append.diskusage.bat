setlocal
for /F "usebackq" %%i in (`hostname`) do set HOSTNAME=%%i
echo %HOSTNAME%
@REM So we know where we are when looking at the HPC Cluster Manager
pushd %WORKDIR%
echo  Disk usage information (in bytes)       >>%LOGFILE%
echo  ---------------------------------------->>%LOGFILE%
forfiles -m %JOBNAME%.%JOBID%.*.diskusage -c "CMD /C TYPE @FILE >> %LOGFILE%"
echo %DATE%>>%LOGFILE%
echo %TIME%>>%LOGFILE%
popd
