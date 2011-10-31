setlocal
for /F "usebackq" %%i in (`hostname`) do set HOSTNAME=%%i
echo %HOSTNAME%
@REM So we know where we are when looking at the HPC Cluster Manager
rmdir /S /Q %LOCALDIR%
