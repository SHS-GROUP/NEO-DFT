#!/bin/bash
#   
#   for execution of GAMESS on the IBM BG/P model
#   
#   Assuming this script is called 'gmsbgp', type-
#    >   gmsbgp   [job]   [nodes]   [minutes]
#   
#   Things to note:
#   1) Single (colon-separated multicomponent) argument in the bash syntax
#   after the '--env' option replaces the 'FILE_GETENV' stuff needed on L.
#   2) names for all the files used by GAMESS are not given below.  Check
#   the standard 'rungms' for a complete listing of all possible files.
#   3) 'ProjectGoesHere' is the name of the project I bill time to.
#   4) IMPORTANT: It is not necessary to write files  
#   OUTPUT,PUNCH,DICTNRY,SOCCDAT to /tmp, this is specific to exploiting  
#   the 'RAMDISK' feature of BG/P (but this critical to performance for  
#   FMO) - this part is a work-in-progress.
#   5) You can use --proccount, but you still need -n and --mode, so it is  
#   largely useless. BG needs the 'mode' so it actually computes the proc  
#   count by ('nodes')x('mode'), ignoring --proccount.
#
cp $1.inp JOB.F05
#
qsub -q prod -t $3 -n $2 --mode vn -A ProjectGoesHere --env \
OUTPUT=/tmp/JOB.F06:\
ERICFMT=/intrepid-fs0/users/fletcher/persistent/gamess/ericfmt.dat:\
IRCDATA=JOB.irc:\
INPUT=JOB.F05:\
PUNCH=/tmp/JOB.dat:\
AOINTS=JOB.F08:\
MOINTS=JOB.F09:\
DICTNRY=/tmp/JOB.F10:\
DRTFILE=JOB.F11:\
CIVECTR=JOB.F12:\
CASINTS=JOB.F13:\
CIINTS=JOB.F14:\
WORK15=JOB.F15:\
WORK16=JOB.F16:\
CSFSAVE=JOB.F17:\
FOCKDER=JOB.F18:\
WORK19=JOB.F19:\
DASORT=JOB.F20:\
DFTINTS=JOB.F21:\
DFTGRID=JOB.F22:\
JKFILE=JOB.F23:\
ORDINT=JOB.F24:\
EFPIND=JOB.F25:\
PCMDATA=JOB.F26:\
PCMINTS=JOB.F27:\
SVPWRK1=JOB.F26:\
SVPWRK2=JOB.F27:\
MLTPL=JOB.F28:\
MLTPLT=JOB.F29:\
DAFL30=JOB.F30:\
SOINTX=JOB.F31:\
SOINTY=JOB.F32:\
SOINTZ=JOB.F33:\
SORESC=JOB.F34:\
SIMEN=JOB.simen:\
SIMCOR=JOB.simcor:\
GCILIST=JOB.F37:\
HESSIAN=JOB.F38:\
SOCCDAT=/tmp/JOB.F40:\
AABB41=JOB.F41:\
BBAA42=JOB.F42:\
BBBB43=JOB.F43:\
MCQD50=JOB.F50:\
MCQD51=JOB.F51:\
MCQD52=JOB.F52:\
MCQD53=JOB.F53:\
MCQD54=JOB.F54:\
MCQD55=JOB.F55:\
MCQD56=JOB.F56:\
MCQD57=JOB.F57:\
MCQD58=JOB.F58:\
MCQD59=JOB.F59:\
MCQD60=JOB.F60:\
MCQD61=JOB.F61:\
MCQD62=JOB.F62:\
MCQD63=JOB.F63:\
MCQD64=JOB.F64:\
NMRINT1=JOB.F61:\
NMRINT2=JOB.F62:\
NMRINT3=JOB.F63:\
NMRINT4=JOB.F64:\
NMRINT5=JOB.F65:\
NMRINT6=JOB.F66:\
DCPHFH2=JOB.F67:\
DCPHF21=JOB.F68:\
BG_MAPPING=XYZT   /intrepid-fs0/users/fletcher/persistent/gamess/gamess.fmo.x
