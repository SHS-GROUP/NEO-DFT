#!/bin/csh
#
#    Set up input file(s) for replica-exchange molecular dynamics.
#    This script is called from 'rungms' during batch execution,
#    for runs with MREMD=1 in $MD inputs groups.      8/2013
#
#    a) fresh starts use a single input file, which must be given to
#       every node to ensure the master process in each subgroup can
#       find it.   This is like a normal GDDI run.
#    b) restarts give a different input file to each subgroup, with
#       every node in the 1st, 2nd, 3rd... subgroups receiving the
#       same input file.  2nd group gets .001, 3rd group .002, and so
#       on, with the primary group receiving an input with no suffix.
#
#   This script uses a number of environment variables set by 'rungms',
#   and needs two ordinary paramters passed during its invocation.
#   There should be no need to customize anything in this script.
#
set TARGET=$1
set nREMDreplica=$2
#
echo " "
echo This is a REMD job, using $nREMDreplica replicas.
#
if ($NNODES < $nREMDreplica) then
   echo The number of compute nodes is smaller than that of replica,
   echo current remd.csh script does not support this situation.
   exit 8
endif
#
#      first we take care of case "b"
#
#      Logic here is not fully general: instead we assume the likely
#      case of an equal number of nodes -Mgroup- in every subgroup.
#      The script does trap such an error, however.
#
if (-e $JOB.inp.001) then
   echo "Your job looks like a REMD restart, using multiple input files."
   echo "REMD restart files should be in `pwd` and be named"
   echo "    $JOB.inp, $JOB.inp.001, $JOB.inp.002, ..."
   echo "with the first processor group's file NOT containing a .000 suffix."
   echo "REMD will now copy your input files to every node..."

   @ Mgroup = $NNODES / $nREMDreplica
   
#       every group must contain the same number of nodes: Mgroup
   @ nnodeschk = $Mgroup * $nREMDreplica
   if ($nnodeschk != $NNODES) then
      echo error, this job is using $NNODES nodes, and $nREMDreplica replicas,
      echo but the latter does not evenly divide into the former.
      echo current remd.csh script does not support this situation.
      exit 8
   endif

#            give first input (not numbered .000!) to the first group of nodes
#            'rungms' already gave the input to the 1st node of 1st group.
   @ m=2
   while ($m <= $Mgroup)
      switch ($TARGET)
         case sockets:
            host=$HOSTLIST[$m]
            host=`echo $host | cut -f 1 -d :` # drop anything behind a colon
            breaksw
         case mpi:
            set host=`sed -n -e "$m p" $HOSTFILE`
            set host=$host[1]
            breaksw
         default:
            echo "unknown execution target in remd.csh"
            exit 4
            breaksw
      endsw
      echo $DDI_RCP $JOB.inp ${host}:$SCR/$JOB.F05
           $DDI_RCP $JOB.inp ${host}:$SCR/$JOB.F05
      @ m++
   end
#            give .001, .002, ... inputs to all nodes in each one's subgroup
   @ mrep=2
   @ mversion=1
   while ($mrep <= $nREMDreplica)
      @ m = ($mrep - 1) * $Mgroup + 1
      @ Nnmul=1
      while ($Nnmul <= $Mgroup)
         switch ($TARGET)
            case sockets:
               host=$HOSTLIST[$m]
               host=`echo $host | cut -f 1 -d :` # drop anything behind a colon
               breaksw
            case mpi:
               set host=`sed -n -e "$m p" $HOSTFILE`
               set host=$host[1]
               breaksw
            default:
               echo "unknown execution target in remd.csh"
               exit 4
               breaksw
         endsw
                             set suffix=00$mversion
         if ($mversion >  9) set suffix=0$mversion
         if ($mversion > 99) set suffix=$mversion
         echo $DDI_RCP $JOB.inp.$suffix ${host}:$SCR/$JOB.F05
              $DDI_RCP $JOB.inp.$suffix ${host}:$SCR/$JOB.F05
         @ m++
         @ Nnmul++
      end
      @ mrep++
      @ mversion++
   end
#
#      while this takes care of the simpler case "a"...
#
else
   echo 'Your job looks like a fresh new REMD run, using one input file.'
   echo 'REMD will now copy your input file to every node...'
   @ m=2   # input has already been copied into the master node, by 'rungms'
   while ($m <= $NNODES)
      switch ($TARGET)
         case sockets:
            host=$HOSTLIST[$m]
            host=`echo $host | cut -f 1 -d :` # drop anything behind a colon
            breaksw
         case mpi:
            set host=`sed -n -e "$m p" $HOSTFILE`
            set host=$host[1]
            breaksw
         default:
            echo "unknown execution target in remd.csh"
            exit 4
            breaksw
      endsw
      echo $DDI_RCP $SCR/$JOB.F05 ${host}:$SCR/$JOB.F05
           $DDI_RCP $SCR/$JOB.F05 ${host}:$SCR/$JOB.F05
      @ m++
   end
endif
#
echo " "
exit 0
