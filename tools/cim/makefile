#------------------------------------------------
# makefile for ompjob and mpijob of CIM programs
#------------------------------------------------
#
  comp = ifort
  ompopt  = -openmp -O2 -ftz -auto -assume byterecl -vec-report0 -w95 -cm
#
# mpicomp = mpif90
# mpiopt  = -i8 -O2 -ftz -auto -assume byterecl -vec-report0 -w95 -cm
#
# comp = gfortran
# ompopt = ???
#
  mpicomp = mpif90
  mpiopt = -fdefault-integer-8 -O2 -std=legacy
#
#---------------------------------------
  srcomp = ompjob.for
  objomp = ompjob.o
  exeomp = ompjob
#
  srcmpi = mpijob.for
  objmpi = mpijob.o
  exempi = mpijob
#---------------------------------------
#
#-----------------------------------------------------------
exelist: ${exeomp} ${exempi}
#-----------------------------------------------------------
#
${exeomp}: ${objomp}
	${comp} ${ompopt} ${objomp} -o $@
#
${exempi}: ${objmpi}
	${mpicomp} ${mpiopt} ${objmpi} -o $@
#
#
${objomp}: ${srcomp}
	${comp} -c ${ompopt} ${srcomp} -o $@
#
${objmpi}: ${srcmpi}
	${mpicomp} -c ${mpiopt} ${srcmpi} -o $@
#
#
clean:
	rm -f *.o
tag:
	ctags -R ${srcomp} ${srcmpi}
# END
