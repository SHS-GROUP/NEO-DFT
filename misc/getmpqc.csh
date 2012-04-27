#!/bin/csh
# This script is designed to download then install libint and mpqc. This script is set to be run in $GMS_BUILD_DIR/misc.
# Once GAMESS has been istalled, the environment variables from install.info are used to gather information about the system math library. If this script does not work as is, the math libraries may have specified below.

chdir ..

if (-e install.info) then
   source install.info
else
   echo "Please run 'config' first, to set up GAMESS compiling information"
   exit 4
endif

set make_j=1
set MATH_LIBS=' '
set MATH_INC_PATH = ' '
set MATH_LIB_PATH = ' '

switch ($GMS_MATHLIB)

  case mkl:
    switch ($GMS_MKL_VERNO)
     case 10:
     case 11:
     case 12:
      set MATH_LIBS='-lmkl_intel_ilp64 -lmkl_sequential -lmkl_core -lm'
      set MATH_LIB_PATH = `dirname $GMS_MATHLIB_PATH`
      set MATH_LIB_PATH = `dirname $MATH_LIB_PATH`
      set MATH_INC_PATH="-I$MATH_LIB_PATH/include"
      set MATH_LIB_PATH="-L$GMS_MATHLIB_PATH"
     breaksw
    endsw
  breaksw

  case atlas:
   set MATH_LIB_PATH = "-L$GMS_MATHLIB_PATH"
   set MATH_LIBS='-lf77blas -latlas'
  breaksw

  case acml:
   set MATH_LIB_PATH = "$GMS_MATHLIB_PATH"
   set MATH_LIBS="$MATH_LIB_PATH/libacml.a -lm"
   set MATH_INC_PATH="-I$MATH_LIB_PATH/include"
  breaksw

  case none:
   echo 'ABORT: sytem math library required'
   exit
  breaksw

endsw

# MPQC is hardwired for 64 bit integers...
set GAMESS64 = '--enable-blas-f77-i8'

####
# pull libint
####
rm -rf libint
curl -L http://sourceforge.net/projects/libint/files/libint-for-gamess%2Bmpqc/libint-2.0.0-stable-gamess%2Bmpqc.tgz/download | tar -xvz
mv libint-2.0.0-stable-gamess+mpqc libint

####
# pull mpqc
####
rm -rf mpqc
curl -L http://mpqc.hg.sourceforge.net/hgweb/mpqc/mpqc/archive/tip.tar.gz | tar -xvzmv
mv mpqc-* mpqc

####
# unpack and compile libint
####
chdir libint
./configure CXX=g++
make -j$make_j
# N.B. this allows libint to be usable without make install
chdir include; ln -s . libint2; chdir ..

####
# unpack and compile mpqc
####
chdir ../mpqc
make configure
# WARNING: need to feed proper fortran compiler and BLAS/LAPACK libraries to MPQC
#./configure CXX=g++ --with-f77=ifort --with-libint2=`pwd`/../libint --with-include="-I`pwd`/../libint/include $MATH_INC_PATH" --with-integrallibint2-normconv=exact $GAMESS64 --with-libs="$MATH_LIBS" 

./configure --with-libint2=`pwd`/../libint --with-include="-I`pwd`/../libint/include $MATH_INC_PATH" --with-integrallibint2-normconv=exact F77=$GMS_FORTRAN CXX=g++ --with-libs="$MATH_LIBS" --with-libdirs="$MATH_LIB_PATH"  $GAMESS64

make -j$make_j
chdir src/bin/pt2r12
make gamesslib

echo DONE

