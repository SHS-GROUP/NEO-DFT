#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(libcchem, 0.1, asadchev@gmail.com)
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([-Wall foreign])
AC_CONFIG_SRCDIR([src/basis/molecule.hpp])
AC_CONFIG_HEADER([config/config.h])

LT_INIT([disable-shared])
AC_CONFIG_MACRO_DIR([config/m4])

# store current user given compiler flags to avoid default setup via AC_PROG_CXX
OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS

# check for C++ preprocessor and compiler and the library compiler
# (might change the compiler flags)
AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP

# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS


# Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([stdint.h stdlib.h sys/time.h])


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Checks for library functions.
#AC_FUNC_ERROR_AT_LINE
#AC_FUNC_MALLOC
#AC_CHECK_FUNCS([gettimeofday sqrt])


# integer size
integer8=""
AC_CHECK_SIZEOF(void*)
if test "x$ac_cv_sizeof_voidp" = "x8"; then
   integer8="yes"
fi


AC_OPENMP
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"

# boost
sinclude(config/ax_boost_base.m4)
AX_BOOST_BASE([])
AM_CPPFLAGS="$AM_CPPFLAGS $BOOST_CPPFLAGS"

dnl # Cuda
sinclude(config/ax_cuda.m4)
AX_CUDA_ROOT()
AX_CUDA()

AM_CPPFLAGS="$AM_CPPFLAGS $CUDA_CPPFLAGS"
dnl echo "$AM_CPPFLAGS"

#      in actuality, the compiler options which are actually used
#      seem to come from the ~/libcchem/rys/configure.ac file
AC_ARG_WITH(gpu,
	AC_HELP_STRING([--with-gpu], [Select GPU: tesla, fermi, or kepler]),
	[with_gpu=$withval], [with_gpu="tesla"])
if test "x$with_gpu" = "xtesla"; then
    CUDAFLAGS="-arch=sm_13 -Xptxas -dlcm=cg $CUDAFLAGS"
elif test "x$with_gpu" = "xfermi"; then
    CUDAFLAGS="-arch=sm_20 -Xptxas -dlcm=cg -maxrregcount=80 $CUDAFLAGS"
elif test "x$with_gpu" = "xkepler"; then
    CUDAFLAGS="-arch=sm_35 -Xptxas -dlcm=cg -maxrregcount=80 $CUDAFLAGS"
else
    echo "*** unknown GPU:" $with_gpu 
fi
CUDAFLAGS="--ptxas-options=-v $CUDAFLAGS"


# VTK
sinclude(config/am_path_vtk.m4)
AM_OPTIONS_VTK
AM_PATH_VTK
AM_CONDITIONAL(HAVE_VTK, test "$with_vtk" = "yes")
if test "$with_vtk" != "no"; then
   AC_DEFINE([HAVE_VTK], [1], [with VTK])   
fi

# DOXYGEN SUPPORT
sinclude(config/ax_prog_doxygen.m4)

# doxygen features
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(OFF)
DX_PS_FEATURE(OFF)

# initialize doxygen
DX_INIT_DOXYGEN($PACKAGE_NAME, [doc/Doxyfile], [doc])

AC_CONFIG_FILES([Makefile
		src/Makefile
		src/basis/Makefile
		src/basis/config.hpp
		src/integrals/Makefile
		src/core/Makefile

		src/hf/Makefile
		src/mp2/Makefile
		src/cc/Makefile

		src/utility/Makefile
		src/adapter/Makefile
		src/matrix/Makefile
		src/array/Makefile
		src/bindings/Makefile
		examples/Makefile])

# debug compilation
AC_ARG_ENABLE(debug,
    AC_HELP_STRING(--enable-debug, [Debug compilation (Default = no)]),
    enable_debug=$enableval, enable_debug=no)
# set compiler flags
if test "$enable_debug" = "yes" ; then
    CFLAGS="$CFLAGS -O0 -g -Wall -Wno-uninitialized"
    CXXFLAGS="$CXXFLAGS -O0 -g -Wall -Wno-uninitialized"
else
    CFLAGS="$CFLAGS -O3"
    CXXFLAGS="$CXXFLAGS -O3"
    CPPFLAGS="$CPPFLAGS -DNDEBUG"
fi


# instrumentation
AC_ARG_ENABLE(instrument,
	AC_HELP_STRING(--enable-instrument, [Enable source instrumentation]),
	enable_instrument=$enableval, enable_instrument=no)
# set compiler flags
if test "$enable_instrument" = "yes"; then
   AM_CPPFLAGS="$AM_CPPFLAGS -DBOOST_PROFILE_ENABLE"
fi


dnl # Rysq settings
dnl AC_ARG_WITH(rysq, AC_HELP_STRING([--with-rysq=PATH], [librysq installation]),
dnl             [with_rysq=$withval], [with_rysq="yes"])
dnl if test "$with_rysq" != "no"; then
dnl    if test "$with_rysq" != "yes"; then
dnl       AM_CPPFLAGS="$AM_CPPFLAGS -I$with_rysq/include"
dnl       LDFLAGS="$LDFLAGS -L$with_rysq/lib"
dnl    fi
dnl #   LIBS="$LIBS -lrysq"
dnl fi  


# hdf5
#AM_LDFLAGS="$AM_LDFLAGS -lhdf5"

## enable global arrays
AC_ARG_WITH(ga, AC_HELP_STRING([--with-ga=PATH], [Global Arrays directory]),
            [with_ga=$withval], [with_ga="no"])
#AM_CONDITIONAL(HAVE_GA, [test "$with_ga" != "no"])
if test "$with_ga" != "no"; then
   if test "$with_ga" != "yes"; then
       AM_CPPFLAGS="$AM_CPPFLAGS -I$with_ga/include"
       AM_LDFLAGS="$AM_LDFLAGS -L$with_ga/lib"
   fi
   AM_LDFLAGS="$AM_LDFLAGS -lga"
   AC_DEFINE(HAVE_GA, 1, [enable ga])
fi

## with mpi
AC_ARG_WITH(mpi, AC_HELP_STRING([--with-mpi], [Global Arrays with MPI]),
            [with_mpi=$withval], [with_mpi="no"])
if test "$with_mpi" != "no"; then
   if test "$with_mpi" != "yes"; then
       AM_CPPFLAGS="$AM_CPPFLAGS -I$with_mpi/include"
   fi
   AC_DEFINE(HAVE_MPI, 1, [enable mpi])
fi


# hdf5, must be preceded by mpi in case of parallel hdf5
sinclude([config/h5cpp.m4])

# MKL
AC_ARG_WITH(blas, AC_HELP_STRING([--with-blas],
		   [Use BLAS, --with-blas=mkl selects mkl blas]),
		 [with_blas=$withval], [with_blas="no"])
#AM_CONDITIONAL(HAVE_MKL, [test "$with_mkl" != "no"])
if test "$with_blas" != "no"; then
   AC_DEFINE(HAVE_BLAS, 1, [with BLAS])
fi
if test "$with_blas" == "mkl"; then
   AC_DEFINE(HAVE_MKL, 1, [with MKL])
fi

# Cublas
AC_ARG_WITH(cublas,
	AC_HELP_STRING([--with-cublas], [Use CUBLAS]), [with_cublas=$withval])
if test "$with_cublas" = ""; then
    if test "$with_cuda" = "no"; then
        with_cublas="no"
    else
        with_cublas="yes"
    fi
fi
AM_CONDITIONAL(HAVE_CUBLAS, [test "$with_cublas" != "no"])
if test "$with_cublas" != "no"; then
    if test "$with_cublas" != "yes"; then
        AM_CPPFLAGS="$AM_CPPFLAGS -I$with_cublas/include"
        #LDFLAGS="$LDFLAGS -L$with_cublas/lib64 -L$with_cublas/lib"
        #LDFLAGS="$LDFLAGS -lcublas -lcudart"
    fi
    AC_DEFINE(HAVE_CUBLAS, 1, [with CUBLAS])
    #m4_warn([all], ["CUBLAS=$with_cublas"])
    #echo "CUBLAS=$with_cublas"
fi


# gamess bindings
AC_ARG_WITH(gamess,
	AS_HELP_STRING([--with-gamess],[Adds GAMESS support to libcchem]),
	[if test "$withval" = "yes"; then
	    AC_DEFINE([LIBCCHEM_WITH_GAMESS], [1], [Compile with GAMESS interface])
	 fi],
        [with_gamess="no"])
AM_CONDITIONAL(LIBCCHEM_WITH_GAMESS, test "$with_gamess" = "yes")

# integer8 Fortran bindings
AC_ARG_WITH(integer8,
	AS_HELP_STRING([--with-integer8],[Assume 8-byte Fortran integer]),
	[if test "x$withval" = "xyes"; then
             integer8="yes";
         fi])

if test "x$integer8" = "xyes"; then
   AC_DEFINE([LIBCCHEM_WITH_INTEGER8], [1], [8-byte Fortran integer])
fi
		

dnl AC_ARG_ENABLE(lmax, AS_HELP_STRING([--enable-lmax[=lmax]],[Max angular momentum]), [RYSQ_LMAX=$enableval],[RYSQ_LMAX=3])
dnl AC_SUBST([RYSQ_LMAX])

dnl AC_ARG_ENABLE(sp, AS_HELP_STRING([--disable-sp],[disables SP orbitals]),[
dnl if test x$enableval = xno; then
dnl    RYSQ_SP=0
dnl else
dnl    RYSQ_SP=1
dnl fi],[RYSQ_SP=1])
dnl AC_SUBST([RYSQ_SP])

dnl AC_ARG_ENABLE(unroll-bra, AS_HELP_STRING([--enable-unroll-bra[=unroll_size]],[Sets the size of Unrolled Bra integral loop]),
dnl 			  [RYSQ_UNROLL_BRA=$enableval],[RYSQ_UNROLL_BRA=10])
dnl AC_SUBST([RYSQ_UNROLL_BRA])

dnl AC_ARG_ENABLE(cuda, AS_HELP_STRING([--enable-cuda],[Enables CUDA code to be compiled]), [
dnl if test x$enableval = xyes; then
dnl    CURYSQ=1
dnl else
dnl 	CURSYQ=0
dnl fi],[CURYSQ=0])
dnl AM_CONDITIONAL(CURYSQ, test x$CURYSQ = x1)

dnl AC_ARG_ENABLE(cu-lmax, AS_HELP_STRING([--enable-cu-lmax[=lmax]],[Max angular momentum for CUDA]), [CURYSQ_LMAX=$enableval],[CURYSQ_LMAX=3])
dnl AC_SUBST([CURYSQ_LMAX])

dnl AC_ARG_WITH(integer8, AS_HELP_STRING([--with-integer8],[With Fortran integer8.]), [
dnl  if test x$withval = xyes; then
dnl     libcchem_with_integer8=1
dnl  else
dnl 	libcchem_with_integer8=0
dnl  fi], [libcchem_with_integer8=0])
dnl AM_CONDITIONAL(RYSQ_WITH_GAMESS, test x$rysq_with_gamess = x1)

#AM_CPPFLAGS="$AM_CPPFLAGS -I $top_srcdir"

     dnl AC_INIT([arm], [1.0])
     dnl AC_CONFIG_AUX_DIR([.])
     dnl AM_INIT_AUTOMAKE
     dnl AC_PROG_CC
     dnl AC_CONFIG_FILES([Makefile])
     # Call hand's ./configure script recursively.


AC_CONFIG_SUBDIRS([rysq])

AM_CPPFLAGS="$AM_CPPFLAGS -I $PWD/src -I $PWD/rysq/src -I $PWD/src/externals"

AC_SUBST(BASIS_SHELL_MAX, [3])
AC_SUBST(BASIS_SHELL_SP, [1])


AC_SUBST([LIBS])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LDFLAGS])
#AC_SUBST([CUDAFLAGS])

AC_OUTPUT

# Report on what we found.
AC_MSG_NOTICE([])
AC_MSG_NOTICE([**************************************************************])
AC_MSG_NOTICE([ $PACKAGE_NAME configuration:])
AC_MSG_NOTICE([**************************************************************])
AC_MSG_NOTICE([])
dnl AC_MSG_NOTICE([              librysq = $with_rysq])
AC_MSG_NOTICE([                Boost = $with_boost])
AC_MSG_NOTICE([                 HDF5 = $with_hdf5])
AC_MSG_NOTICE([                 BLAS = $with_blas])
AC_MSG_NOTICE([        Global Arrays = $with_ga])
AC_MSG_NOTICE([                  MPI = $with_mpi])
AC_MSG_NOTICE([                 CUDA = $with_cuda])
AC_MSG_NOTICE([               CUBLAS = $with_cublas])
AC_MSG_NOTICE([                  GPU = $with_gpu])
AC_MSG_NOTICE([               GAMESS = $with_gamess])
AC_MSG_NOTICE([            INTEGER*8 = $integer8])
AS_IF([test "x$enable_f77"  =  xyes], [
AC_MSG_NOTICE([                  F77 = $F77])
AC_MSG_NOTICE([               FFLAGS = $FFLAGS])
])
AC_MSG_NOTICE([                  CXX = $CXX])
AC_MSG_NOTICE([             CXXFLAGS = $CXXFLAGS])
AC_MSG_NOTICE([                  CPP = $CPP])
AC_MSG_NOTICE([             CPPFLAGS = $CPPFLAGS])
AC_MSG_NOTICE([              LDFLAGS = $LDFLAGS])
AC_MSG_NOTICE([                    AR = $AR])
AC_MSG_NOTICE([             AR_FLAGS = $AR_FLAGS])
AC_MSG_NOTICE([])

