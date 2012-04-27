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


sinclude([config/h5cpp.m4])


# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST
AC_C_INLINE
AC_TYPE_INT8_T
AC_TYPE_SIZE_T
AC_HEADER_TIME

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday sqrt])

AC_OPENMP
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"

# boost
sinclude(config/ax_boost_base.m4)
AX_BOOST_BASE([])
AM_CPPFLAGS="$AM_CPPFLAGS $BOOST_CPPFLAGS"

# Cuda
sinclude(config/ax_cuda.m4)
AX_CUDA_ROOT()
AX_CUDA()

AM_CPPFLAGS="$AM_CPPFLAGS $CUDA_CPPFLAGS"
dnl echo "$AM_CPPFLAGS"

AC_ARG_WITH(gpu,
	AC_HELP_STRING([--with-gpu], [Select GPU: tesla or fermi]),
	[with_gpu=$withval], [with_gpu="tesla"])
if test "x$with_gpu" = "xtesla"; then
    CUDAFLAGS="-arch=sm_13 -Xptxas -dlcm=cg $CUDAFLAGS"
elif test "x$with_gpu" = "xfermi"; then
    CUDAFLAGS="-arch=sm_21 -Xptxas -dlcm=cg -maxrregcount=80 $CUDAFLAGS"
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
		src/core/Makefile
		src/hf/Makefile
		src/mp2/Makefile
		src/cc/Makefile
		src/dft/Makefile
		src/utility/Makefile
		src/adapter/Makefile
		src/matrix/Makefile
		src/multi_array/Makefile
		src/bindings/Makefile
		src/phoenix/Makefile
		src/contrib/Makefile])

# debug compilation
AC_ARG_ENABLE(debug,
    AC_HELP_STRING(--enable-debug, [Debug compilation (Default = no)]),
    enable_debug=$enableval, enable_debug=no)

if test "$enable_debug" = "yes" ; then
    CFLAGS="$CFLAGS -O0 -g -Wall -Wno-uninitialized"
    CXXFLAGS="$CXXFLAGS -O0 -g -Wall -Wno-uninitialized"
else
    CPPFLAGS="$CPPFLAGS -O3 -DNDEBUG"
fi


## enable ga
AC_ARG_WITH(ga, AC_HELP_STRING([--with-ga=PATH], [where ga is installed]),
            [with_ga=$withval], [with_ga="no"])
AM_CONDITIONAL(_HAVE_GA, [test "$with_ga" != "no"])
if test "$with_ga" != "no"; then
   if test "$with_ga" != "yes"; then
       AM_CPPFLAGS="$AM_CPPFLAGS -I$with_ga/include"
       #LDFLAGS="$LDFLAGS -L$with_ga/lib"
   fi
   #LIBS="$LIBS -lglobal -lga++"
   with_ga=1
   #AC_DEFINE(_HAVE_GA, 1, [enable ga])
else
   with_ga=0
fi
AC_SUBST(_HAVE_GA, ["$with_ga"])

# MKL
AC_ARG_WITH(mkl, AC_HELP_STRING([--with-mkl], [Use MKL Blas]),
		 [with_mkl=$withval], [with_mkl="no"])
AM_CONDITIONAL(HAVE_MKL, [test "$with_mkl" != "no"])
if test "$with_mkl" != "no"; then
    AC_DEFINE(HAVE_MKL, 1, [with MKL])
fi
    

# MKL
AC_ARG_WITH(cublas,
	AC_HELP_STRING([--with-cublas], [Use Cuda Blas]),
	[with_cublas=$withval], [with_cublas=""])
if test "$HAVE_CUDA" != "0"; then
    echo "$HAVE_CUDA"
    if test "$with_cublas" != "no"; then
        with_cublas="$CUDA_ROOT"
    fi
else
    with_cublas="no"
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
    echo "CUBLAS=$with_cublas"
fi



# Rysq settings
AC_ARG_WITH(rysq, 
            AC_HELP_STRING([--with-rysq=PATH], [where rysq is installed]))

if test "$with_rysq" != "no"; then
   if test "$with_rysq" != "yes"; then
      AM_CPPFLAGS="$AM_CPPFLAGS -I$with_rysq/include"
      LDFLAGS="$LDFLAGS -L$with_rysq/lib"
   fi
#   LIBS="$LIBS -lrysq"
fi  


# gamess bindings
AC_ARG_WITH(gamess,
	AS_HELP_STRING([--with-gamess],[Adds GAMESS support to libcchem]),
	[if test "$withval" = "yes"; then
	    AC_DEFINE([LIBCCHEM_WITH_GAMESS], [1], [Compile with GAMESS interface])
	 fi])
AM_CONDITIONAL(LIBCCHEM_WITH_GAMESS, test "$with_gamess" = "yes")

# gamess limits
AC_ARG_WITH(gamess-mxao,
	AS_HELP_STRING([--with-gamess-mxao],[GAMESS MXAO parameter]),
        AC_DEFINE([LIBCCHEM_WITH_GAMESS_MXAO], [$withval], [GAMESS MXAO parameter]),
        AC_DEFINE([LIBCCHEM_WITH_GAMESS_MXAO], [8192], [GAMESS MXAO parameter]))

# integer8 Fortran bindings
AC_ARG_WITH(integer8,
	AS_HELP_STRING([--with-integer8],[Adds Integer8 Fortran support]),
	[if test x$withval = xyes; then
	    AC_DEFINE([LIBCCHEM_WITH_INTEGER8], [1], [8-byte Fortran integer])
	 fi])

		

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
