#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ(2.61)
AC_INIT(librysq, 0.1, BUG-REPORT-ADDRESS)
AC_CONFIG_AUX_DIR([config])
AM_INIT_AUTOMAKE([foreign -Wall])
AC_CONFIG_HEADER([config/config.h])

LT_INIT([disable-shared])
AC_CONFIG_MACRO_DIR([config/m4])

# store current user given compiler flags to avoid default setup via AC_PROG_CXX
OLD_CXXFLAGS=$CXXFLAGS
OLD_CFLAGS=$CFLAGS

# check for C++ preprocessor and compiler and the library compiler
# (might change the compiler flags)
AM_DISABLE_SHARED
AM_PROG_LIBTOOL

AC_PROG_CC
AC_PROG_CPP
AC_PROG_CXX
AC_PROG_CXXCPP

# other programs
AC_CHECK_PROG(CHEETAH, [cheetah], [cheetah], [])

# reset compiler flags to initial flags
CXXFLAGS=$OLD_CXXFLAGS
CFLAGS=$OLD_CFLAGS


# Checks for libraries.
# FIXME: Replace `main' with a function in `-lcudart':
#AC_CHECK_LIB([cudart], [main])
# FIXME: Replace `main' with a function in `-lrysq':
#AC_CHECK_LIB([rysq], [main])

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
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_CHECK_FUNCS([gettimeofday sqrt])

AC_OPENMP
CFLAGS="$CFLAGS $OPENMP_CFLAGS"
CXXFLAGS="$CXXFLAGS $OPENMP_CFLAGS"

# Cuda
sinclude(config/ax_cuda.m4)
AX_CUDA_ROOT()
AX_CUDA()

#
#        some compiler options affecting math results are
#        -ftz=true -prec-div=false -prec-sqrt=false --fmad=false
#
AC_ARG_WITH(gpu,
	AC_HELP_STRING([--with-gpu], [Select GPU: tesla, fermi, or kepler]),
	[with_gpu=$withval], [with_gpu="tesla"])
if test "x$with_gpu" = "xtesla"; then
    CUDAFLAGS="-arch=sm_13  $CUDAFLAGS"
elif test "x$with_gpu" = "xfermi"; then
    CUDAFLAGS="-arch=sm_20 $CUDAFLAGS"
elif test "x$with_gpu" = "xkepler"; then
    CUDAFLAGS="-arch=sm_35 $CUDAFLAGS"
else
    echo "*** unknown GPU:" $with_gpu 
fi
CUDAFLAGS="--ptxas-options=-v $CUDAFLAGS"


# DOXYGEN SUPPORT
sinclude(config/ax_prog_doxygen.m4)

# doxygen features
DX_HTML_FEATURE(ON)
DX_CHM_FEATURE(OFF)
DX_CHI_FEATURE(OFF)
DX_MAN_FEATURE(OFF)
DX_RTF_FEATURE(OFF)
DX_XML_FEATURE(OFF)
DX_PDF_FEATURE(ON)
DX_PS_FEATURE(OFF)

# initialize doxygen
DX_INIT_DOXYGEN($PACKAGE_NAME, [doc/Doxyfile], [doc])

# boost
sinclude(config/ax_boost_base.m4)
AX_BOOST_BASE([])
AM_CPPFLAGS="$AM_CPPFLAGS $BOOST_CPPFLAGS"

AC_CONFIG_SRCDIR([src/rysq.hpp.in])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/roots/Makefile
		 src/kernel/Makefile
                 src/bindings/Makefile
		 src/externals/Makefile
                 src/cuda/Makefile
		 src/cuda/kernel/Makefile
                 src/rysq-core.hpp
                 src/rysq.hpp])

# debug compilation
AC_ARG_ENABLE(debug,
    AC_HELP_STRING(--enable-debug, [Debug compilation (Default = no)]),
    enable_debug=$enableval, enable_debug=no)
# set compiler flags
if test "$enable_debug" = "yes" ; then
    CFLAGS="$CFLAGS  -g -O0 -Wall -Wno-uninitialized"
    CXXFLAGS="$CXXFLAGS -g -O0 -Wall -Wno-uninitialized"
    # CUDAFLAGS="-g -G0 $CUDAFLAGS"
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

build_sources=0

## maximum angular momentum
AC_ARG_WITH(lmax, AS_HELP_STRING([--with-lmax],[Max angular momentum]),
              [with_lmax=$withval], [with_lmax=yes])
# substitute in make
if test "$with_lmax" != "yes"; then
    build_sources=1
    echo "$with_lmax" > config/lmax
    echo "*** Lmax = $with_lmax ***" 
else
    with_lmax=
fi
#AC_SUBST(RYSQ_LMAX, ["$with_lmax"])
AC_SUBST_FILE([RYSQ_LMAX])
RYSQ_LMAX=$srcdir/config/lmax


dnl ## disable SP functions
dnl AC_ARG_ENABLE(sp, AS_HELP_STRING([--disable-sp],[disables SP orbitals]),
dnl               [disable_sp=$enableval], [disable_sp=])
dnl # substitute in make
dnl if test "$disable_sp" = "yes"; then
dnl    build_sources=1
dnl    echo "*** No SP functions ***"
dnl fi
dnl AC_SUBST(RYSQ_DISABLE_SP, ["$disable_sp"])

sources="src/rysq-types.hpp src/normalize.h src/constants.h"
sources="$sources src/kernel/quadrature1-impl.hpp src/kernel/quadrature2-impl.hpp"
touch $sources # prevent unnecessary rebuilds
if test "$build_sources" = "1"; then
    echo "*** Templates will be rebuilt: $sources"
    rm -f $sources
    if test -z "$CHEETAH"; then
        AC_MSG_ERROR([cheetah not found])
    fi
fi
AM_CONDITIONAL(RYSQ_BUILD_SOURCES, test "$rysq_build_sources" = "1")


dnl with_sse=no
with_kernel_block=10
with_int2d_block=1024

## with SSE instructions
AC_ARG_WITH(sse, AS_HELP_STRING([--with-sse],[with SSE]),
              [with_sse=$withval], [with_sse=no])

dnl ## kernel blocking
dnl AC_ARG_WITH(kernel-block, AS_HELP_STRING([--with-kernel-block], [kernel block]),
dnl               [with_kernel_block=$withval], [with_kernel_block=10])

dnl ## int2d-block ing
dnl AC_ARG_WITH(int2d-block, AS_HELP_STRING([--with-int2d-block], [2d ints blocking]),
dnl               [with_int2d_block=$withval], [with_int2d_block=1024])

if test "$with_sse" = "yes" ; then
    AC_DEFINE(RYSQ_WITH_SSE, 1, [with SSE])
fi
AC_SUBST(RYSQ_WITH_KERNEL_BLOCK, ["$with_kernel_block"])
AC_DEFINE_UNQUOTED(RYSQ_WITH_INT2D_BLOCK, $with_int2d_block, [int2d blocking])


AC_SUBST_FILE([RYSQ_CORE_HPP])
RYSQ_CORE_HPP=$srcdir/src/rysq-core.hpp

AC_SUBST_FILE([RYSQ_ERI_HPP])
RYSQ_ERI_HPP=$srcdir/src/rysq-eri.hpp

AC_SUBST_FILE([RYSQ_FOCK_HPP])
RYSQ_FOCK_HPP=$srcdir/src/rysq-fock.hpp


AM_CPPFLAGS="$AM_CPPFLAGS -I\$(top_srcdir)/src/externals"

FLYMAKE_CPPFLAGS="\$(DEFAULT_INCLUDES) \$(INCLUDES) \$(AM_CPPFLAGS) \$(CPPFLAGS) \$(AM_CXXFLAGS) \$(CXXFLAGS)"

if test -n "$CUDA_ROOT"; then
    FLYMAKE_CPPFLAGS="$FLYMAKE_CPPFLAGS -I $CUDA_ROOT/include"
fi
AC_SUBST([FLYMAKE_CPPFLAGS])

AC_SUBST([LIBS])
AC_SUBST([AM_CPPFLAGS])
AC_SUBST([AM_LDFLAGS])

## output
AC_OUTPUT

