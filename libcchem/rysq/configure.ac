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
AC_CHECK_PROG(CHEETAH,cheetah,cheetah)
dnl AC_CHECK_PROG(NVCC,nvcc,nvcc)

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


# Cuda
sinclude(config/ax_cuda.m4)
AX_CUDA_ROOT()
AX_CUDA()

AC_ARG_WITH(gpu,
	AC_HELP_STRING([--with-gpu], [Select GPU: tesla or fermi]),
	[with_gpu=$withval], [with_gpu="tesla"])
if test "x$with_gpu" = "xtesla"; then
    CUDAFLAGS="-arch=sm_13 -Xptxas -dlcm=cg $CUDAFLAGS"
elif test "x$with_gpu" = "xfermi"; then
    CUDAFLAGS="-arch=sm_21 -Xptxas -dlcm=ca -ftz=true $CUDAFLAGS"
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

AC_CONFIG_SRCDIR([src/rysq.hpp])

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 src/roots/Makefile
		 src/kernel/Makefile
                 src/bindings/Makefile
		 src/externals/Makefile
                 src/cuda/Makefile
		 src/cuda/kernel/Makefile
		 src/rysq/Makefile
		 src/rysq/config.hpp])

# debug compilation
AC_ARG_ENABLE(debug,
    AC_HELP_STRING(--enable-debug, [Debug compilation (Default = no)]),
    enable_debug=$enableval, enable_debug=no)
# set compiler flags
if test "$enable_debug" = "yes" ; then
    CFLAGS="$CFLAGS  -g -O0 -Wall -Wno-uninitialized"
    CXXFLAGS="$CXXFLAGS -g -O0 -Wall -Wno-uninitialized"
else
    CPPFLAGS="$CPPFLAGS -DNDEBUG"
fi


## maximum angular momentum
AC_ARG_WITH(lmax, AS_HELP_STRING([--with-lmax],[Max angular momentum]),
              [with_lmax=$withval],[with_lmax=4])
# substitute in make
echo "*** Lmax = $with_lmax ***" 
AC_SUBST(RYSQ_LMAX,["$with_lmax"])


## disable SP functions
AC_ARG_ENABLE(sp, AS_HELP_STRING([--disable-sp],[disables SP orbitals]),
              [disable_sp=$enableval], [disable_sp=no])
# substitute in make
AC_SUBST(RYSQ_DISABLE_SP, ["$disable_sp"])


## with SSE instructions
AC_ARG_WITH(sse, AS_HELP_STRING([--with-sse],[with SSE]),
              [with_sse=$withval], [with_sse=no])
# define macro
if test "$with_sse" = "yes" ; then
   AC_DEFINE(RYSQ_WITH_SSE, 1, [with SSE])
fi

## kernel blocking
AC_ARG_WITH(kernel-block, AS_HELP_STRING([--with-kernel-block], [kernel block]),
              [with_kernel_block=$withval], [with_kernel_block=10])
# substitute in make
AC_SUBST(RYSQ_WITH_KERNEL_BLOCK, ["$with_kernel_block"])

## int2d-block ing
AC_ARG_WITH(int2d-block, AS_HELP_STRING([--with-int2d-block], [2d ints blocking]),
              [with_int2d_block=$withval], [with_int2d_block=1024])
# substitute in make
AC_DEFINE_UNQUOTED(RYSQ_WITH_INT2D_BLOCK, $with_int2d_block, [int2d blocking])

AC_SUBST_FILE([CUDAXX_HPP])
CUDAXX_HPP=$srcdir/src/externals/cuda/cudaxx.hpp



## External bindings			 

# with gamess
AC_ARG_WITH(gamess, AS_HELP_STRING([--with-gamess],[With GAMESS bindings]),
	    [with_gamess=$withval], [with_gamess=no])
# define make conditional
AM_CONDITIONAL(RYSQ_WITH_GAMESS, test "$with_gamess" = "yes")

dnl with fortran
AC_ARG_WITH(fortran, AS_HELP_STRING([--with-fortran],[With FORTRAN bindings]),
            [with_fortran=$withval], [with_fortran=no])
# define make conditional
AM_CONDITIONAL(RYSQ_WITH_FORTRAN, test "$with_gamess" = "yes" -o "$with_fortran" = "yes")

# with integer8
AC_ARG_WITH(integer8, AS_HELP_STRING([--with-integer8],[Compile to use INTEGER(8)]),
            [with_integer8=$withval], [with_integer8=no])
# define macro
if test "$with_integer8" = "yes" ; then
   AC_DEFINE(RYSQ_WITH_INTEGER8, 1, [with INTEGER8])
fi

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

