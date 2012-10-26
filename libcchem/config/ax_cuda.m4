
AC_DEFUN([AX_CUDA_ROOT],
[
AC_ARG_WITH([cuda],
            [AS_HELP_STRING([--with-cuda@<:@=DIR@:>@],
                            [use cuda (default is yes) - it is possible to specify
 			     the root directory for cuda (optional)])],
            [if test "$withval" = "no"; then
                 CUDA_ROOT="/bin/false"
             elif test "$withval" != "yes"; then
                 CUDA_ROOT="$withval"
             fi])
])


AC_DEFUN([AX_CUDA],
[

AM_PROG_CUDA($1)

# AC_ARG_WITH(cuda-emulation,
#             AS_HELP_STRING([--with-cuda-emulation],[CUDA device emulation]),
#             [],[withval="no"])
# if test "$withval" = "yes"; then
#    CUDAFLAGS="$CUDAFLAGS --device-emulation --device-compilation='C++'"
# fi

# ## Cuda architecture
# AC_ARG_WITH(cuda-arch,
#             AS_HELP_STRING([--with-cuda-arch[=arch]],[CUDA device arch]),
#             [with_cuda_arch=$withval],[with_cuda_arch=""])
# # Set compiler flags
# if test "x$with_cuda_arch" != "x"; then
#     CUDAFLAGS="$CUDAFLAGS -arch=$with_cuda_arch"
# fi

## verbos if debugging enable
if test "$enable_debug" = "yes" ; then
    CUDAFLAGS="$CUDAFLAGS -g -O0 --ptxas-options='-v'"
fi

CUDA_CPPFLAGS=""
if test -d "$CUDA_ROOT"; then
    CUDA_CPPFLAGS="-I$CUDA_ROOT/include"
fi

if test -n "$CUDA"; then
    if test -z "$with_cuda"; then
        with_cuda="yes";
    fi
    AC_DEFINE([HAVE_CUDA], [1], [have CUDA])
    AC_DEFINE([HAVE_GPU], [1], [with GPU])
else
    with_cuda="no"
fi
AM_CONDITIONAL(HAVE_CUDA, test "$with_cuda" != "no")
AC_SUBST([HAVE_CUDA], `expr "$with_cuda" != "no"`)


])
