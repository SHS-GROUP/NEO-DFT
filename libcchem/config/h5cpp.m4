dnl ######################################################################
dnl
dnl Purpose: Determine the locations of C++ Interface for 
dnl          the hdf5 includes and libraries
dnl
dnl Version: $Id: hdf5_cpp.m4,v 1.1 2002/08/13 14:01:30 cary Exp $
dnl
dnl NOTE: This is a modified version of the hdf5.m4 file.
dnl
dnl Tech-X configure system
dnl
dnl Copyright Tech-X Corporation
dnl
dnl ######################################################################

dnl ######################################################################
dnl
dnl Allow the user to specify an overall hdf5 directory.  If specified,
dnl we look for include and lib under this.
dnl
dnl ######################################################################

AC_ARG_WITH(hdf5,[  --with-hdf5=<location of hdf5 installation> ],
            HDF5="$withval", HDF5="")

dnl ######################################################################
dnl
dnl Find hdf5 includes - looking in include location if present,
dnl otherwise in dir/include if present, otherwise in default locations,
dnl first parallel, then serial.
dnl
dnl ######################################################################


if test -n "$HDF5"; then
   CPPFLAGS="$CPPFLAGS -I $HDF5/include $AM_CPPFLAGS"
   export CPPFLAGS
fi
AC_LANG_PUSH(C) 
AC_CHECK_HEADERS([hdf5.h], [HAVE_HDF5_H="yes"])
AC_LANG_POP(C)

if test "x$HAVE_HDF5_H" != "xyes"; then
   AC_MSG_ERROR(*** hdf5.h not found)
else
   AC_DEFINE(HAVE_HDF5, [], [HDF5])
fi

dnl ######################################################################
dnl
dnl See if built parallel
dnl
dnl ######################################################################

# if test $ac_cv_have_hdf5 = yes; then
#   if test -f $HDF5_INCDIR/H5config.h; then
#     hdf5par=`grep "HAVE_PARALLEL 1" $HDF5_INCDIR/H5config.h`
#   elif test -f $HDF5_INCDIR/H5pubconf.h; then
#     hdf5par=`grep "HAVE_PARALLEL 1" $HDF5_INCDIR/H5pubconf.h`
#   fi
# fi
