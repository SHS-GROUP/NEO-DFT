include $(top_srcdir)/config/flymake.mk

lib_LTLIBRARIES = libcchem.la

libcchem_la_SOURCES = runtime.cpp

libcchem_la_LIBADD =
libcchem_la_LIBADD += basis/libbasis.la
libcchem_la_LIBADD += utility/libutility.la
libcchem_la_LIBADD += bindings/libbindings.la

libcchem_la_LIBADD += integrals/libintegrals.la
libcchem_la_LIBADD += hf/libhf.la
libcchem_la_LIBADD += mp2/libmp2.la
libcchem_la_LIBADD += cc/libcc.la
#libcchem_la_LIBADD += ci/libci.la
#libcchem_la_LIBADD += dft/libdft.la
libcchem_la_LIBADD += matrix/libmatrix.la
libcchem_la_LIBADD += array/libarray.la

SUBDIRS = core basis
SUBDIRS += integrals
SUBDIRS += hf
SUBDIRS += mp2
SUBDIRS += cc
#SUBDIRS += dft
#SUBDIRS += ci
SUBDIRS += utility bindings
SUBDIRS += matrix array adapter

# libcchem_la_LIBADD += contrib/libcontrib.la
# SUBDIRS += contrib
