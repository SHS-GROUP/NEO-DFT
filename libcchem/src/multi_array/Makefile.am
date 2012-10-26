AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@

#LIBS += -lboost_thread

include $(top_srcdir)/config/flymake.mk

noinst_LTLIBRARIES = libmulti_array.la

libmulti_array_la_SOURCES = link.cpp

# if HAVE_CUDA
# libmulti_array_la_SOURCES += gpu/algorithm.cu
# endif

# subdirs = gpu

#
#noinst_PROGRAMS = test
#test_SOURCES = test.cpp
