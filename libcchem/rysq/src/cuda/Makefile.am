include $(top_srcdir)/config/flymake.am

SUBDIRS = kernel

AM_CPPFLAGS =  -I$(top_srcdir)/src -I$(CUDA_ROOT)/include @AM_CPPFLAGS@

# CLEANFILES = $(BUILT_SOURCES)

noinst_LTLIBRARIES = libcuda.la

libcuda_la_SOURCES = runtime.cpp
libcuda_la_SOURCES += fock.cpp
libcuda_la_SOURCES += eri.cpp
libcuda_la_SOURCES += shell.cpp
#libcuda_la_SOURCES += thread.cpp

libcuda_la_LIBADD = kernel/libkernel.la
libcuda_la_LIBADD += -lcudart

