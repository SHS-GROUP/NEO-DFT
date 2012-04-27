AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@

include $(top_srcdir)/config/flymake.mk

noinst_LTLIBRARIES = libcc.la
libcc_la_SOURCES = triples.cpp

#libcc_la_SOURCES += triples/host.cpp
if HAVE_CUDA
    libcc_la_SOURCES += triples/device.cu
endif
