AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@

include $(top_srcdir)/config/flymake.mk

noinst_LTLIBRARIES = libdft.la
libdft_la_SOURCES = dft.cpp
libdft_la_SOURCES += xc/xc.cpp
libdft_la_SOURCES += grid/evaluate.cpp

