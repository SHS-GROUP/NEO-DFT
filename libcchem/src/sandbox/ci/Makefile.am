AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@

include $(top_srcdir)/config/flymake.mk

noinst_LTLIBRARIES = libci.la
libci_la_SOURCES = fci.cpp

LDADD = $(top_srcdir)/src/libcchem.la

EXTRA_PROGRAMS = test
test_SOURCES = test.cpp

