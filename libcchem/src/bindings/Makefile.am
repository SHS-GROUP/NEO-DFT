AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@
include $(top_srcdir)/config/flymake.mk

noinst_LTLIBRARIES = libbindings.la
libbindings_la_SOURCES = bindings.cpp fortran.cpp
libbindings_la_SOURCES += hf.cpp
libbindings_la_SOURCES += dft.cpp
libbindings_la_SOURCES += cc.cpp
# libbindings_la_SOURCES += bindings.cpp 

if LIBCCHEM_WITH_GAMESS
lib_LTLIBRARIES = libcchem_gamess.la
libcchem_gamess_la_SOURCES = gamess.cpp
endif
