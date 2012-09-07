AM_CPPFLAGS = -I$(top_srcdir)/src @AM_CPPFLAGS@

include $(top_srcdir)/config/flymake.mk

noinst_LTLIBRARIES = libcc.la
libcc_la_SOURCES = cc.cpp
libcc_la_SOURCES += triples.cpp
libcc_la_SOURCES += direct.cpp
libcc_la_SOURCES += integrals.cpp
libcc_la_SOURCES += mp2.cpp
libcc_la_SOURCES += sd.cpp
libcc_la_SOURCES += doubles.cpp

#libcc_la_SOURCES += triples/host.cpp
# if HAVE_CUDA
#    libcc_la_SOURCES += triples/device.cu
# endif
# 

dgemm: dgemm.cpp
	$(CXX) -I../../config $(DEFS) $(AM_CPPFLAGS) $(CPPFLAGS) $(AM_CXXFLAGS) $(CXXFLAGS) $< -o $@
