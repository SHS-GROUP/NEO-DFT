include $(top_srcdir)/config/flymake.am
include $(top_srcdir)/config/cheetah.mk

BUILT_SOURCES = rysq-types.hpp

lib_LTLIBRARIES = librysq.la
librysq_la_SOURCES = rysq.cpp
librysq_la_SOURCES += instantiate.cpp
librysq_la_SOURCES += shell.cpp quartet.cpp
librysq_la_SOURCES += eri.cpp fock.cpp

nobase_include_HEADERS = rysq.hpp

SUBDIRS = externals roots kernel
librysq_la_LIBADD = roots/libroots.la

if HAVE_CUDA
SUBDIRS += cuda
librysq_la_LIBADD += cuda/libcuda.la
nobase_include_HEADERS += rysq-cuda.hpp
endif

