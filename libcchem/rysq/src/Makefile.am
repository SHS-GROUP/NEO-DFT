include $(top_srcdir)/config/flymake.am
include $(top_srcdir)/config/cheetah.mk

BUILT_SOURCES = normalize.h constants.h

lib_LTLIBRARIES = librysq.la $(librysq_gamess_la)
librysq_la_SOURCES = rysq.cpp
librysq_la_SOURCES += instantiate.cpp
librysq_la_SOURCES += shell.cpp quartet.cpp
librysq_la_SOURCES += eri.cpp fock.cpp

if HAVE_CUDA
SUBDIRS_CUDA = cuda
libcuda_la = cuda/libcuda.la
rysq_cuda_hpp = rysq/cuda.hpp rysq/cuda-memory.hpp rysq/cuda-matrix.hpp
endif


SUBDIRS = rysq externals roots kernel $(SUBDIRS_CUDA)
librysq_la_LIBADD = roots/libroots.la
librysq_la_LIBADD += $(libcuda_la)

nobase_include_HEADERS = rysq.hpp rysq/config.hpp rysq/core.hpp
nobase_include_HEADERS += rysq/eri.hpp rysq/fock.hpp
nobase_include_HEADERS += $(rysq_cuda_hpp)

rysq/config.hpp: rysq/config.hpp.in


# BUILT_SOURCES += kernel/quadrature1-impl.hpp kernel/quadrature2-impl.hpp

# kernel/quadrature2.py: $(top_srcdir)/lib/python/rysq-generate.py
# 	$(top_srcdir)/lib/python/rysq-generate.py > $@

# kernel/quadrature2-impl.hpp: kernel/quadrature2.py
