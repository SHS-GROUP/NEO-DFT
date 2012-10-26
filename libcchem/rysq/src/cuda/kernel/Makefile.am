#include $(top_srcdir)/config/cheetah.mk
include $(top_srcdir)/config/flymake.am

AM_CPPFLAGS =  -I$(top_srcdir)/src -I$(CUDA_ROOT)/include @AM_CPPFLAGS@

# BUILT_SOURCES = specialization.hpp
CLEANFILES = $(BUILT_SOURCES)

noinst_LTLIBRARIES = libkernel.la

libkernel_la_SOURCES = link.cpp
libkernel_la_SOURCES += eri.cu fock.cu

libkernel_la_SOURCES += kernel-sx.cu
libkernel_la_SOURCES += \
kernel-ps.cu  \
kernel-sps.cu  \
kernel-ds.cu  \
kernel-fs.cu  \
kernel-pp.cu  \
kernel-spp.cu  \
kernel-dp.cu  \
kernel-fp.cu  \
kernel-psp.cu  \
kernel-spsp.cu  \
kernel-dsp.cu  \
kernel-fsp.cu  \
kernel-pd.cu  \
kernel-spd.cu  \
kernel-dd.cu  \
kernel-fd.cu  \
kernel-pf.cu  \
kernel-spf.cu  \
kernel-df.cu  \
kernel-ff.cu


#libkernel_la_SOURCES += reduce.cu

headers = device.hpp quartet.hpp

# .cu.o:
# 	$(NVCC) -o $@ -c $< $(NVCCFLAGS) $(AM_CPPFLAGS)
# #	source='$<' object='$@' libtool=no \
# #	DEPDIR=$(DEPDIR) $(CXXDEPMODE) $(depcomp) \
# #	$(CXXCOMPILE) -c -o $@ $<

# .cu.lo: $(headers)
# 	$(top_srcdir)/config/cudalt.py $@ $(NVCC) -c $(NVCCFLAGS) $(AM_CPPFLAGS) $<