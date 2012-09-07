include $(top_srcdir)/config/flymake.mk

EXTRA_PROGRAMS = cublas
cublas_SOURCES = cublas.cpp
cublas_LDADD = -L/usr/local/cuda/lib64 -lcublas -lcudart 

# bench: prog$(EXEEXT)
#         execute prog$(EXEEXT)...
# .PHONY: bench 