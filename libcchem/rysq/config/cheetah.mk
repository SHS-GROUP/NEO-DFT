#PYTHONPATH := $(top_srcdir)/lib/python:$(PYTHONPATH)
export PYTHONPATH:=$(top_srcdir)/lib/python:$(PYTHONPATH)

export RYSQ_LMAX=@RYSQ_LMAX@
export RYSQ_DISABLE_SP=@RYSQ_DISABLE_SP@
export RYSQ_WITH_KERNEL_BLOCK=@RYSQ_WITH_KERNEL_BLOCK@
export RYSQ_WITH_KERNEL_UNROLL=@RYSQ_WITH_KERNEL_UNROLL@
export RYSQ_ENABLE_CUDA=@RYSQ_ENABLE_CUDA@
export RYSQ_ENABLE_CUDA_LMAX=@RYSQ_ENABLE_CUDA_LMAX@


%:: %.tmpl
	cheetah fill --env --oext fill $<
	@mv $@.fill $@
