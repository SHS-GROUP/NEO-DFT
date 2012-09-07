#PYTHONPATH := $(top_srcdir)/lib/python:$(PYTHONPATH)
export PYTHONPATH:=$(top_srcdir)/lib/python:$(PYTHONPATH)

# @RYSQ_LMAX@ may be AC_SUBST_FILE
export RYSQ_LMAX=\
@RYSQ_LMAX@

export RYSQ_DISABLE_SP=@RYSQ_DISABLE_SP@
export RYSQ_WITH_KERNEL_BLOCK=@RYSQ_WITH_KERNEL_BLOCK@
export RYSQ_WITH_KERNEL_UNROLL=@RYSQ_WITH_KERNEL_UNROLL@
export HAVE_CUDA=@HAVE_CUDA@

%:: %.tmpl
	cheetah fill --env --oext fill $<
	@mv $@.fill $@
