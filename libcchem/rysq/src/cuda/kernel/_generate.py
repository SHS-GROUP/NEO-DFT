#!/usr/bin/python

shells = ("s", "p", "sp", "d", "f")#, "g")
files = []

for b in shells:
    for a in shells[1:]:
	f = open("kernel-%s.cu" % (a+b), 'w')
	A = "rysq::%s" % a.upper()
	B = "rysq::%s" % b.upper()
	f.write("#define RYSQ_CUDA_KERNEL_SHELL_A (%s)\n" % A)
	f.write("#define RYSQ_CUDA_KERNEL_SHELL_B (%s)\n" % B)
	f.write('#include "cuda/kernel/kernel-xx.hpp"\n')
	f.write("#undef RYSQ_CUDA_KERNEL_SHELL_A // (%s)\n" % A)
	f.write("#undef RYSQ_CUDA_KERNEL_SHELL_B // (%s)\n" % B)

	files.append(f.name)

print "libkernel_la_SOURCES += \\"
print "  \\\n".join(files)
	
