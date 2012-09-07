#define RYSQ_CUDA_KERNEL_SHELL_A (rysq::F)
#define RYSQ_CUDA_KERNEL_SHELL_B (rysq::S)
#include "cuda/kernel/kernel-xx.hpp"
#undef RYSQ_CUDA_KERNEL_SHELL_A // (rysq::F)
#undef RYSQ_CUDA_KERNEL_SHELL_B // (rysq::S)
