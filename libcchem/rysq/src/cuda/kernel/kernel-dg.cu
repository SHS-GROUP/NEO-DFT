#define RYSQ_CUDA_KERNEL_SHELL_A (rysq::D)
#define RYSQ_CUDA_KERNEL_SHELL_B (rysq::G)
#include "cuda/kernel/kernel-xx.hpp"
#undef RYSQ_CUDA_KERNEL_SHELL_A // (rysq::D)
#undef RYSQ_CUDA_KERNEL_SHELL_B // (rysq::G)
