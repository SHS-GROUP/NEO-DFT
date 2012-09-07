#define RYSQ_CUDA_KERNEL_SHELL_A (rysq::SP)
#define RYSQ_CUDA_KERNEL_SHELL_B (rysq::G)
#include "cuda/kernel/kernel-xx.hpp"
#undef RYSQ_CUDA_KERNEL_SHELL_A // (rysq::SP)
#undef RYSQ_CUDA_KERNEL_SHELL_B // (rysq::G)
