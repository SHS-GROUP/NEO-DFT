#define RYSQ_CUDA_KERNEL_SHELL_A (rysq::P)
#define RYSQ_CUDA_KERNEL_SHELL_B (rysq::G)
#include "cuda/kernel/kernel-xx.hpp"
#undef RYSQ_CUDA_KERNEL_SHELL_A // (rysq::P)
#undef RYSQ_CUDA_KERNEL_SHELL_B // (rysq::G)
