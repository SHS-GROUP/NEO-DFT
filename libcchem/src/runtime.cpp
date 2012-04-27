#include "runtime.hpp"

#include "config.h"
#if defined(HAVE_CUDA) || defined(HAVE_CUBLAS)
#warning "Compiling BOOST CUDA implementation"
#include "boost/cuda/detail/implementation.hpp"
#endif

std::auto_ptr<runtime> runtime::rt_;
