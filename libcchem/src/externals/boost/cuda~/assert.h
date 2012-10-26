#ifndef _CUDA_ASSERT_H_
#define _CUDA_ASSERT_H_

#include <cuda.h>
#include <stdio.h>

#define cuda_assert(...) {					\
	cudaError_t status = cudaGetLastError();		\
	if (status != cudaSuccess) {				\
	    fprintf(stderr, "%s:%i: %s\n", __FILE__, __LINE__,	\
		    cudaGetErrorString(status));		\
	    assert(status == cudaSuccess);			\
	}							\
    }

#endif /* _CUDA_ASSERT_H_ */
