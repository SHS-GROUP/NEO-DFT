#ifndef BOOST_NUMERIC_BINDINGS_CUBLAS_EXCEPTION_HPP
#define BOOST_NUMERIC_BINDINGS_CUBLAS_EXCEPTION_HPP

#include <cuda_runtime.h>
#include <cublas.h>
#include <stdexcept>

namespace boost {
namespace numeric {
namespace bindings {
namespace cublas {

    struct cublas_runtime_error : std::runtime_error {
	cublas_runtime_error(const char *what)
	    : std::runtime_error(what) {}
    };

    struct cublas_status_init_failed : cublas_runtime_error {
	cublas_status_init_failed(const char *what)
	    : cublas_runtime_error(what) {}
    };

    struct cublas_status_not_initialized : cublas_runtime_error {
	cublas_status_not_initialized()
	    : cublas_runtime_error("cublas_status_not_initialized") {}
    };
    struct cublas_status_alloc_failed : cublas_runtime_error {
	cublas_status_alloc_failed()
	    : cublas_runtime_error("cublas_status_alloc_failed") {}
    };
    struct cublas_status_invalid_value : cublas_runtime_error {
	cublas_status_invalid_value()
	    : cublas_runtime_error("cublas_status_invalid_value") {}
    };
    struct cublas_status_arch_mismatch : cublas_runtime_error {
	cublas_status_arch_mismatch()
	    : cublas_runtime_error("cublas_status_arch_mismatch") {}
    };
    struct cublas_status_mapping_error : cublas_runtime_error {
	cublas_status_mapping_error()
	    : cublas_runtime_error("cublas_status_mapping_error") {}
    };
    struct cublas_status_execution_failed : cublas_runtime_error {
	cublas_status_execution_failed()
	    : cublas_runtime_error("cublas_status_execution_failed") {}
    };
    struct cublas_status_internal_error : cublas_runtime_error {
	cublas_status_internal_error()
	    : cublas_runtime_error("cublas_status_execution_failed") {}
    };

    cublasStatus status(){ return cublasGetError(); }

    void check_status() {
	cublasStatus status = cublas::status();
	if (status == CUBLAS_STATUS_SUCCESS) return;
	if (status == CUBLAS_STATUS_NOT_INITIALIZED)
	    throw cublas_status_not_initialized();
	if (status == CUBLAS_STATUS_ALLOC_FAILED)
	    throw cublas_status_alloc_failed();
	if (status == CUBLAS_STATUS_INVALID_VALUE)
	    throw cublas_status_invalid_value();
	if (status == CUBLAS_STATUS_ARCH_MISMATCH)
	    throw cublas_status_arch_mismatch();
	if (status == CUBLAS_STATUS_MAPPING_ERROR)
	    throw cublas_status_mapping_error();
	if (status == CUBLAS_STATUS_EXECUTION_FAILED)
	    throw cublas_status_execution_failed();
	if (status == CUBLAS_STATUS_INTERNAL_ERROR)
	    throw cublas_status_internal_error();
    }

}
}
}
}

#endif // BOOST_NUMERIC_BINDINGS_CUBLAS_EXCEPTION_HPP
