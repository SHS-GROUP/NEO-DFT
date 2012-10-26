#ifndef CUBLAS_HPP
#define CUBLAS_HPP

#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>

#include "cuda.hpp"
#include <cublas_v2.h>

namespace cublas {

    typedef cudaStream_t stream;
    typedef cublasHandle_t handle_t;

    struct status {
    private:
	cublasStatus value_;
    public:
	status(cublasStatus status = cublasGetError())
	    : value_(status) {
	    //std::cout << status << std::endl;
	}
	operator cublasStatus() const { return value_; }
	std::string str() const {
	    if (value_ == CUBLAS_STATUS_SUCCESS)
		return "CUBLAS_STATUS_SUCCESS";
	    if (value_ == CUBLAS_STATUS_NOT_INITIALIZED)
		return "CUBLAS_STATUS_NOT_INITIALIZED";
	    if (value_ == CUBLAS_STATUS_ALLOC_FAILED)
		return "CUBLAS_STATUS_ALLOC_FAILED";
	    if (value_ == CUBLAS_STATUS_INVALID_VALUE)
		return "CUBLAS_STATUS_INVALID_VALUE";
	    if (value_ == CUBLAS_STATUS_ARCH_MISMATCH)
		return "CUBLAS_STATUS_ARCH_MISMATCH";
	    if (value_ == CUBLAS_STATUS_MAPPING_ERROR)
		return "CUBLAS_STATUS_MAPPING_ERROR";
	    if (value_ == CUBLAS_STATUS_EXECUTION_FAILED)
		return "CUBLAS_STATUS_EXECUTION_FAILED";
	    if (value_ == CUBLAS_STATUS_INTERNAL_ERROR)
		return "CUBLAS_STATUS_INTERNAL_ERROR";
	    return "Unknown CUBLAS error";
	}
    };

    struct error : std::runtime_error {
	static std::string str(int line) {
	    std::stringstream ss;
	    ss << __FILE__ << ":" << line;
	    return ss.str();
	}
	error(const status &s, int line)
	    : std::runtime_error(str(line) + ": " + s.str()) {}
	error(const cudaError_t &s, int line)
	    : std::runtime_error(str(line) + ": " +
				 cudaGetErrorString(s)) {}
    };

#define CUBLAS_CHECK_STATUS(expr) {					\
	status s__ = (expr);						\
	if (s__ != CUBLAS_STATUS_SUCCESS) throw error(s__, __LINE__);	\
    }

    inline void init() {
	CUBLAS_CHECK_STATUS(cublasInit());
    }

    inline void shutdown() {
	CUBLAS_CHECK_STATUS(cublasShutdown());
    }

    template<typename T>
    T* alloc(int n) {
	void *ptr;
	CUBLAS_CHECK_STATUS(cublasAlloc(n, sizeof(T), &ptr));
	return static_cast<T*>(ptr);
    }

    inline void free(void *ptr) {
	CUBLAS_CHECK_STATUS(cublasFree(ptr));
    }

    inline void set_stream(const cuda::stream &s) {
	CUBLAS_CHECK_STATUS(cublasSetKernelStream(s));
    }

    inline void set_stream(cublas::handle_t h, const cuda::stream &s) {
	CUBLAS_CHECK_STATUS(cublasSetStream(h, s));
	    // (cublasSetKernelStream(s));
    }

    template<typename T>
    void set_vector(int n, const T *x, int incx, T *y, int incy,
		    const cuda::stream &s = cuda::stream(NULL)) {
	CUBLAS_CHECK_STATUS
	    ((s != cuda::stream(NULL))
	     ? cublasSetVectorAsync(n, sizeof(T), x, incx, y, incy, s)
	     : cublasSetVector(n, sizeof(T), x, incx, y, incy));
    }

    template<typename T>
    void set_vector(int n, const T *x, T *y,
		    const cuda::stream &s = cuda::stream(NULL)) {
	set_vector(n, x, 1, y, 1, s);
    }

    template<typename T>
    void get_vector(int n, const T *x, int incx, T *y, int incy,
		    const cuda::stream &s = cuda::stream(NULL)) {
	CUBLAS_CHECK_STATUS
	    ((s != cuda::stream(NULL))
	     ? cublasGetVectorAsync(n, sizeof(T), x, incx, y, incy, s)
	     : cublasGetVector(n, sizeof(T), x, incx, y, incy));
    }

    template<typename T>
    void get_vector(int n, const T *x, T *y,
		    const cuda::stream &s = cuda::stream(NULL)) {
	get_vector(n, x, 1, y, 1, s);
    }

    template<typename T>
    void clear(int n, T *x) {
	cuda::memset(x, char(0), n*sizeof(T));
    }

    inline void axpy(cublasHandle_t h, int n, double alpha,
		     const double *x, double *y) {
	cublasDaxpy(h, n, &alpha, x, 1, y, 1);
	CUBLAS_CHECK_STATUS(cublasGetError());
    }

    inline void copy(cublasHandle_t h, int n, const double *x, double *y) {
	cublasDcopy(h, n, x, 1, y, 1);
	CUBLAS_CHECK_STATUS(cublasGetError());
    }

#undef CUBLAS_CHECK_STATUS

}

#endif // CUBLAS_HPP
