#ifndef CUDA_HPP
#define CUDA_HPP

#include <string>
#include <stdexcept>
#include <iostream>
#include <sstream>

#include <cuda_runtime.h>
#include <cuda.h>

namespace cuda {

#define CUDA_VERIFY(expr) ::cuda::assert_(__FILE__, __LINE__, (expr));

    struct error : std::runtime_error {
	static std::string str(const char *file, int line) {
	    std::stringstream ss;
	    ss << file << ":" << line;
	    return ss.str();
	}
	error(const cudaError_t &s, const char *file, int line)
	    : std::runtime_error(str(file, line) + ": " +
				 cudaGetErrorString(s)) {}
    };

    struct error_no_device : error {
	explicit error_no_device(const char *file, int line)
	    : error(cudaErrorNoDevice, file, line) {}
    }; 

    inline void assert_(const char *file, int line, const cudaError_t &s) {
	if (s != cudaSuccess) {
	    if (s == cudaErrorNoDevice) throw error_no_device(file, line);
	    throw error(s, file, line);
	}
    }

    struct stream;
    struct event;

    struct stream {
	stream() {}
	stream(const cudaStream_t &s) : data_(s) {}
	void create() {
	    CUDA_VERIFY( cudaStreamCreate(&data_) );
	}
	void destroy() {
	    CUDA_VERIFY( cudaStreamDestroy(data_) );
	}
	void synchronize() {
	    CUDA_VERIFY( cudaStreamSynchronize(data_) );
	}
	bool query() {
	    cudaError_t err = cudaStreamQuery(data_);
	    if (err == cudaSuccess) return true;
	    if (err == cudaErrorNotReady) return false;
	    CUDA_VERIFY( err );
	    return false;
	}
	void wait(const event &e);
	operator cudaStream_t() const { return data_; }
    private:
	cudaStream_t data_;
    };


    struct event {
	void create() {
	    CUDA_VERIFY( cudaEventCreate(&data_) );
	}
	void destroy() {
	    CUDA_VERIFY( cudaEventDestroy(data_) );
	}
	void synchronize() {
	    CUDA_VERIFY( cudaEventSynchronize(data_) );
	}
	void record(const stream &s = NULL) {
	    CUDA_VERIFY( cudaEventRecord(data_, s) );
	}
	operator cudaEvent_t() const { return data_; }
    private:
	cudaEvent_t data_;
    };

    inline void stream::wait(const event &e) {
	CUDA_VERIFY( cudaStreamWaitEvent(this->data_, e, 0) );
    }


    typedef cudaDeviceProp device_prop;

    inline void synchronize() {
	CUDA_VERIFY(cudaThreadSynchronize());
    }

    template<typename T>
    inline T* malloc(size_t size) {
	T *ptr = 0;
	CUDA_VERIFY(cudaMalloc(&ptr, size*sizeof(T)));
	return ptr;
    }

    inline void free(void *ptr) {
	CUDA_VERIFY(cudaFree(ptr));
    }

    template<typename T>
    inline T* malloc_host(size_t size) {
	T *ptr = 0;
	CUDA_VERIFY(cudaMallocHost(&ptr, size*sizeof(T)));
	return ptr;
    }

    inline void free_host(void *ptr) {
	CUDA_VERIFY(cudaFreeHost(ptr));
    }

#if (CUDA_VERSION >= 4000)
    inline void host_register(const void *ptr, size_t size) {
	//std::cout <<  "host_register " << ptr << std::endl;
	CUDA_VERIFY(cudaHostRegister((void*)ptr, size, 0));
    }
#endif

#if (CUDA_VERSION >= 4000)
    inline void host_unregister(const void *ptr) {
	//std::cout <<  "host_unregister " << ptr << std::endl;
	CUDA_VERIFY(cudaHostUnregister((void*)ptr));
    }
#endif

    inline void memset(void *ptr, char value, size_t count) {
	CUDA_VERIFY(cudaMemset(ptr, value, count));
    }

    inline int get_device_count() {
	int count;
	CUDA_VERIFY(cudaGetDeviceCount(&count));
	return count;
    }

    inline device_prop get_device_properties(int device) {
	cudaDeviceProp prop;
	CUDA_VERIFY(cudaGetDeviceProperties(&prop, device));
	return device_prop(prop);
    }

    inline void set_device(int device) {
	CUDA_VERIFY(cudaSetDevice(device));
    }

    inline void thread_exit() {
	CUDA_VERIFY(cudaThreadExit());
    }

    typedef enum {
	host_to_host = cudaMemcpyHostToHost,
	host_to_device = cudaMemcpyHostToDevice,
	device_to_host = cudaMemcpyDeviceToHost,
	device_to_device = cudaMemcpyDeviceToDevice
    } memcpy_kind;

    template<typename T>
    struct element {
	static const size_t size = sizeof(T);
    };

    template<>
    struct element<void> {
	static const size_t size = 1;
    };

    template<typename T>
    inline void copy(size_t size, const T *input, T *output,
		     memcpy_kind kind) {
	CUDA_VERIFY(cudaMemcpy(output, input, size*element<T>::size,
			       cudaMemcpyKind(kind)));
    }

    template<typename T>
    inline void copy(size_t size, const T *input, T *output,
		     memcpy_kind kind, const stream &s) {
	CUDA_VERIFY(cudaMemcpyAsync(output, input, size*element<T>::size,
				    cudaMemcpyKind(kind), s));
    }

}

#endif // CUDA_HPP
