#ifndef GPU_ALLOCATOR_HPP
#define GPU_ALLOCATOR_HPP

#include "gpu/copy.hpp"
#include "gpu/runtime.hpp"
#include <cuda.h>
#include <cuda_runtime.h>

namespace gpu {

    template<typename T>
    struct allocator {

	typedef T* pointer;
	typedef const T* const_pointer;

	static pointer malloc(size_t size) {
	    void *data = pointer();
	    if (size) {
		// std::cout << "malloc" << std::endl;
		throw_(cudaMalloc(&data, sizeof(T)*size));
		// std::cout << "malloc" << data << std::endl;
	    }
	    return static_cast<pointer>(data);
	}

	static void free(void *data) {
	    if (data) {
		// std::cout << "free " << data << std::endl;
		throw_(cudaFree(data));
	    }
	}

	allocator(size_t size = 0) : capacity_(), data_() {
	    resize(size);
	}
	allocator(const allocator &a) : capacity_(), data_() {
	    operator=(a);
	}
	allocator& operator=(const allocator &a) {
	    resize(a.size());
	    copy<device_to_device>(a.begin(), a.end(), data_);
	    return *this;
	}
	~allocator() {
	    free(data_);
	}

	size_t size() const { return size_; }
	void reserve(size_t capacity) {
	    if (capacity_ < capacity) {
		//std::cout << "reserve " << capacity_ << " " << capacity << std::endl;
		if (data_) free(data_);
		data_ = malloc(capacity);
		capacity_ = capacity;
	    }
	}
	void resize(size_t size) {
	    reserve(size);
	    size_ = size;
	}
	void resize(size_t size,  char value) { 
	    resize(size);
	    if (size) throw_(cudaMemset(data_, value, size*sizeof(T)));
	}	    
	void clear() { size_ = 0; }
	pointer begin() { return data_; }
	pointer end() { return data_ + size_; }
	const_pointer begin() const { return data_; }
	const_pointer end() const { return data_ + size_; }

    private:
	size_t size_;
	size_t capacity_;
	pointer data_;
    };
    
}

#endif // GPU_ALLOCATOR_HPP
