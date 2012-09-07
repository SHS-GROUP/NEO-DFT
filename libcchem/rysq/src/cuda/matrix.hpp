#ifndef RYSQ_CUDA_MATRIX_HPP
#define RYSQ_CUDA_MATRIX_HPP

#include <stdlib.h>

#if !defined(__host__)  && !defined(__device__)
#define __host__
#define __device__
#endif

namespace rysq {

    struct matrix_layout {
	matrix_layout() {}
	matrix_layout(size_t size1, size_t size2, size_t ld = 0)
	    : size1(size1),size2(size2), ld_((ld) ? ld : size1), size(ld_*size2) {}
	__host__ __device__
	int element_at(int i,int j)  const { return i + j* ld_;}
	size_t size1, size2, ld_, size;
    };

    template<typename T>
    class matrix_data_array {
	matrix_layout layout_;
	T *data_;
	size_t size_;
    public:
	matrix_data_array() {}
	matrix_data_array(const matrix_layout &layout, T *data, size_t size)
	    : layout_(layout), data_(data), size_(size){}
	__host__ __device__
	const matrix_layout& layout() const { return layout_; }
	__host__ __device__
	T* operator[](int i) { return data_ + i*layout_.size; }
	__host__ __device__
	const T* operator[](int i) const { return data_ + i* layout_.size; }
	__host__ __device__
	size_t size()const { return size_; }
    };

    template<typename T>
    class matrix {
	matrix_layout layout_;
	T *data_;
    public:
	matrix() {}
	matrix(size_t size1, size_t size2, T *data,size_t ld = 0)
	    : layout_(size1, size2, ld), data_(data) {}
	size_t size1 () const { return layout_.size1; }
	size_t size2 () const { return layout_.size2; }
	__host__ __device__
	T& operator()(int i,int j) { return data_[layout_.element_at(i,j)]; }
	__host__ __device__
	const T& operator()(int i,int j) const { return data_[layout_.element_at(i,j)]; }
    };



}

#endif /* RYSQ_CUDA_MATRIX_HPP */
