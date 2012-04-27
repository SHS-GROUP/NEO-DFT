#ifndef GPU_MULTI_ARRAY_HPP
#define GPU_MULTI_ARRAY_HPP

#include <boost/multi_array.hpp>
#include "gpu/allocator.hpp"

namespace gpu {

    template<typename T, size_t N>
    struct multi_array : boost::multi_array_ref<T,N> {
	struct gpu_tag;
	typedef boost::multi_array_ref<T,N> base_type;

	typedef T value_type;
	typedef T& reference;
	typedef const T& const_reference;

	multi_array() : base_type(NULL, extents())
	{
	    // std::cout << "create" << std::endl;
	}
	template<class A>
	multi_array(const A &dims)
	    : base_type(NULL, extents())
	{
	    //std::cout << "create" << std::endl;
	    resize(dims);
	}

	template<typename U>
	void resize(const U (&dims)[N]) {
	    boost::array<U,N> dims_;
	    std::copy(dims, dims + N, dims_.begin());
	    resize(dims_);
	}

	template<typename U>
	void resize(const boost::array<U,N> &dims) {
	    size_t size = 1;
	    boost::array<size_t,N> shape;
	    for (size_t i = 0; i < N; ++i)  {
		size *= dims[i];
		shape[N-(i+1)] = dims[i];
	    }
	    data_.clear();
	    data_.resize(size, 0);
	    // update base_type parent
	    set_base_ptr(&data_.begin()[0]);
	    this->num_elements_ = size;
	    reshape(shape);
	}
	size_t size() const { return data_.size(); }
	size_t size(size_t i) const {
	    if (!(i < N)) throw std::range_error("invalid tensor dimension");
	    return this->shape()[N-(i+1)];
	}
	template<size_t I> size_t size() const { return size(I); }

	// multi_array& fill(const T &value) {
	//     std::fill(data_.begin(), data_.end(), value);
	//     return *this;
	// }

	// typedef boost::multi_array_types::index_range index_range;
	// boost::detail::multi_array::index_genarray<index_range,N> indices() const {
	// 	boost::array<index_range,N> indices;
	// 	for (int i = 0; i < N; ++i) {
	// 	    indices[i] = 
	// 	}

    private:
	typedef boost::detail::multi_array::extent_gen<N> extents;
	gpu::allocator<T> data_;
    };


}


#endif // GPU_MULTI_ARRAY_HPP
