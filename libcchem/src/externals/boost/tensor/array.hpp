#ifndef TENSOR_ARRAY_HPP
#define TENSOR_ARRAY_HPP

#include <boost/multi_array.hpp>

template<size_t N, typename T>
struct tensor_array : boost::multi_array_ref<T,N> {
    typedef boost::multi_array_ref<T,N> base_type;

    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;

    tensor_array() : base_type(NULL, extents())
    {
	// std::cout << "create" << std::endl;
    }
    template<class A>
    tensor_array(const A &dims)
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
	set_base_ptr(&data_[0]);
	this->num_elements_ = size;
	reshape(shape);
    }
    size_t size() const { return data_.size(); }
    size_t size(size_t i) const {
	if (!(i < N)) throw std::range_error("invalid tensor dimension");
	return this->shape()[N-(i+1)];
    }
    template<size_t I> size_t size() const { return size(I); }

    tensor_array& fill(const T &value) {
	std::fill(data_.begin(), data_.end(), value);
	return *this;
    }

    // typedef boost::multi_array_types::index_range index_range;
    // boost::detail::multi_array::index_genarray<index_range,N> indices() const {
    // 	boost::array<index_range,N> indices;
    // 	for (int i = 0; i < N; ++i) {
    // 	    indices[i] = 
    // 	}

    void put(const T *buffer, const size_t *start, const size_t *stop) {
	boost::array<size_t,N> dims;
	static const size_t K = N-1;
	for (size_t i = 0; i < N; ++i) {
	    dims[K-i] = stop[i] - start[i];
	}
	boost::const_multi_array_ref<T,N> B(buffer, dims);
	base_type::operator[](indices(start+K, stop+K, boost::indices)) = B;
    }

    void get(T *buffer, const size_t *start, const size_t *stop) const {
	boost::array<size_t,N> dims;
	static const size_t K = N-1;
	for (size_t i = 0; i < N; ++i) {
	    dims[K-i] = stop[i] - start[i];
	}
	boost::multi_array_ref<T,N> B(buffer, dims);
	B = base_type::operator[](indices(start+K, stop+K, boost::indices));
    }

private:
    typedef boost::detail::multi_array::extent_gen<N> extents;
    std::vector<T> data_;
private:
    template<typename I, class G>
    static boost::detail::multi_array::index_gen<N,N>
    indices(I start, I stop, const G &gen) {
	boost::multi_array_types::index_range range(*start--, *stop--);
	return indices(start, stop, gen[range]);
    }
    template<typename I>
    static boost::detail::multi_array::index_gen<N,N>
    indices(I start, I stop,
	    const boost::detail::multi_array::index_gen<N,N> &indices) {
	return indices;
    }
    
};


#endif // TENSOR_ARRAY_HPP
