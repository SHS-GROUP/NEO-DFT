#ifndef MULTI_ARRAY_DETAIL_HPP
#define MULTI_ARRAY_DETAIL_HPP

#include <cstdlib>
#include <cassert>
#include <algorithm>

#include <boost/mpl/if.hpp>
#include <boost/utility/enable_if.hpp>

namespace detail {

    template<typename T, size_t N>
    BOOST_GPU_ENABLED
    T multiply(const T (&array)[N]) {
	T value = 1;
	for (int i = 0; i < int(N); ++i) {
	    value *= array[i];
	}
	return value;
    }

    template<size_t N, typename T, typename U>
    bool check(const T *dims, const U *strides) {
	U stride = strides[N-1];
	bool check = (stride == 1 && dims[N-1] > 0);
	for (int i = int(N)-2; i >= 0;  -- i) {
	    check = check && (dims[i] > 0); 
	    stride *= dims[i+1];
	    check = check && (stride <= strides[i]);
	}
	return check;
    }

    template<typename T, size_t N>
    struct multi_array_base {
	typedef size_t size_type;
	
	struct {
	    BOOST_GPU_ENABLED
	    operator const size_type*() const { return data_; }
	    // const size_type& operator[](int i) { return data_[i]; }
	private:
	    friend class multi_array_base;
	    size_type data_[N];
	} size;

	BOOST_GPU_ENABLED
	multi_array_base(const size_type *dims, T *data) {
	    size_type strides[N];
	    strides[N-1] = 1;
	    for (int i = N-2; i >= 0; --i) {
		strides[i] = strides[i+1]*dims[i+1];
	    }
	    initialize(dims, data, strides);
	}
	template<typename U>
	BOOST_GPU_ENABLED
	multi_array_base(const size_type *dims, T *data,
			 const U *strides = NULL) {
	    if (strides) initialize(dims, data, strides);
	    else initialize(dims, data);
	}
    protected:
	T *data_;
	size_type strides_[N];
    private:
	template<typename U>
	BOOST_GPU_ENABLED
	void initialize(const size_type *dims, T *data,
			const U *strides) {
	    //assert(check<N>(dims, strides));
	    for (int i = 0; i < int(N); ++i) {
		size.data_[i] = dims[i];
		strides_[i] = strides[i];
	    }
	    data_ = data;
	}
	BOOST_GPU_ENABLED
	void initialize(const size_type *dims, T *data)	{
	    size_type strides[N];
	    strides[N-1] = 1;
	    for (int i = int(N)-2; i >= 0; --i) {
		strides[i] = strides[i+1]*dims[i+1]; 
	    }
	    initialize(dims, data, strides);
	}
    };
	
    template<typename T, size_t N>
    struct multi_array_ref : multi_array_base<T,N> {
	typedef size_t size_type;
	typedef multi_array_ref<T,N-1> reference;
	template<typename U>
	BOOST_GPU_ENABLED
	multi_array_ref(const size_type *dims, T *data,
			const U *strides = NULL)
	    : multi_array_base<T,N>(dims, data, strides) {}
	BOOST_GPU_ENABLED
	reference operator[](int i ) {
	    return reference(this->size+1,
			     this->data_ + i*stride(),
			     this->strides_+1);
	}
    private:
	BOOST_GPU_ENABLED
	size_type stride() const {
	    return this->strides_[0];
	}
    };

    template<typename T>
    struct multi_array_ref<T,1> : multi_array_base<T,1> {
	typedef size_t size_type;
	typedef T* iterator;
	typedef const T* const_iterator;
	template<typename U>
	BOOST_GPU_ENABLED
	multi_array_ref(const size_type *dims, T *data,
			const U *strides = NULL)
	    : multi_array_base<T,1>(dims, data, strides) {}
	BOOST_GPU_ENABLED
	iterator begin() { return this->data_; }
	BOOST_GPU_ENABLED
	iterator end() { return begin() + this->size[0]; }
	BOOST_GPU_ENABLED
	const_iterator begin() const { return this->data_; }
	BOOST_GPU_ENABLED
	const_iterator end() const { return begin() + this->size[0]; }
	BOOST_GPU_ENABLED
	T& operator[](int i) { return this->data_[i]; }
    };

}


#endif // MULTI_ARRAY_DETAIL_HPP
