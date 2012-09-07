#ifndef ARRAY_REF_HPP
#define ARRAY_REF_HPP

#include <cstdlib>

template<typename T, size_t N = 1, typename TP = const T*>
struct const_array_ref {

    template<typename U, size_t N_, size_t disable = N_>
    struct traits {
	typedef U* reference;
    };

    template<typename U, size_t N_>
    struct traits<U,N_,1> {};

    template<typename U>
    struct traits<U,1> {
	typedef U& reference;
    };

    typedef size_t size_type;
    typedef const T value_type;
    typedef TP pointer;
    typedef const T* const_pointer;
    typedef typename traits<T,N>::reference reference;
    typedef typename traits<const T,N>::reference const_reference;
    const_array_ref(size_t size, TP data, size_t stride)
	: size_(size), data_(data), stride_(stride) {}
    const_pointer data() const { return data_; }
    size_type size() const { return size_; }
    const_reference operator[](int i) const {
	return access(data_ + i*stride_, traits<const T, N>());
    }
protected:
    template<typename U>
    static typename traits<U,1>::reference
    access(U *data, traits<U,1>) {
	return *data;
    }
    template<typename U, size_t N_>
    static typename traits<U,N_>::reference
    access(U *data, traits<U,N_>) {
	return data;
    }
    size_t size_;
    TP data_;
    size_t stride_;
};

template<typename T, size_t N = 1>
struct array_ref : const_array_ref<T,N,T*> {
    typedef const_array_ref<T,N,T*> base_type;
    typedef typename base_type::pointer pointer;
    typedef typename base_type::reference reference;
    array_ref(T *data, size_t size)
	: base_type(data, size) {}
    pointer data() { return this->data_; }
    reference operator[](int i) {
	typename base_type::template traits<T,N> traits;
	return base_type::access(this->data_ + i*this->stride_, traits);
    }
};

#endif // ARRAY_REF_HPP
