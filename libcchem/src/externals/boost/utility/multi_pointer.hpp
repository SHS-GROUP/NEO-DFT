#ifndef BOOST_MULTI_POINTER_HPP
#define BOOST_MULTI_POINTER_HPP

#include <cstddef>
#include <boost/mpl/if.hpp>

namespace boost {

    template<typename T, size_t N>
    struct multi_pointer {
	typedef ptrdiff_t difference_type;
	typedef typename boost::mpl::if_c<(N > 1),
	    multi_pointer<T,N-1>, T*>::type reference;
	multi_pointer(T *data, const size_t (&size)[N], 
		      const difference_type (&base)[N]) {
	    data_ = data;
	    for (int i = 0; i < N; ++i) {
	        size.data_[i] = size[i];
	        base.data_[i] = base[i];
	    }
	}
	reference operator[](size_t i) const {
	    return make_reference(data_ + index(i)*size[i], *this);
	}
    private:
	template<typename U>
	struct array {
	    typedef U type[N];
	    operator const type&() const { return data_; }
	private:
	    friend class multi_pointer;
	    type data_;
	};
    public:
	array<size_t> size;
	array<difference_type> base;
    private:
	T *data_;
	difference_type index(size_t i) const {
	    return (difference_type(i) - base[i]);
	}
	template<size_t> struct dim {};
	template<size_t N_>
	static multi_pointer<T,N_-1>
	make_reference(T *data, const multi_pointer<T,N_> &pointer) {
	    return multi_pointer<T,N_-1>(data, pointer.size, pointer.base);
	}
	static T* make_reference(T *data, const multi_pointer<T,1> &pointer) {
	    return data;
	}
    };

}

#endif // BOOST_MULTI_POINTER_HPP
