#ifndef MULTI_DETAIL_BASE_HPP
#define MULTI_DETAIL_BASE_HPP

#include <boost/config.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/int.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

namespace array {
namespace detail {

    template<class A>
    struct reference {
	typedef typename boost::mpl::if_<
	    boost::is_const<A>,
	    typename A::const_reference,
	    typename A::reference
	    >::type type;
    };

    template<class A>
    BOOST_GPU_ENABLED
    typename reference<A>::type
    access(A &a, int i, boost::mpl::int_<1>) {
	return a.data()[i];
    }

    template<class A, class K>
    BOOST_GPU_ENABLED
    typename reference<A>::type
    access(A &a, int i, K) {
	BOOST_STATIC_ASSERT((K::value == A::dimensionality));
	static const size_t N = A::dimensionality;
	typedef typename reference<A>::type R;
	return R(a.size+1, a.data() + i*a.strides()[N-1]);
    }

    template<typename T, size_t N, typename P = const T*>
    struct const_adapter_base {

	typedef T element;
	static const size_t dimensionality = N;
	typedef size_t size_type;
	
	struct {
	    typedef size_type A[N];
	    BOOST_GPU_ENABLED
	    //operator const size_type*() const { return data_; }
	    operator const A&() const { return data_; }
	private:
	    friend class const_adapter_base;
	    A data_;
	} size;

	typedef typename boost::mpl::if_c<
	    (N == 1), const T&, const_adapter_base<T,N-1>
	    >::type const_reference;
	typedef const_reference reference;

	typedef const T* const_iterator;

	template<typename U>
	BOOST_GPU_ENABLED
	const_adapter_base(const U *dims, P data) {
	    size_type strides[N];
	    strides[N-1] = 1;
	    for (int i = int(N)-2; i >= 0; --i) {
		strides_[i] = strides[i+1]*dims[i+1];
	    }
	    for (int i = 0; i < int(N); ++i) {
		size.data_[i] = dims[i];
	    }
	    data_ = data;
	}

	BOOST_GPU_ENABLED
	P data() { return data_; }
	BOOST_GPU_ENABLED
	const T* data() const { return data_; }

	BOOST_GPU_ENABLED
	const size_type (&strides() const)[N] {
	    return strides_;
	}

	BOOST_GPU_ENABLED
	const_reference operator[](int i) const {
	    return access(*this, i, boost::mpl::int_<N>());
	}

	BOOST_GPU_ENABLED
	const_iterator begin() const {
	    BOOST_STATIC_ASSERT((N == 1));
	    return this->data_;
	}
	BOOST_GPU_ENABLED
	const_iterator end() const {
	    return begin() + this->size[0];
	}

    protected:
	P data_;
	size_type strides_[N];
    };


    template<typename T, size_t N>
    struct adapter_base : const_adapter_base<T,N,T*> {
	typedef const_adapter_base<T,N,T*> base;

	typedef typename boost::mpl::if_c<
	    (N == 1), T&, adapter_base<T,N-1>
	    >::type reference;
	typedef T* iterator;

	template<typename U>
	BOOST_GPU_ENABLED
	adapter_base(const U *dims, T* data)
	    : base(dims, data) {}

	using base::operator[];

	BOOST_GPU_ENABLED
	reference operator[](int i) {
	    return access(*this, i, boost::mpl::int_<N>());
	}

	BOOST_GPU_ENABLED
	iterator begin() {
	    BOOST_STATIC_ASSERT((N == 1));
	    return this->data_;
	}
	BOOST_GPU_ENABLED
	iterator end() {
	    return begin() + this->size[0];
	}
    };


}
}

#endif // MULTI_DETAIL_BASE_HPP
