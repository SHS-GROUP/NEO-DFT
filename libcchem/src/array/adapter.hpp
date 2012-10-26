#ifndef ARRAY_ADAPTER_HPP
#define ARRAY_ADAPTER_HPP

#include <boost/config.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/print.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/assert.hpp>
#include <boost/static_assert.hpp>
#include <boost/type_traits.hpp>

#include "array/forward.hpp"
#include "array/assign.hpp"

namespace array {
namespace detail {

    template<class A>
    struct rank {
	static const size_t value = A::dimensionality;
    };

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
	//BOOST_ASSERT((0 <= i && i < int(a.shape()[A::dimensionality-1])));
	return a.data()[i];
    }

    template<class A, class K>
    BOOST_GPU_ENABLED
    typename reference<A>::type
    access(A &a, int i, K) {
	BOOST_STATIC_ASSERT((K::value == int(A::dimensionality)));
	//BOOST_ASSERT((0 <= i && i < int(a.shape()[A::dimensionality-1])));
	typedef typename reference<A>::type R;
	return R(a.data() + i*a.stride(), a.shape()+1);
    }

}
}

namespace array {

    template<typename T, size_t N, class Tag, typename P>
    struct const_adapter {

	typedef T element;
	static const size_t dimensionality = N;
	typedef size_t size_type;

	typedef typename boost::mpl::if_c<
	    (N == 1), const T&, const_adapter<T,N-1>
	    >::type const_reference;
	typedef const_reference reference;

	typedef const T* const_iterator;

	BOOST_GPU_ENABLED
	const_adapter() {
	    size_t dims[N] = { 0 };
	    initialize(0, dims);
	}

	template<typename U>
	BOOST_GPU_ENABLED
	const_adapter(P data, const U *dims) {
	    initialize(data, dims);
	}

	BOOST_GPU_ENABLED
	P data() { return data_; }
	BOOST_GPU_ENABLED
	const T* data() const { return data_; }

	BOOST_GPU_ENABLED
	P origin() { return data_; }
	BOOST_GPU_ENABLED
	const T* origin() const { return data_; }

	BOOST_GPU_ENABLED
	size_type stride() const {
	    return stride_;
	}

	BOOST_GPU_ENABLED
	const size_type (&shape() const)[N] {
	    return shape_;
	}

	size_type num_elements() const {
	    return num_elements_;
	}

	BOOST_GPU_ENABLED
	const_reference operator[](int i) const {
	    return detail::access(*this, i, boost::mpl::int_<N>());
	}

	BOOST_GPU_ENABLED
	const_iterator begin() const {
	    BOOST_STATIC_ASSERT((N == 1));
	    return this->data_;
	}
	BOOST_GPU_ENABLED
	const_iterator end() const {
	    return begin() + this->shape()[0];
	}

	void swap(const_adapter &other) {
	    swap<const_adapter>(other);
	}

    protected:
	P data_;
	size_type stride_, num_elements_;
	size_type shape_[N];

	template<class A>
	void swap(A &other) {
	    std::swap(data_, other.data_);
	    std::swap(stride_, other.stride_);
	    std::swap(num_elements_, other.num_elements_);
	    for (size_t i = 0; i < N; ++i) {
		std::swap(shape_[i], other.shape_[i]);
	    }
	}

    private:
	template<typename U>
	BOOST_GPU_ENABLED
	void initialize(P data, const U *dims) {
	    stride_ = 1;
	    num_elements_ = 1;
	    for (size_t i = 1; i < N; ++i) {
		stride_ *= dims[N-i];
	    }
	    for (size_t i = 0; i < N; ++i) {
		shape_[i] = dims[i];
		num_elements_ *= shape_[i];
	    }
	    data_ = data;
	}
    };


    template<typename T, size_t N, class Tag>
    struct adapter : const_adapter<T, N, Tag, T*> {
	typedef const_adapter<T, N, Tag, T*> base;

	typedef typename boost::mpl::if_c<
	    (N == 1), T&, adapter<T,N-1>
	    >::type reference;
	typedef T* iterator;

	BOOST_GPU_ENABLED
	adapter() : base() {}

	template<typename U>
	BOOST_GPU_ENABLED
	adapter(T* data, const U *dims)
	    : base(data, dims) {}

	template<typename U, class O>
	BOOST_GPU_ENABLED
	adapter& operator=(const adapter<U,N,O> &a) {
	    assign(*this, a);
	}

	using base::operator[];

	BOOST_GPU_ENABLED
	reference operator[](int i) {
	    return detail::access(*this, i, boost::mpl::int_<N>());
	}

	BOOST_GPU_ENABLED
	iterator begin() {
	    BOOST_STATIC_ASSERT((N == 1));
	    return this->data_;
	}
	BOOST_GPU_ENABLED
	iterator end() {
	    return begin() + this->shape()[0];
	}

	// static void clear(iterator begin, iterator end) {
	//     std::fill(begin, end, 0);
	// }

	void swap(adapter &other) {
	    base::template swap<adapter>(other);
	}
    };


}

#endif // ARRAY_ADAPTER_HPP
