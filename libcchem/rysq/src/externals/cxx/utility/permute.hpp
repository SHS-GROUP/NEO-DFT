#ifndef _CXX_UTILITY_PERMUTE_HPP
#define _CXX_UTILITY_PERMUTE_HPP

#include <boost/config.hpp>
#include <boost/array.hpp>
#include <boost/iterator/iterator_traits.hpp>
#include "cxx/array.hpp"

namespace cxx {
namespace utility {

    namespace permutation_ {

	template<class Q>
	struct value_{
	    typedef typename Q::value_type type;
	};

	template< typename T>
	struct value_<T[4]> { typedef T type; };
	  
	BOOST_GPU_ENABLED
	static int index(int index, int mask) {
	    return (mask>>(index*2) & 0x03);
	}

    }

    template<class Q>
    struct permutation {
	typedef typename  permutation_::value_<Q>::type value_type;
	typedef const value_type& const_reference;
	BOOST_GPU_ENABLED
	permutation(const Q &quartet, int mask) : data_(quartet), mask_(mask) {}
	BOOST_GPU_ENABLED
	const_reference operator[](size_t i) const {
	    return data_[permutation_::index(i, mask_)];
	}
    private:
	const Q &data_;
	const int mask_;

	template<class A>
	BOOST_GPU_ENABLED
	static
	const_reference element(const A &data, size_t i) {
	    return data[i];
	}

	template<typename T>
	BOOST_GPU_ENABLED
	static
	const_reference element(const boost::array<T,4> &data, size_t i) {
	    return data.elems[i];
	}

    };


    template<class T>
    struct permutation< boost::array<T,4> > {
	typedef boost::array<T,4> Q;
	typedef typename  permutation_::value_<Q>::type value_type;
	BOOST_GPU_ENABLED
	permutation(const Q &quartet, int mask) : data_(quartet), mask_(mask) {}
	BOOST_GPU_ENABLED
	const value_type& operator[](size_t i) const {
	    return data_.elems[permutation_::index(i, mask_)];
	}
    private:
	const Q &data_;
	const int mask_;
    };

    template<class Q>
    BOOST_GPU_ENABLED
    permutation<Q> make_permutation(const Q &quartet, int mask) {
	return permutation<Q>(quartet, mask);
    }

    template <typename T>
    BOOST_GPU_ENABLED
    void permute(T &i, T &j, T &k, T &l, int mask) {
	T q[4] = {i, j, k, l};
	permutation<T[4]> p(q, mask);
	i = p[0];
	j = p[1];
	k = p[2];
	l = p[3];
    }

    template<typename T>
    BOOST_GPU_ENABLED
    T permute(int i, const T *array, int mask) {
	return array[permutation_::index(i, mask)];
    }

    template <typename T>
    BOOST_GPU_ENABLED
    boost::array<T,4> permute(const boost::array<T,4> &array, int mask) {
	permutation<boost::array<T,4> > p(array, mask);
	boost::array<T,4> q = {{ p[0], p[1], p[2], p[3] }};
	// std::cout << mask << q << array << std::endl;
	return q;
    }

    template< class Q>
    BOOST_GPU_ENABLED
    Q permute(const Q &array, int mask) {
	Q q;
	permutation<Q> p(array, mask);
	for (int i = 0; i < 4; ++i) {
	    q[i] = p[i];
	}
	return q;
    }



}
}

#endif /* _CXX_UTILITY_PERMUTE_HPP */
