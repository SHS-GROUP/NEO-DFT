#ifndef CC_TENSOR_HPP
#define CC_TENSOR_HPP

#include "array/permute.hpp"

#include <boost/multi_array.hpp>
#include <boost/numeric/ublas/adaptor.hpp>
#include <boost/mpl/vector_c.hpp>

namespace cc {


    template<size_t N>
    struct Tensor : boost::multi_array<double,N> {
	typedef boost::multi_array<double,N> Array;
	typedef typename boost::multi_array<
	    double,N+1>::reference reference;
	typedef typename boost::multi_array<
	    double,N+1>::const_reference const_reference;
	Tensor() {}
	template<class E>
	Tensor(const E &extents) : Array(extents) {}
	template<class E>
	void resize(const E &extents) {
	    Array::resize(extents);
	}
    };
    

    template<class T, size_t M, size_t N, class R>
    T as_matrix_impl(R a) {
	size_t size1 = 1;
	size_t size2 = 1;
	BOOST_AUTO(shape, a.shape()+M+N);
	for (size_t i = 0; i < M; ++i) {
	    size1 *= *(--shape);
	}
	for (size_t i = M; i < M+N; ++i) {
	    size2 *= *(--shape);
	}
	boost::numeric::ublas::column_major O;
	return boost::numeric::ublas::make_matrix(size1, size2, a.origin(), O);
    }

    template<size_t M, size_t N>
    typename boost::numeric::ublas::matrix_adaptor<
	double, boost::numeric::ublas::column_major>::type
    as_matrix(typename Tensor<M+N>::reference a) {
	typedef typename boost::numeric::ublas::matrix_adaptor<
	double, boost::numeric::ublas::column_major>::type T;
	return as_matrix_impl<T,M,N>(a);
    }

    template<size_t M, size_t N>
    typename boost::numeric::ublas::matrix_adaptor<
	const double, boost::numeric::ublas::column_major>::type
    as_matrix(typename Tensor<M+N>::const_reference a) {
	typedef typename boost::numeric::ublas::matrix_adaptor<
	double, boost::numeric::ublas::column_major>::type T;
	return as_matrix_impl<T,M,N>(a);
    }

}


namespace cc {
namespace detail {

    template<int I, int J, int K>
    void permute(Tensor<3>::reference a,
		 boost::mpl::int_<I>,
		 boost::mpl::int_<J>,
		 boost::mpl::int_<K>) {
	array::permute<I,J,K>(a);
    }

    template<int I, int J, int K>
    void permute(Tensor<4> &a,
		 boost::mpl::int_<I>,
		 boost::mpl::int_<J>,
		 boost::mpl::int_<K>,
		 boost::mpl::int_<3>) {
	using boost::mpl::int_;
	for (size_t i = 0; i < a.shape()[0]; ++i) {
	    permute(a[i], int_<I>(), int_<J>(), int_<K>());
	}
    }

} // namespace detail
} // namespace cc

namespace cc {

    template<int I, int J, int K>
    void permute(Tensor<3>::reference a) {
	using boost::mpl::int_;
	detail::permute(a, int_<I>(), int_<J>(), int_<K>());
    }

    template<int I, int J, int K, int L>
    void permute(Tensor<4> &a) {
	using boost::mpl::int_;
	detail::permute(a, int_<I>(), int_<J>(), int_<K>(), int_<L>());
    }


} // namespace cc

#endif /* CC_TENSOR_HPP */
