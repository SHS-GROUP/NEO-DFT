#ifndef TENSOR_DETAIL_CONTRACT_HPP
#define TENSOR_DETAIL_CONTRACT_HPP

#include "tensor/forward.hpp"

#include <cstdlib>
#include <boost/type_traits.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/mpl/assert.hpp>

#include <boost/mpl/equal.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/fusion/include/mpl.hpp>
#include <boost/fusion/include/intrinsic.hpp>

#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/bindings/blas.hpp>
#include <boost/numeric/bindings/ublas.hpp>

namespace tensor {
namespace detail {

    template<class S, class I>
    struct is_first : boost::mpl::bool_<
	(boost::mpl::front<S>::type::value == I::value)> {};

    template<class S, class I>
    struct is_last : boost::mpl::bool_<
	(boost::mpl::back<S>::type::value == I::value)> {};

#define TENSOR_ASSERT_INDEX_IS_FIRST_OR_LAST(A,I)				\
    BOOST_MPL_ASSERT_MSG((is_first<typename A::indices_type, I>::value ||	\
			  is_last<typename A::indices_type, I>::value),		\
			 TENSOR_INDEX_IS_NOT_FIRST_OR_LAST, (A,I))

    template<class A>
    struct blas_view {

	typedef typename A::value_type value_type;
	typedef boost::numeric::ublas::array_adaptor<
	    typename A::value_type> array_type;
	typedef boost::numeric::ublas::matrix<
	    value_type, boost::numeric::ublas::column_major, array_type>
	matrix_type;
	typedef typename boost::numeric::ublas::matrix_range<
	    matrix_type> matrix_range;
	typedef boost::numeric::ublas::range range;

	template<size_t L>
	blas_view(A &a, boost::mpl::size_t<L>) {
	    initialize<L,A::rank-L>(a.data());
	}

	struct trans {
	    static blas_view::matrix_range matrix_range();
	    typedef BOOST_TYPEOF(boost::numeric::ublas::trans(matrix_range())) type;
	};

	matrix_range operator()(boost::mpl::false_ = boost::mpl::false_()) {
	    return boost::numeric::ublas::project(matrix_, r1_, r2_);
	}

	typename trans::type operator()(boost::mpl::true_) {
	    return boost::numeric::ublas::trans(this->operator()());
	}

    private:
	matrix_type matrix_;
	range r1_, r2_;

	void resize(size_t size1, size_t size2, typename A::array_type &a) {
	    array_type data(size1*size2, storage::data(a));
	    matrix_.data().swap(data);
	    matrix_.resize(size1, size2, false);
	}

	template<size_t M, size_t N>
	void initialize(typename A::array_type &a) {
	    bool status = (storage::stride(a,0) == 1);
	    int dim[M+N], size[M+N];
	    for (size_t i = 0; i < (M+N-1); ++i) {
		div_t v = div(storage::stride(a,i+1), storage::stride(a,i));
		//std::cout << v.quot << std::endl;
		dim[i] = v.quot;
		size[i] = storage::size(a,i);
		status *= (v.rem == 0);
	    }
	    {
		bool valid_matrix_layout = status;
		BOOST_ASSERT(valid_matrix_layout);
		dim[M+N-1] = storage::size(a,M+N-1);
		size[M+N-1] = storage::size(a,M+N-1);
	    }

	    size_t msize[2];
	    for (int k = 0, j = 0, n = M; k < 2; ++k) {
		msize[k] = size[j];
		// std::cout << j << "+" << size[j] << " " << dim[j] << std::endl;
		for (int i = j+1; i < n; ++i) {
		    status *= (dim[i-1] == size[i-1]);
		    msize[k] *= size[i];
		    // std::cout << i << "-" << size[i] << " " << dim[i] << std::endl;
		    //std::cout << i << " " << status << std::endl;
		}
		bool valid_matrix_layout = status;
		BOOST_ASSERT(valid_matrix_layout);
		j += M;
		n += N;
	    }

	    resize(storage::stride(a,M), msize[1], a);
	    r1_ = range(0, msize[0]);
	    r2_ = range(0, msize[1]);
	}

    };

    template<typename Alpha, class A, class B, typename Beta, class C>
    void contract_impl(Alpha alpha, A a, B b, Beta beta, C c) {

	namespace fusion = boost::fusion;
	namespace mpl = boost::mpl;

	typedef typename fusion::result_of::intersection<
	typename A::keys_type, typename B::keys_type>::type keys_type;

	BOOST_MPL_ASSERT_MSG((mpl::size<keys_type>::value == 1),
			     TENSOR_CONTRACT_MUST_HAVE_ONE_CONTRACTED_INDEX,
			     (keys_type));
	typedef typename mpl::front<keys_type>::type K;

	TENSOR_ASSERT_INDEX_IS_FIRST_OR_LAST(A,K);
	TENSOR_ASSERT_INDEX_IS_FIRST_OR_LAST(B,K);

	typedef typename mpl::erase_key<typename A::indices_type, K>::type IA;
	typedef typename mpl::erase_key<typename B::indices_type, K>::type IB;
	typedef typename fusion::result_of::join<IA,IB>::type IC;

	BOOST_MPL_ASSERT_MSG((mpl::equal<IC, typename C::indices_type>::value),
			     TENSOR_INVALID_CONTRACTED_INDICES,
			     (IC, typename C::indices_type));
    
	typedef is_first<typename A::indices_type, K> transa;
	typedef is_last<typename B::indices_type, K> transb;

	typename mpl::size_t<(transa::value ? 1 : A::rank-1)> LA;
	typename mpl::size_t<(!transb::value ? 1 : B::rank-1)> LB;
	typename mpl::size_t<(A::rank-1)> LC;	

	blas_view<A> a_(a, LA);
	blas_view<B> b_(b, LB);
	blas_view<C> c_(c, LC);

	{
	    BOOST_AUTO(a, a_(transa()));
	    BOOST_AUTO(b, b_(transb()));
	    BOOST_AUTO(c, c_());

	    BOOST_ASSERT(a.size1() == c.size1());
	    BOOST_ASSERT(b.size2() == c.size2());
	    BOOST_ASSERT(a.size2() == b.size1());

	    namespace blas = boost::numeric::bindings::blas;
	    blas::gemm(alpha, a, b, beta, c);
	}
    }

}
}

#endif // TENSOR_DETAIL_CONTRACT_HPP
