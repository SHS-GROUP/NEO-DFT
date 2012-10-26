#ifndef TENSOR_PRODUCT_HPP
#define TENSOR_PRODUCT_HPP

#include "tensor/forward.hpp"
#include "tensor/expression.hpp"
#include "tensor/lambda.hpp"
#include "tensor/index.hpp"
#include "tensor/view.hpp"
#include "tensor/functional.hpp"
#include "tensor/assert.hpp"

#include <boost/mpl/contains.hpp>
#include <boost/fusion/include/map.hpp>
#include <boost/fusion/include/transformation.hpp>
#include <boost/fusion/include/mpl.hpp>


namespace tensor {
namespace detail {


    template<class A, class B>
    struct index_product {
	typedef typename fusion::result_of::remove_if<
	    A, boost::mpl::contains<B, boost::mpl::_> >::type A_;
	typedef typename fusion::result_of::remove_if<
	    B, boost::mpl::contains<A, boost::mpl::_> >::type B_;
	typedef typename fusion::result_of::as_set<
	    typename fusion::result_of::join<A_, B_>::type
	    >::type keys_type;

	template<class K>
	struct key_as_pair {
	    typedef index<K::value> type;
	};

	typedef typename fusion::result_of::as_map<
	    typename boost::mpl::transform<
		keys_type, key_as_pair<boost::mpl::_> >::type
	    >::type type;
	// typedef typename boost::mpl::print< type>::type _;
	static const size_t rank = boost::mpl::size<type>::value;
    };


    template<class A, class B>
    struct product_expression_base {
	typedef detail::index_product<
	    typename A::keys_type, typename B::keys_type> index_product;
	typedef typename index_product::keys_type keys_type;
	typedef typename index_product::type indices_type;
	static const size_t rank = index_product::rank;

	typedef typename A::value_type value_type;
	typedef value_type const_result_type;

	product_expression_base(const expression<A> &a, const expression<B> &b)
	    : a_(a()), b_(b()) {}

	const indices_type& indices() const;

	template<class S>
	const_result_type operator[](const S &index) const {
	    value_type v(0);
	    using namespace lambda;
	    contract((ref(v) += _1*_2), a_, b_, index);
	    return v;
	}
    private:
	A a_;
	B b_;
    };

}
}

namespace tensor {


    template<class A, class B, 
	     size_t N = detail::product_expression_base<A,B>::rank>
    struct product_expression :
	detail::product_expression_base<A,B>,
	expression<product_expression<A,B> >
    {
	typedef detail::product_expression_base<A,B> base_type;
	product_expression(const expression<A> &a, const expression<B> &b)
	    : base_type(a,b) {}
	using base_type::operator[];
    };

    template<class A, class B>
    struct product_expression<A,B,0> :
	detail::product_expression_base<A,B>,
	expression<product_expression<A,B> >
    {
	typedef detail::product_expression_base<A,B> base_type;
	typedef typename base_type::value_type value_type;
	product_expression(const expression<A> &a, const expression<B> &b)
	    : base_type(a,b) {}
	using base_type::operator[];
	operator value_type() const;
    };


    template<class A, class B>
    product_expression<A,B>
    operator*(const expression<A> &a, const expression<B> &b) {
	return product_expression<A,B>(a,b);
    }

} // namespace tensor


#endif // TENSOR_PRODUCT_HPP

