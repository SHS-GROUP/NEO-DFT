#ifndef TENSOR_ORDER_HPP
#define TENSOR_ORDER_HPP

#include <boost/tuple/tuple.hpp>

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/at.hpp>
#include <boost/mpl/copy.hpp>
#include <boost/mpl/iterator_range.hpp>
#include <boost/mpl/find_if.hpp>
#include <boost/mpl/count_if.hpp>
#include <boost/mpl/comparison.hpp>
#include <boost/mpl/back_inserter.hpp>
#include <boost/mpl/assert.hpp>

#include <boost/preprocessor/inc.hpp>
#include <boost/preprocessor/enum.hpp>
#include <boost/preprocessor/enum_params.hpp>
#include <boost/preprocessor/enum_params_with_a_default.hpp>


namespace tensor {
namespace detail {

    namespace mpl = boost::mpl;

#define TENSOR_LIMIT 4

    template<BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(TENSOR_LIMIT, int O, -1)>
    struct order {
	typedef mpl::vector_c<int, BOOST_PP_ENUM_PARAMS(TENSOR_LIMIT, O)> O;
	typedef mpl::_1 _1;
	typedef mpl::_2 _2;
	typedef mpl::int_<0> zero;
	typedef mpl::less<_1, zero> negative;
	typedef typename mpl::find_if<O, negative>::type last;
	typedef typename mpl::iterator_range<
	    mpl::begin<O>, last>::type range;
	typedef typename mpl::copy<
	    range, mpl::back_inserter<mpl::vector<> >
	    >::type vector;
	static const size_t rank = mpl::size<vector>::value;
    private:
	static const int index_sum = mpl::fold<
	vector, zero, mpl::plus<_1,_2> >::type::value;
	BOOST_MPL_ASSERT_RELATION(index_sum, ==, (rank*rank - rank)/2);

	typedef typename mpl::count_if<
	    typename mpl::iterator_range<last, mpl::end<O> >::type,
	    mpl::greater_equal<_1, zero>
	    >::type must_be_zero;
	BOOST_MPL_ASSERT_NOT((must_be_zero));

    public:
	template<BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(TENSOR_LIMIT, class A,
						     boost::tuples::null_type)>
	struct result {
	    typedef boost::tuple<
		BOOST_PP_ENUM_BINARY_PARAMS(TENSOR_LIMIT, const A,
					    & BOOST_PP_INTERCEPT)> data;
	    template<int I>
	    struct element {
		static const int index = boost::mpl::at_c<vector,I>::type::value;
		typedef typename boost::tuples::element<index, data>::type type;
	    };
#define TENSOR_ORDER_TUPLE_TYPE(z, index, data) typename element<index>::type
	    typedef boost::tuple<BOOST_PP_ENUM(TENSOR_LIMIT,
					       TENSOR_ORDER_TUPLE_TYPE, nil)> type;
#undef TENSOR_ORDER_TUPLE_TYPE
	    static type tuple(const data &A) {
#define TENSOR_ORDER_TUPLE_GET(z, index, data) boost::get<index>(A)
		return type(BOOST_PP_ENUM(TENSOR_LIMIT,
					  TENSOR_ORDER_TUPLE_GET, nil));
#undef TENSOR_ORDER_TUPLE_GET
	    }
	};
    };

} // namespace detail


    template<BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(BOOST_PP_INC(TENSOR_LIMIT),
						 int O, -1)>
    struct Order;

#define TENSOR_ORDER(N)							\
    template<BOOST_PP_ENUM_PARAMS(N, int O)>				\
    struct Order<BOOST_PP_ENUM_PARAMS(N, O)>				\
	: detail::order<BOOST_PP_ENUM_PARAMS(N, O)>			\
    {									\
	typedef detail::order<BOOST_PP_ENUM_PARAMS(N, O)> base;		\
	template<BOOST_PP_ENUM_PARAMS(N, class A)>			\
	    static typename base::template result<			\
	BOOST_PP_ENUM_PARAMS(N,A)>::type				\
	    tuple(BOOST_PP_ENUM_BINARY_PARAMS(N, const A, &a)) {	\
	    typedef typename base::template result<			\
	    BOOST_PP_ENUM_PARAMS(N,A)>::tuple R;			\
	    return R(boost::tie(BOOST_PP_ENUM_PARAMS(N, a)));		\
	}								\
    }

    TENSOR_ORDER(1);
    TENSOR_ORDER(2);
    TENSOR_ORDER(3);
    TENSOR_ORDER(4);

#undef TENSOR_ORDER


} // namespace tensor


    


#endif // TENSOR_ORDER_HPP
