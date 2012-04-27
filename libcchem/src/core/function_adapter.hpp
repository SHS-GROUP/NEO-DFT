#ifndef FUNCTION_ADAPTER_HPP
#define FUNCTION_ADAPTER_HPP

#include <boost/spirit/home/phoenix/function/function.hpp>

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/mpl/bool.hpp>

namespace boost { namespace phoenix {

namespace ublas {

template<class E, class Enable = void>
struct is_vector_expression : boost::mpl::false_ { };

template<class E, class Enable = void>
struct is_matrix_expression : boost::mpl::false_ { };

template<class E>
struct is_vector_expression<E, typename boost::enable_if<
				   boost::is_same<
				       typename E::type_category,
				       boost::numeric::ublas::vector_tag>
				   >::type
			    > : boost::mpl::true_{};

template<class E>
struct is_matrix_expression<E, typename boost::enable_if<
				   boost::is_same<
				       typename E::type_category,
				       boost::numeric::ublas::matrix_tag>
				   >::type
			    > : boost::mpl::true_{};

template<class T, class Enable = void>
struct function_adapter;

template<class E>
struct function_adapter<E, typename boost::enable_if<
			       is_matrix_expression<E>
			       >::type> {
    typedef typename E::reference reference;
    template<typename Arg> struct result { typedef reference type; };
    function_adapter(E &e): e_(e) {}
    template<typename T1, typename T2>
    reference& operator()(const T1 &_1, const T2 &_2) const {
	return e_(_1,_2);
    }
    E &e_;
};

template<class E>
struct function : phoenix::function<function_adapter<E> > {
    function(E &e) : phoenix::function<function_adapter<E> >(e) {}
};

} } }

#endif // FUNCTION_ADAPTER_HPP
