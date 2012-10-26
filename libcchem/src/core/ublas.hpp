#ifndef UBLAS_EXPRESSION_TRAITS_HPP
#define UBLAS_EXPRESSION_TRAITS_HPP

#include <boost/utility/enable_if.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/numeric/ublas/expression_types.hpp>
#include <boost/mpl/bool.hpp>

namespace boost { namespace numeric { namespace ublas {

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

} } }

#endif // UBLAS_EXPRESSION_TRAITS_HPP
