#ifndef PHOENIX_ARRAY_HPP
#define PHOENIX_ARRAY_HPP

#include <boost/spirit/include/phoenix_function.hpp>
#include <boost/typeof/typeof.hpp>
#include <boost/array.hpp>

#include <boost/preprocessor/repetition/enum.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_params_with_a_default.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>

namespace boost { namespace phoenix {

	template<typename T, size_t N>
	struct array_type : boost::array<T,N> {
	    typedef boost::array<T,N> base_type;
	    typedef T elements[N];
	    array_type(const boost::array<T,N> &array)
		: boost::array<T,N>(array) {}
	    operator const elements&() const { return this->elems; }
	};

	template<typename T>
	struct array_impl {
	    template<BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(PHOENIX_LIMIT,
							 typename A, void)>
	    struct result { typedef void type; };

#define PHOENIX_ARRAY_RESULT(z, N, text)				\
	    template<BOOST_PP_ENUM_PARAMS(N, typename A)>		\
	    struct result<BOOST_PP_ENUM_PARAMS(N, A)> {			\
		typedef array_type<T,(N)> type; };
	    BOOST_PP_REPEAT_FROM_TO(1, PHOENIX_LIMIT,
				    PHOENIX_ARRAY_RESULT, nil)
#undef PHOENIX_ARRAY_RESULT

#define PHOENIX_ARRAY_OPERATOR(z, N, text)				\
	    array_type<T,(N)>						\
	    operator()(BOOST_PP_ENUM_PARAMS(N, const T &t)) const {	\
		boost::array<T,(N)> array = {{				\
			BOOST_PP_ENUM_PARAMS(N, t) }};			\
		return array;						\
	    }
	    BOOST_PP_REPEAT(PHOENIX_LIMIT, PHOENIX_ARRAY_OPERATOR, nil)
#undef PHOENIX_ARRAY_OPERATOR
	};

	template<typename T>
	boost::phoenix::function<array_impl<T> > array() {
	    return array_impl<T>();
	}
	

    }				 

}

#endif // PHOENIX_ARRAY_HPP
