#ifndef PHOENIX_BIND_NAME_HPP
#define PHOENIX_BIND_NAME_HPP

#include <boost/spirit/home/phoenix/core.hpp>
#include <boost/spirit/home/phoenix/function.hpp>
#include <boost/typeof/typeof.hpp>

#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum_params_with_a_default.hpp>

#define PHOENIX_BIND_NAME_DEFAULT_PARAMS(T)				\
    BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(PHOENIX_COMPOSITE_LIMIT, T, void)


namespace boost {
namespace phoenix {


    template<class F, PHOENIX_BIND_NAME_DEFAULT_PARAMS(class A)>
    struct function_result;

    template<class RT, class C>
    struct function_result<RT(C::*)()> { typedef RT type; };

    template<class RT, class C>
    struct function_result<RT(C::*)() const> { typedef RT type; };

#define PHOENIX_FUNCTION_RESULT(z, N, CQ)				\
    template<class RT, class C, BOOST_PP_ENUM_PARAMS(N, class A)>	\
    struct function_result<RT(C::*)(BOOST_PP_ENUM_PARAMS(N, A)) CQ> {	\
	typedef RT type;						\
    };
    BOOST_PP_REPEAT_FROM_TO(1, PHOENIX_COMPOSITE_LIMIT,
			    PHOENIX_FUNCTION_RESULT, /*const*/);
    BOOST_PP_REPEAT_FROM_TO(1, PHOENIX_COMPOSITE_LIMIT,
			    PHOENIX_FUNCTION_RESULT, const);
#undef PHOENIX_FUNCTION_RESULT


#define PHOENIX_BIND_NAME_OPERATOR(z, N, NAME)			\
    template<class C, BOOST_PP_ENUM_PARAMS(N, class A)>			\
    typename result<C>::type						\
    operator()(C &object,						\
	       BOOST_PP_ENUM_BINARY_PARAMS(N, A, const& a)) const {	\
 	return object.NAME(BOOST_PP_ENUM_PARAMS(N, a));		\
    }

#define PHOENIX_BIND_NAME(NAME)					\
    struct bind_ ## NAME ## _eval {						\
	template<class C, PHOENIX_BIND_NAME_DEFAULT_PARAMS(class A)>	\
	struct result {							\
	    typedef typename						\
	    function_result<BOOST_TYPEOF(&C::NAME)>::type type;	\
	};								\
	template<class C>						\
	typename result<C>::type					\
	operator()(C &object) const {					\
	    return object.NAME();					\
	}								\
	BOOST_PP_REPEAT_FROM_TO(1, PHOENIX_COMPOSITE_LIMIT,		\
				PHOENIX_BIND_NAME_OPERATOR, NAME)	\
    };									\
    const function<bind_ ## NAME ## _eval> NAME =			\
	function<bind_ ## NAME ## _eval>();


    PHOENIX_BIND_NAME(data);

    PHOENIX_BIND_NAME(start);
    PHOENIX_BIND_NAME(stop);

    PHOENIX_BIND_NAME(get);
    PHOENIX_BIND_NAME(put);


} // namespace phoenix
} // namespace boost

#endif // PHOENIX_BIND_NAME_HPP
