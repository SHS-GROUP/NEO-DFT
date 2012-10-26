#ifndef UTILITY_ARRAY
#define UTILITY_ARRAY

#include <boost/array.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/repeat.hpp>

namespace utility {
	
    template<typename T, size_t N>
    struct array : boost::array<T,N> {
	typedef boost::array<T,N> base_type;
	typedef T elements[N];
	template<typename U>
	array(const boost::array<U,N> & array)
	    : base_type(array) {}
	operator const elements&() const { return this->elems; }
	operator elements&() { return this->elems; }
	template<typename U>
	operator array<U,N>() const { return *this; }
    };

#define UTILITY_MAKE_ARRAY(z, N, text)					\
    template<typename T>						\
    array<T,(N)>							\
    make_array(BOOST_PP_ENUM_PARAMS(N, const T &t)) {			\
	boost::array<T,(N)> array = {{ BOOST_PP_ENUM_PARAMS(N, t) }};	\
	return array;							\
    }

    BOOST_PP_REPEAT(6, UTILITY_MAKE_ARRAY, nil)

#undef UTILITY_MAKE_ARRAY

}

#endif // UTILITY_ARRAY
