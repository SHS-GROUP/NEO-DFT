#ifndef TRANSFORM_LAMBDA_HPP
#define TRANSFORM_LAMBDA_HPP

#include "core/transform/transform.hpp"

#include <boost/spirit/home/phoenix/function/function.hpp>
#include <boost/preprocessor/repetition/repeat_from_to.hpp>
#include <boost/preprocessor/repetition/enum_params.hpp>
#include <boost/preprocessor/repetition/enum_binary_params.hpp>
#include <boost/preprocessor/repetition/enum_params_with_a_default.hpp>

namespace transform {

    namespace detail {

	template<class F>
	struct lambda {
	    lambda(F &transform) : transform_(transform) {}

	    template<
		BOOST_PP_ENUM_PARAMS_WITH_A_DEFAULT(PHOENIX_COMPOSITE_LIMIT,
						    class A, void)>
	    struct result { typedef void type; };

#define TRANSFORM_LAMBDA_OPERATOR(z, N, data)				\
	template<BOOST_PP_ENUM_PARAMS(N, class A)>			\
	void operator()(BOOST_PP_ENUM_BINARY_PARAMS(N, A, &a)) const { \
	    transform_(BOOST_PP_ENUM_PARAMS(N, a));			\
	}

	    BOOST_PP_REPEAT_FROM_TO(1, PHOENIX_COMPOSITE_LIMIT,
				    TRANSFORM_LAMBDA_OPERATOR, ());
	private:
	    F &transform_;
	};
    }


    template<bool Symmetry>
    struct Lambda<Transform<Symmetry> >
	: boost::phoenix::function<detail::lambda<Transform<Symmetry> > >
    {
	typedef Transform<Symmetry> F;
	typedef detail::lambda<F> lambda;
	Lambda(F &transform)
	    : boost::phoenix::function<lambda>(lambda(transform)) {}
    };
    
    template<bool Symmetry>
    Lambda<Transform<Symmetry> >
    Transform<Symmetry>::operator()() {
	return Lambda<Transform>(*this);
    }

}

#endif // TRANSFORM_LAMBDA_HPP
