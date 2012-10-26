#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/seq/for_each.hpp>

#ifndef TYPE__
#error "type is not defined"
#endif

namespace BOOST_PP_CAT(TYPE__ , __) {

#include "roots/roots0.hpp"
#include "roots/roots1.hpp"
#include "roots/roots2.hpp"
#include "roots/roots3.hpp"
#include "roots/roots4.hpp"
#include "roots/roots5.hpp"

}

template<>
struct roots<0,TYPE__> {
    BOOST_GPU_ENABLED
    static void evaluate(TYPE__ X, TYPE__ (&t2)[0], TYPE__ (&W)[1]) {
	BOOST_PP_CAT(TYPE__ , __)::roots0::evaluate(X, t2, W);
    }
    template<class T>
    BOOST_GPU_ENABLED
    static void evaluate(TYPE__ X, TYPE__ *t2, TYPE__ *W,
			 const T &thread) {
	BOOST_PP_CAT(TYPE__ , __)::roots0::evaluate(X, t2, W, thread);
    }
};

#define RYSQ_ROOTS_N(R, DATA, N)					\
    template<>								\
    struct roots<N,TYPE__> {						\
    BOOST_GPU_ENABLED							\
    static void evaluate(TYPE__ X, TYPE__ (&t2)[N], TYPE__ (&W)[N]) {	\
	BOOST_PP_CAT(TYPE__ , __)::					\
	    BOOST_PP_CAT(roots, N)::evaluate(X, t2, W);			\
    }									\
    template<class T>							\
    BOOST_GPU_ENABLED							\
    static void evaluate(TYPE__ X, TYPE__ *t2, TYPE__ *W,		\
			 const T &thread) {				\
    BOOST_PP_CAT(TYPE__ , __)::						\
	BOOST_PP_CAT(roots, N)::evaluate(X, t2, W, thread);		\
    }									\
};

BOOST_PP_SEQ_FOR_EACH(RYSQ_ROOTS_N, (), (1)(2)(3)(4)(5))

// template<>
// struct roots<2,TYPE__> {
//     template<class T>
//     BOOST_GPU_ENABLED
//     static void evaluate(TYPE__ X, TYPE__ (&t2)[2], TYPE__ (&W)[2],
// 			 const T &thread) {
// 	BOOST_PP_CAT(TYPE__ , __)::roots2::evaluate(X, t2, W, thread);
//     }
// };

// template<>
// struct roots<3,TYPE__> {
//     template<class T>
//     BOOST_GPU_ENABLED
//     static void evaluate(TYPE__ X, TYPE__ (&t2)[3], TYPE__ (&W)[3],
// 			 const T &thread) {
// 	BOOST_PP_CAT(TYPE__ , __)::roots3::evaluate(X, t2, W, thread);
//     }
// };

// template<>
// struct roots<4,TYPE__> {
//     template<class T>
//     BOOST_GPU_ENABLED
//     static void evaluate(TYPE__ X, TYPE__ (&t2)[4], TYPE__ (&W)[4],
// 			 const T &thread) {
// 	BOOST_PP_CAT(TYPE__ , __)::roots4::evaluate(X, t2, W, thread);
//     }
// };

// template<>
// struct roots<5,TYPE__> {
//     template<class T>
//     BOOST_GPU_ENABLED
//     static void evaluate(TYPE__ X, TYPE__ (&t2)[5], TYPE__ (&W)[5],
// 			 const T &thread) {
// 	BOOST_PP_CAT(TYPE__ , __)::roots5::evaluate(X, t2, W, thread);
//     }
// };
