#ifndef _CXX_SUGAR_PP_DIMENSION_HPP_
#define _CXX_SUGAR_PP_DIMENSION_HPP_


#include <boost/preprocessor/cat.hpp>
#include <boost/preprocessor/tuple/elem.hpp>
#include "cxx/sugar/detail/narg.hpp"

#define DIMENSION_1_INDEX(N0,_0) (_0)

#define DIMENSION_1(A, dims, index)				\
    (A)[ DIMENSION_1_INDEX(BOOST_PP_TUPLE_ELEM(1,0,dims),	\
			   BOOST_PP_TUPLE_ELEM(1,0,index)) ]


#define DIMENSION_2_INDEX(N0,N1,_0,_1) (_0) + ((_1)*(N0))

#define DIMENSION_2(A, dims, index)				\
    (A)[ DIMENSION_2_INDEX(BOOST_PP_TUPLE_ELEM(2,0,dims),	\
			   BOOST_PP_TUPLE_ELEM(2,1,dims),	\
			   BOOST_PP_TUPLE_ELEM(2,0,index),	\
			   BOOST_PP_TUPLE_ELEM(2,1,index)) ]


#define DIMENSION_3_INDEX(N0,N1,N2,_0,_1,_2)	\
    (_0) + ((_1) + (_2)*(N1))*(N0)

#define DIMENSION_3(A, dims, index)				\
    (A)[ DIMENSION_3_INDEX(BOOST_PP_TUPLE_ELEM(3,0,dims),	\
			   BOOST_PP_TUPLE_ELEM(3,1,dims),	\
			   BOOST_PP_TUPLE_ELEM(3,2,dims),	\
			   BOOST_PP_TUPLE_ELEM(3,0,index),	\
			   BOOST_PP_TUPLE_ELEM(3,1,index),	\
			   BOOST_PP_TUPLE_ELEM(3,2,index)) ]




#define DIMENSION(A, dims, index)					\
    BOOST_PP_CAT(DIMENSION_, SUGAR_PP_NARG dims) (A, dims, index)
 

#endif /* _CXX_SUGAR_PP_DIMENSION_HPP_ */
