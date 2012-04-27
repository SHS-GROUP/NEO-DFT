#ifndef _CXX_NAMESPACE_HPP_
#define _CXX_NAMESPACE_HPP_

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/tuple/to_seq.hpp>
#include <boost/preprocessor/comma_if.hpp>
#include "cxx/sugar/detail/narg.hpp"

#define BEGIN_NAMESPACE_1(r, data, elem) namespace elem {
#define END_NAMESPACE_1(r, data, elem) } /* namespace elem */

#define BEGIN_NAMESPACE(...)						\
    BOOST_PP_SEQ_FOR_EACH(BEGIN_NAMESPACE_1, (nil),			\
			  BOOST_PP_TUPLE_TO_SEQ(SUGAR_PP_NARG(__VA_ARGS__), \
						  (__VA_ARGS__)))

#define END_NAMESPACE(...)						\
    BOOST_PP_SEQ_FOR_EACH(END_NAMESPACE_1, (nil),			\
			  BOOST_PP_TUPLE_TO_SEQ(SUGAR_PP_NARG(__VA_ARGS__), \
						  (__VA_ARGS__)))

#endif /* _CXX_NAMESPACE_HPP_ */
