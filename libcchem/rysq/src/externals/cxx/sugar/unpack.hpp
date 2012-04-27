#ifndef _SUGAR_PP_UNPACK_HPP_
#define _SUGAR_PP_UNPACK_HPP_

#include <boost/preprocessor/seq/for_each_i.hpp>
#include <boost/preprocessor/tuple/to_seq.hpp>
#include <boost/preprocessor/comma_if.hpp>
#include "cxx/sugar/detail/narg.hpp"

#define SUGAR_PP_UNPACK_I(r, data, i, elem) BOOST_PP_COMMA_IF(i) elem = (data)[i]

#define UNPACK(u, p)							\
    BOOST_PP_SEQ_FOR_EACH_I(SUGAR_PP_UNPACK_I, p,			\
			    BOOST_PP_TUPLE_TO_SEQ(SUGAR_PP_NARG u, u))


#endif /* _SUGAR_PP_UNPACK_HPP_ */
