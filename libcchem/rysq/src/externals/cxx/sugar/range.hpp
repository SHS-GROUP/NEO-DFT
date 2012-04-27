#ifndef _SUGAR_PP_RANGE_HPP_
#define _SUGAR_PP_RANGE_HPP_

//#include "util/unpack.hpp"
#include <boost/typeof/typeof.hpp>
#include <boost/preprocessor/cat.hpp>
#include "cxx/sugar/detail/narg.hpp"

#define SUGAR_PP_RANGE_2(it, begin, end) BOOST_TYPEOF((begin)) it = (begin); \
  it != (end); ++it

#define SUGAR_PP_RANGE_1(it, range) SUGAR_PP_RANGE_2(it, (range).begin(), (range).end())

#define RANGE(it, ...)  BOOST_PP_CAT(SUGAR_PP_RANGE_, SUGAR_PP_NARG(__VA_ARGS__)) \
  (it, __VA_ARGS__)



#endif /* _SUGAR_PP_RANGE_HPP_ */
