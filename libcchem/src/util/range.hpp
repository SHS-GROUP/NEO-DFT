#ifndef _UTIL_RANGE_HPP_
#define _UTIL_RANGE_HPP_

#include "util/unpack.hpp"
#include <boost/typeof/typeof.hpp>
#include <boost/preprocessor/cat.hpp>

#define PP_RANGE_2(it, begin, end) BOOST_TYPEOF((begin)) it = (begin); \
    it != (end); ++it

#define PP_RANGE_1(it, range)  PP_RANGE_2(it, (range).begin(), (range).end())

#define RANGE(it, ...)  BOOST_PP_CAT(PP_RANGE_, PP_NARG(__VA_ARGS__)) (it, __VA_ARGS__)

#define XRANGE(i, N) BOOST_TYPEOF((N)) i = 0; i < (N); ++i


#endif /* _UTIL_RANGE_HPP_ */
