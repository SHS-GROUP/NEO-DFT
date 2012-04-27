#ifndef _UTIL_ASSERT_HPP_
#define _UTIL_ASSERT_HPP_

#include <cassert>

#define ASSERT_RANGE(i, lower, upper) assert(int(lower) <= (i) && (i) < int(upper))

#endif /* _UTIL_ASSERT_HPP_ */
