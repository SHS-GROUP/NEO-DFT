#ifndef _SUGAR_PP_XRANGE_HPP_
#define _SUGAR_PP_XRANGE_HPP_

#include <boost/typeof/typeof.hpp>

#define XRANGE(i, N) BOOST_TYPEOF((N)) i = 0; i < (N); ++i

#endif /* _SUGAR_PP_XRANGE_HPP_ */
